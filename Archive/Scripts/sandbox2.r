

library(dplyr)
library(tidyr)
library(ggplot2)
library(here)
library(rstan)
library(biscale)

dataset <- "Full"

load(here("Rout", paste0("fit_", dataset, ".RData")))



region_id <- stan_data$region_id

site_levels
region_levels

eval_seasonality_curves <- function(
  fit,
  H,
  P,
  S,
  entity = c("province", "region"),
  region_id = NULL,
  region_labels = NULL,
  month_grid = seq(1, 12, length.out = 241),
  probs = c(0.025, 0.5, 0.975)
) {
  entity <- match.arg(entity)

  dr <- rstan::extract(fit, pars = c("a_cos", "a_sin"))
  a_cos <- dr$a_cos  # iter x H x P x S
  a_sin <- dr$a_sin  # iter x H x P x S

  stopifnot(length(dim(a_cos)) == 4, dim(a_sin)[1] == dim(a_cos)[1])

  n_iter <- dim(a_cos)[1]
  theta <- 2 * pi * month_grid / 12

  # Precompute basis: (grid x H)
  cos_basis <- sapply(seq_len(H), \(h) cos(h * theta))
  sin_basis <- sapply(seq_len(H), \(h) sin(h * theta))

  # Result: iter x grid x P x S
  f <- array(0, dim = c(n_iter, length(month_grid), P, S))
for (p in seq_len(P)) {
  for (s in seq_len(S)) {
    A_cos <- a_cos[, , p, s, drop = FALSE][, , 1, 1]  # n_iter x H
    A_sin <- a_sin[, , p, s, drop = FALSE][, , 1, 1]  # n_iter x H

    f[, , p, s] <- A_cos %*% t(cos_basis) + A_sin %*% t(sin_basis)
  }
}

  mult <- exp(f) # iter x grid x P x S

  # Optionally aggregate provinces -> regions for each draw
  if (entity == "region") {
    stopifnot(!is.null(region_id))
    R <- max(region_id)
    mult_reg <- array(NA_real_, dim = c(n_iter, length(month_grid), P, R))
    for (r in seq_len(R)) {
      s_idx <- which(region_id == r)
      mult_reg[, , , r] <- apply(mult[, , , s_idx, drop = FALSE], c(1, 2, 3), mean)
    }
    mult <- mult_reg
    S <- dim(mult)[4]

    if (is.null(region_labels)) region_labels <- paste0("Region ", seq_len(S))
    entity_names <- region_labels
  } else {
    entity_names <- paste0("Province ", seq_len(S))
  }

  # Summarize across iterations
  out <- tidyr::crossing(
    x = month_grid,
    period = factor(seq_len(P), labels = paste0("Period ", seq_len(P))),
    entity = factor(entity_names, levels = entity_names)
  ) |>
    mutate(q025 = NA_real_, q500 = NA_real_, q975 = NA_real_)

  # Create index matching
  out <- out |>
    mutate(
      idx = match(
        paste(x, period, entity, sep = "|"),
        paste(month_grid[match(x, month_grid)], period, entity, sep = "|")
      )
    )
  
  for (e in seq_len(S)) {
    for (p in seq_len(P)) {
      idx_mask <- out$period == levels(out$period)[p] & out$entity == entity_names[e]
      qs_all <- apply(mult[, , p, e], 2, \(col) stats::quantile(col, probs = probs, names = FALSE))
      out$q025[idx_mask] <- qs_all[1, ]
      out$q500[idx_mask] <- qs_all[2, ]
      out$q975[idx_mask] <- qs_all[3, ]
    }
  }

  out
}

H <- stan_data$H
P <- stan_data$P
S <- length(stan_data$region_id)

season_prov_smooth <- eval_seasonality_curves(
  fit = fit, H = H, P = P, S = S,
  entity = "province",
  region_id = region_id,
  month_grid = seq(1, 12, length.out = 241)
)


seasonality_indices <- season_prov_smooth |>
  rename(region = entity, month = x, mult = q500) |>
  group_by(region, period) |>
  summarize(
    peak_month   = month[which.max(mult)], # Month of peak seasonality
    trough_month = month[which.min(mult)], # Month of trough seasonality
    peak_mult    = max(mult), # Peak seasonal multiplier
    trough_mult  = min(mult), # Trough seasonal multiplier
    amp_ratio    = peak_mult / trough_mult, # Amplitude ratio (peak/trough)
    amp_log      = log(amp_ratio), # Log amplitude (log(peak/trough))
    .groups = "drop"
  )



wss <- sapply(1:10, function(k) {
  kmeans(seasonality_indices |> 
  select(-region, -period, -clusters, -trough_month, -trough_mult, -amp_ratio), centers = k, nstart = 25)$tot.withinss
})

plot(1:10, wss, type = "b", pch = 19,
     xlab = "Number of clusters (k)",
     ylab = "Total within-cluster sum of squares",
     main = "Elbow Plot for K-means Clustering")

library(factoextra)
k <- 3 # Start with 3 archetypes
seasonality_indices$clusters <- kmeans(seasonality_indices |> 
  select(-region, -period, -clusters, -trough_month, -trough_mult, -amp_ratio),centers = k, nstart = 25)$cluster


# Summarize median peak month and amplitude ratio per cluster
cluster_summaries <- seasonality_indices |>
  group_by(clusters) |>
  summarize(
    median_peak_month = median(peak_month),
    median_amp_ratio = median(amp_ratio),
    .groups = "drop"
  ) |>
  mutate(
    label = paste0(
      "Peak month: ", formatC(median_peak_month, format = "f", digits = 1), "\n",
      "Amplitude-to-through ratio: ", formatC(median_amp_ratio, format = "f", digits = 1)
    ),
    x = 1,   # adjust x/y for label placement as needed
    y = 2.0
  )



season_prov_smooth  %>%
  left_join(seasonality_indices %>% select(region, period, clusters), by = c("entity" = "region", "period")) %>%
  group_by(x, clusters) %>%
  summarize(q_med = median(q500),
            q_low = quantile(q500, 0.25),
            q_up = quantile(q500, 0.77)) %>%
  ggplot(aes(x =x)) +
  geom_ribbon(aes(ymin = q_low, ymax = q_up, fill = as.factor(clusters)), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.5, aes(y = q_med, color = as.factor(clusters))) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_wrap(~ clusters, ncol = 3) +
  geom_text(data = cluster_summaries, aes(x = x, y = y, label = label),
    inherit.aes = FALSE, hjust = 0, vjust = 1, size = 4) +
  theme_bw(base_size = 14) +
  scale_x_continuous(breaks = seq(1, 12, by = 3), 
                     minor_breaks = NULL, 
                     labels =  month.abb[seq(1, 12, by = 3)]) +
  scale_y_continuous(breaks = seq(0.5, 3, by = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none") +
  xlab("") +
  ylab("Pattern mean seasonal multiplier") +
  scale_color_manual(values = c("#E97F4A", "#628141", "#0992C2"))
 


p_prov <- ggplot(season_prov_smooth, aes(x = x, y = q500, colour = period, fill = period)) +
  geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.10, linewidth = 0, show.legend = FALSE) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ entity, ncol = 4) +
  scale_x_continuous(breaks = 1:12, minor_breaks = NULL) +
  labs(x = "Month", y = "Seasonal multiplier $\\exp(f_{seasonal}(m))$") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), legend.position = "bottom")

p_prov


