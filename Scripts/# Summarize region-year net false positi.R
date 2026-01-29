# Summarize region-year net false positive counts (total per year)
summarize_region_fp_net <- function(probs = c(0.025, 0.5, 0.975)) {
  stopifnot(exists("fit"), exists("stan_data"))
  Y <- stan_data$Y; R <- stan_data$R; S <- stan_data$S; T <- stan_data$T; A <- stan_data$A
  year_id  <- as.integer(stan_data$year_id)   # length T
  region_id <- as.integer(stan_data$region_id) # length S
  N_tas <- stan_data$N_tas                    # [T, A, S]

  dr <- rstan::extract(fit, pars = c("alpha_FP_yr", "delta_FP_s"), permuted = TRUE)
  alpha <- dr$alpha_FP_yr   # [draw, Y, R]
  d_s   <- dr$delta_FP_s    # [draw, S]
  n_draws <- dim(alpha)[1]

  # Precompute population totals per (year, province): sum over months and ages
  # pop_y_s[y, s] = sum_{t in year y} sum_a (N_tas[t,a,s] + 1e-9)
  pop_y_s <- matrix(0.0, nrow = Y, ncol = S)
  for (y in seq_len(Y)) {
    t_idx <- which(year_id == y)
    for (s in seq_len(S)) {
      pop_y_s[y, s] <- sum(N_tas[t_idx, , s] + 1e-9)
    }
  }

  # Compute net FP per draw for each (y, r)
  FP_net <- matrix(0.0, nrow = n_draws, ncol = Y * R)
  col_map <- expand.grid(y = seq_len(Y), r = seq_len(R))
  for (i in seq_len(n_draws)) {
    for (k in seq_len(nrow(col_map))) {
      y <- col_map$y[k]; r <- col_map$r[k]
      s_in_r <- which(region_id == r)
      # Sum over provinces in region with province deviations
      FP_net[i, k] <- sum(pop_y_s[y, s_in_r] * exp(alpha[i, y, r] + d_s[i, s_in_r]))
    }
  }

  qfun <- function(x) stats::quantile(x, probs = probs, names = FALSE, type = 7)

  out <- purrr::map_dfr(seq_len(nrow(col_map)), function(k) {
    v <- FP_net[, k]
    qs <- qfun(v)
    tibble::tibble(
      r_idx = col_map$r[k],
      y_idx = col_map$y[k],
      FP_net_mean = mean(v),
      FP_net_q025 = qs[1],
      FP_net_q500 = qs[2],
      FP_net_q975 = qs[3]
    )
  })

  # Attach year labels if available
  if (exists("calendar", inherits = TRUE)) {
    year_map <- calendar |>
      dplyr::distinct(y_idx, year) |>
      dplyr::arrange(y_idx)
    out <- dplyr::left_join(out, year_map, by = "y_idx")
  }

  # Region labels if available
#  if (exists("region_levels", inherits = TRUE)) {
#    out <- out |>
#      dplyr::mutate(region = factor(r_idx, levels = seq_len(R), labels = region_levels))
#  } else {
    out <- out |>
      dplyr::mutate(region = factor(r_idx, levels = seq_len(R), labels = paste0("Region ", seq_len(R))))
#  }

  out
}

# Plot net FP count shift by year per region (median + 95% ribbon)
plot_region_fp_net <- function(fp_df) {
  x_var <- if ("year" %in% names(fp_df)) "year" else "y_idx"
  ggplot(fp_df, aes(x = .data[[x_var]], y = FP_net_q500, color = region)) +
    geom_ribbon(aes(ymin = FP_net_q025, ymax = FP_net_q975, fill = region),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    labs(x = "Year", y = "False positives (net yearly count)",
         title = "Region false positive counts by year") +
    theme_bw(base_size = 16)
}

# Example
FP_region_net <- summarize_region_fp_net()
G_fp_net <- plot_region_fp_net(FP_region_net)
G_fp_net
