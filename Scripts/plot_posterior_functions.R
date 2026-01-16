library(here)
library(tidyverse)
library(rstan)

load(here("Rout", "fit_01.RData"))  # loads 'fit', 'stan_data', 'site_levels'

s_idx <- 1  # Example province index for testing
S <- stan_data$S
T <- stan_data$T


ppc_cases_by_province <- function(s_idx) {
  stopifnot(exists("fit"), exists("calendar"), exists("S"), exists("T"), exists("site_levels"))
  if (!(s_idx %in% seq_len(S))) stop("s_idx out of range 1..S")

  # Extract draws
  mu_mat <- as.matrix(fit, pars = "mu_case_ts")
  phi_vec <- as.numeric(as.matrix(fit, pars = "phi_c"))
  if (is.null(colnames(mu_mat))) stop("mu_case_ts columns have no names; re-run sampling with pars='mu_case_ts'.")

  # Parse indices from column names, handle either [t,s] or [s,t]
  idx <- stringr::str_match(colnames(mu_mat), "mu_case_ts\\[(\\d+),(\\d+)\\]")
  if (any(is.na(idx))) stop("Failed to parse mu_case_ts indices from column names.")
  i1 <- as.integer(idx[, 2]); i2 <- as.integer(idx[, 3])

  if (max(i1, na.rm = TRUE) == T && max(i2, na.rm = TRUE) == S) {
    t_idx_col <- i1; s_idx_col <- i2
  } else if (max(i1, na.rm = TRUE) == S && max(i2, na.rm = TRUE) == T) {
    s_idx_col <- i1; t_idx_col <- i2
  } else {
    stop("mu_case_ts index sizes do not match expected T and S.")
  }

  param_map <- tibble::tibble(par = colnames(mu_mat), t_idx = t_idx_col, s = s_idx_col) |>
    dplyr::filter(s == s_idx)

  mu_long <- tibble::as_tibble(mu_mat[, param_map$par, drop = FALSE]) |>
    dplyr::mutate(.draw = dplyr::row_number()) |>
    tidyr::pivot_longer(cols = - .draw, names_to = "par", values_to = "mu") |>
    dplyr::left_join(param_map, by = "par")

  phi_df <- tibble::tibble(.draw = seq_len(nrow(mu_mat)), phi = phi_vec)

  # Summarize PPC for this province over time
  pred_summ <- mu_long |>
    dplyr::left_join(phi_df, by = ".draw") |>
    dplyr::group_by(t_idx) |>
    dplyr::summarize(
      mu_mean   = mean(mu),
      mu_median = median(mu),
      mu_q025   = stats::quantile(mu, 0.025, names = FALSE),
      mu_q975   = stats::quantile(mu, 0.975, names = FALSE),
      yrep_q025 = stats::quantile(stats::rnbinom(n = dplyr::n(), size = phi, mu = mu), 0.025, names = FALSE),
      yrep_q500 = stats::quantile(stats::rnbinom(n = dplyr::n(), size = phi, mu = mu), 0.50, names = FALSE),
      yrep_q975 = stats::quantile(stats::rnbinom(n = dplyr::n(), size = phi, mu = mu), 0.975, names = FALSE),
      .groups = "drop"
    )

  # Observed cases for this province
  obs <- calendar |>
    dplyr::filter(s == s_idx) |>
    dplyr::select(t_idx, Area, year, Month, observed_case = case)

  # Final PPC dataset for the province
  obs |>
    dplyr::left_join(pred_summ, by = "t_idx") |>
    dplyr::arrange(year, Month) |>
    dplyr::mutate(province_index = s_idx, province = site_levels[s_idx]) |>
    dplyr::mutate(date = as.Date(sprintf("%d-%02d-01", year, Month)))
}
######################################################################################
# List available time-site parameters from the fit (use this to choose par_name explicitly)
list_ts_params <- function(fit, T, S, A = NULL) {
  avail <- fit@sim$fnames_oi
  base <- sub("\\[.*\\]$", "", avail)
  tibble(name = unique(base)) |>
    mutate(
      arity = case_when(
        any(grepl(paste0("^", name, "\\[\\d+,\\d+\\]$"), avail)) ~ 2L,
        any(grepl(paste0("^", name, "\\[\\d+,\\d+,\\d+\\]$"), avail)) ~ 3L,
        TRUE ~ NA_integer_
      )
    ) |>
    filter(!is.na(arity)) |>
    # keep only those whose indices match T,S,(A) by checking a few columns
    rowwise() |>
    mutate(matches_dims = {
      cols <- grep(paste0("^", name, "\\["), avail, value = TRUE)
      idx <- stringr::str_match(cols, "\\[(\\d+)(?:,(\\d+))?(?:,(\\d+))?\\]")
      dims <- suppressWarnings(apply(idx[, -1, drop = FALSE], 2, function(x) max(as.integer(x), na.rm = TRUE)))
      if (arity == 2) identical(sort(na.omit(dims)), sort(c(T, S))) else identical(sort(na.omit(dims)), sort(c(T, A, S)))
    }) |>
    ungroup() |>
    filter(matches_dims) |>
    arrange(arity, name)
}######################################################################################
list_ts_params(fit, T = T, S = S)


















plot_data <- ppc_cases_by_province(s_idx = 4)
for(s_idx in 1:S){
   G <- ppc_cases_by_province(s_idx = s_idx) |>
     ggplot(aes(x = date)) +
        geom_ribbon(aes(ymin = yrep_q025, ymax = yrep_q975), fill = "steelblue", alpha = 0.20) +
        geom_line(aes(y = yrep_q500), color = "steelblue", linewidth = 1) +
        geom_point(aes(y = observed_case), color = "black", size = 1.6, alpha = 0.6) +
        labs(x = "", y = "Reported cases",
            title = "Prediction vs observed reported cases",
            subtitle = unique(plot_data$province)) +
        theme_bw(base_size = 18) 
     print(G)
}



ggplot(plot_data, aes(x = date)) +
  geom_ribbon(aes(ymin = yrep_q025, ymax = yrep_q975), fill = "steelblue", alpha = 0.20) +
  geom_line(aes(y = yrep_q500), color = "steelblue", linewidth = 1) +
  geom_point(aes(y = observed_case), color = "black", size = 1.6, alpha = 0.6) +
  labs(x = "", y = "Reported cases",
       title = "Prediction vs observed reported cases",
       subtitle = unique(plot_data$province)) +
  theme_bw(base_size = 18)


# ------------------------------------------------------------------------------
# Monthly totals of expected new clinical infections (I_total[t,s])
# Derived from saved draws without needing age-specific parameters.
# Uses identity: mu_case_ts[t,s] = I_total[t,s] * Se[t,s] + FP_total[t,s]
# where Se[t,s] = logistic(alpha_Se_yr[y(t), r(s)] + delta_Se_s[s]) and
# FP_total[t,s] = sum_a exp(log(N_tas[t,a,s] + 1e-9) + alpha_FP_yr[y(t), r(s)] + delta_FP_s[s])
# Returns posterior summaries by month, optionally joined to calendar.
derive_I_clinical_totals <- function(s_idx, probs = c(0.025, 0.5, 0.975)) {
  stopifnot(exists("fit"), exists("stan_data"))
  S <- stan_data$S; T <- stan_data$T; Y <- stan_data$Y; R <- stan_data$R; A <- stan_data$A
  if (!(s_idx %in% seq_len(S))) stop("s_idx out of range 1..S")

  y_of_t <- as.integer(stan_data$year_id)   # length T, values 1..Y
  r_of_s <- as.integer(stan_data$region_id)[s_idx] # single region index for s_idx
  N_tas  <- stan_data$N_tas                 # [T, A, S]

  dr <- rstan::extract(
    fit,
    pars = c("mu_case_ts", "alpha_Se_yr", "delta_Se_s", "alpha_FP_yr", "delta_FP_s"),
    permuted = TRUE
  )

  # Shapes:
  # mu_case_ts: [draw, T, S]
  # alpha_Se_yr: [draw, Y, R]
  # delta_Se_s: [draw, S]
  # alpha_FP_yr: [draw, Y, R]
  # delta_FP_s: [draw, S]
  mu_ts   <- dr$mu_case_ts[, , s_idx, drop = FALSE][, , 1]
  se_y_r  <- dr$alpha_Se_yr
  dse_s   <- dr$delta_Se_s[, s_idx]
  afp_y_r <- dr$alpha_FP_yr
  dfp_s   <- dr$delta_FP_s[, s_idx]

  n_draws <- nrow(mu_ts)
  Se_draws <- matrix(NA_real_, nrow = n_draws, ncol = T)
  FP_draws <- matrix(0.0, nrow = n_draws, ncol = T)

  # Precompute population sums over age per month for this province
  N_sum_t <- vapply(seq_len(T), function(t) sum(N_tas[t, , s_idx] + 1e-9), numeric(1))

  for (t in seq_len(T)) {
    y <- y_of_t[t]
    r <- r_of_s
    # Se[t,s] shared across ages
    Se_draws[, t] <- plogis(se_y_r[, y, r] + dse_s)
    # FP_total[t,s] = exp(alpha_FP[y,r] + delta_FP[s]) * sum_a (N_tas[t,a,s] + 1e-9)
    FP_draws[, t] <- exp(afp_y_r[, y, r] + dfp_s) * N_sum_t[t]
  }

  # I_total draws; clamp negatives to 0 for numerical safety
  I_tot <- pmax((mu_ts - FP_draws) / Se_draws, 0)

  summarize_t <- function(v) {
    qs <- stats::quantile(v, probs = probs, names = FALSE, type = 7)
    tibble::tibble(I_mean = mean(v), I_q025 = qs[1], I_q500 = qs[2], I_q975 = qs[3])
  }

  out <- purrr::map_dfr(seq_len(T), function(t) {
    summarize_t(I_tot[, t]) |> dplyr::mutate(t_idx = t, .before = 1)
  }) |> dplyr::mutate(s = s_idx)

  # Attach calendar if available
  if (exists("calendar", inherits = TRUE)) {
    cal_s <- calendar |>
      dplyr::filter(s == s_idx) |>
      dplyr::select(t_idx, Area, year, Month) |>
      dplyr::mutate(date = as.Date(sprintf("%d-%02d-01", year, Month)))
    out <- dplyr::left_join(out, cal_s, by = "t_idx")
  }

  out
}


s_toplot <- 8
I_monthly <- derive_I_clinical_totals(s_idx = s_toplot)
plot_data <- ppc_cases_by_province(s_idx = s_toplot)

ggplot(plot_data, aes(x = date)) +
  geom_ribbon(aes(ymin = yrep_q025, ymax = yrep_q975), fill = "steelblue", alpha = 0.20) +
  geom_line(aes(y = yrep_q500), color = "steelblue", linewidth = 1) +
  geom_point(aes(y = observed_case), color = "black", size = 1.6, alpha = 0.6) +
  #Overlay expected new clinical infections (monthly totals)
  geom_ribbon(
    data = I_monthly,
    aes(x = date, ymin = I_q025, ymax = I_q975),
    fill = "firebrick",
    alpha = 0.15,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = I_monthly,
    aes(x = date, y = I_q500),
    color = "firebrick",
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  labs(x = "", y = "Reported cases",
       title = "Prediction vs observed reported cases",
       subtitle = unique(plot_data$province)) +
  theme_bw(base_size = 18)

