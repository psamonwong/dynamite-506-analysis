# Decompose reported totals into Se*I (true reported), (1-Se)*I (missed), and FP per region-year
summarize_region_report_decomposition <- function(probs = c(0.025, 0.5, 0.975)) {
  stopifnot(exists("fit"), exists("stan_data"))
  Y <- stan_data$Y; R <- stan_data$R; S <- stan_data$S; T <- stan_data$T; A <- stan_data$A
  year_id  <- as.integer(stan_data$year_id)    # [T]
  region_id <- as.integer(stan_data$region_id) # [S]
  N_tas <- stan_data$N_tas                     # [T, A, S]

  dr <- rstan::extract(fit,
                       pars = c("mu_case_ts", "alpha_Se_yr", "delta_Se_s",
                                "alpha_FP_yr", "delta_FP_s"),
                       permuted = TRUE)

  mu <- dr$mu_case_ts           # [draw, T, S]
  alphaSe <- dr$alpha_Se_yr     # [draw, Y, R]
  dSe_s   <- dr$delta_Se_s      # [draw, S]
  alphaFP <- dr$alpha_FP_yr     # [draw, Y, R]
  dFP_s   <- dr$delta_FP_s      # [draw, S]

  n_draws <- dim(mu)[1]

  # Precompute per-month, per-province sums of population over ages
  N_sum_t_s <- matrix(0.0, nrow = T, ncol = S)
  for (t in seq_len(T)) {
    for (s in seq_len(S)) {
      N_sum_t_s[t, s] <- sum(N_tas[t, , s] + 1e-9)
    }
  }

  # Containers for yearly-region totals per draw
  true_rep <- matrix(0.0, nrow = n_draws, ncol = Y * R)  # Se * I
  missed   <- matrix(0.0, nrow = n_draws, ncol = Y * R)  # (1-Se) * I
  fp_tot   <- matrix(0.0, nrow = n_draws, ncol = Y * R)  # FP

  col_map <- expand.grid(y = seq_len(Y), r = seq_len(R))
  for (i in seq_len(n_draws)) {
    # Per-month/province Se and FP (depend on year and province/region)
    Se_ts <- matrix(NA_real_, nrow = T, ncol = S)
    FP_ts <- matrix(0.0, nrow = T, ncol = S)
    for (s in seq_len(S)) {
      r <- region_id[s]
      se_y <- plogis(alphaSe[i, , r] + dSe_s[i, s])        # [Y]
      fp_y <- exp(alphaFP[i, , r] + dFP_s[i, s])           # [Y]
      for (t in seq_len(T)) {
        y <- year_id[t]
        Se_ts[t, s] <- se_y[y]
        FP_ts[t, s] <- fp_y[y] * N_sum_t_s[t, s]
      }
    }

    # Solve I[t,s] per month-province
    I_ts <- pmax((mu[i, , ] - FP_ts) / Se_ts, 0.0)  # [T, S]

    # Aggregate to yearly-region totals
    for (k in seq_len(nrow(col_map))) {
      y <- col_map$y[k]; r <- col_map$r[k]
      t_idx <- which(year_id == y)
      s_in_r <- which(region_id == r)
      # Sums over months in year and provinces in region
      true_rep[i, k] <- sum(Se_ts[t_idx, s_in_r] * I_ts[t_idx, s_in_r])
      missed[i, k]   <- sum((1 - Se_ts[t_idx, s_in_r]) * I_ts[t_idx, s_in_r])
      fp_tot[i, k]   <- sum(FP_ts[t_idx, s_in_r])
    }
  }

  qfun <- function(x) stats::quantile(x, probs = probs, names = FALSE, type = 7)

  summarize_vec <- function(v) {
    qs <- qfun(v)
    tibble::tibble(mean = mean(v), q025 = qs[1], q500 = qs[2], q975 = qs[3])
  }

  # Build tidy output
  df_true <- purrr::map_dfr(seq_len(nrow(col_map)), function(k) {
    summarize_vec(true_rep[, k]) |>
      dplyr::mutate(r_idx = col_map$r[k], y_idx = col_map$y[k], component = "True reported")
  })
  df_miss <- purrr::map_dfr(seq_len(nrow(col_map)), function(k) {
    summarize_vec(missed[, k]) |>
      dplyr::mutate(r_idx = col_map$r[k], y_idx = col_map$y[k], component = "Missed reports")
  })
  df_fp <- purrr::map_dfr(seq_len(nrow(col_map)), function(k) {
    summarize_vec(fp_tot[, k]) |>
      dplyr::mutate(r_idx = col_map$r[k], y_idx = col_map$y[k], component = "False positives")
  })

  out <- dplyr::bind_rows(df_true, df_miss, df_fp)

  # Attach year labels if available
  if (exists("calendar", inherits = TRUE)) {
    year_map <- calendar |>
      dplyr::distinct(y_idx, year) |>
      dplyr::arrange(y_idx)
    out <- dplyr::left_join(out, year_map, by = "y_idx")
  }

  # Region labels if available
  if (exists("region_levels", inherits = TRUE)) {
    out <- out |>
      dplyr::mutate(region = factor(r_idx, levels = seq_len(R), labels = region_levels))
  } else {
    out <- out |>
      dplyr::mutate(region = factor(r_idx, levels = seq_len(R), labels = paste0("Region ", seq_len(R))))
  }

  out
}

# Stacked bar plot of components using posterior medians (q500)
plot_region_report_decomposition <- function(decomp_df) {
  x_var <- if ("year" %in% names(decomp_df)) "year" else "y_idx"
  ggplot(decomp_df, aes(x = .data[[x_var]], y = q500, fill = component)) +
    geom_col(position = "stack") +
    facet_grid(rows = vars(region), scales = "free_y") +
    scale_fill_manual(values = c(
      "True reported" = "#1f78b4",
      "Missed reports" = "#a6cee3",
      "False positives" = "#e31a1c"
    )) +
    labs(x = "Year",
         y = "Counts per year",
         title = "Decomposition of reported totals by year and region",
         fill = NULL) +
    theme_bw(base_size = 16)
}

# Example
decomp <- summarize_region_report_decomposition()
G_decomp <- plot_region_report_decomposition(decomp)
G_decomp
