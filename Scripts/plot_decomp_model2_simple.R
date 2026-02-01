# Decompose reported totals into Se*I (true reported), (1-Se)*I (missed), and FP per province-year
# Reconstructs I_ts from C_ts using: I = (C - FP) / Se, where FP = (1-Sp)*N
summarize_province_report_decomposition_model2 <- function(probs = c(0.025, 0.5, 0.975)) {
  stopifnot(exists("fit"), exists("stan_data"))
  T <- stan_data$T; Y <- stan_data$Y; S <- stan_data$S
  year_id  <- as.integer(stan_data$year_id)    # [T]
  N_tas    <- stan_data$N_tas # %>% apply(c(1,3), sum)

  # Extract C_ts (expected observed cases), Se_y, Sp_y
  dr <- rstan::extract(fit, pars = c("C_ts", "Se_y", "Sp_y"), permuted = TRUE)

  # Draw arrays
  C_draw   <- dr$C_ts     # [draw, T, S]
  Se_draw  <- dr$Se_y     # [draw, Y]
  Sp_draw  <- dr$Sp_y     # [draw, Y]
  n_draws  <- dim(C_draw)[1]

  # Reconstruct I_ts from C_ts: C = I*Se + (1-Sp)*N => I = (C - (1-Sp)*N) / Se
  I_draw <- array(0.0, dim = c(n_draws, T, S))
  for (i in seq_len(n_draws)) {
    for (t in seq_len(T)) {
      y <- year_id[t]
      for (s in seq_len(S)) {
        FP_ts <- (1 - Sp_draw[i, y]) * N_tas[t, s]
        I_draw[i, t, s] <- pmax((C_draw[i, t, s] - FP_ts) / (Se_draw[i, y] + 1e-9), 0)
      }
    }
  }

  # True reported per month-province: Se_y[year] * I_ts
  true_rep_month <- array(0.0, dim = c(n_draws, T, S))
  FP_draw <- array(0.0, dim = c(n_draws, T, S))
  
  for (i in seq_len(n_draws)) {
    for (t in seq_len(T)) {
      y <- year_id[t]
      for (s in seq_len(S)) {
        true_rep_month[i, t, s] <- I_draw[i, t, s] * Se_draw[i, y]
        # FP per month-province: (1 - Sp_y) * N
        FP_draw[i, t, s] <- (1 - Sp_draw[i, y]) * N_tas[t, s]
      }
    }
  }

  # Missed per month-province: (1 - Se_y[year]) * I_ts = ((1-Se)/Se) * true_rep
  missed_month <- array(0.0, dim = c(n_draws, T, S))
  for (i in seq_len(n_draws)) {
    for (t in seq_len(T)) {
      y <- year_id[t]
      ratio <- (1 - Se_draw[i, y]) / (Se_draw[i, y] + 1e-9)
      for (s in seq_len(S)) {
        missed_month[i, t, s] <- pmax(ratio, 0) * true_rep_month[i, t, s]
      }
    }
  }

  # Aggregate per year-province totals
  col_map <- expand.grid(y = seq_len(Y), s = seq_len(S))
  true_rep <- matrix(0.0, nrow = n_draws, ncol = Y * S)
  missed   <- matrix(0.0, nrow = n_draws, ncol = Y * S)
  fp_tot   <- matrix(0.0, nrow = n_draws, ncol = Y * S)

  for (k in seq_len(nrow(col_map))) {
    y <- col_map$y[k]; s <- col_map$s[k]
    t_idx <- which(year_id == y)
    for (i in seq_len(n_draws)) {
      true_rep[i, k] <- sum(true_rep_month[i, t_idx, s])
      missed[i,   k] <- sum(missed_month[i, t_idx, s])
      fp_tot[i,   k] <- sum(FP_draw[i, t_idx, s])
    }
  }

  qfun <- function(x) stats::quantile(x, probs = probs, names = FALSE, type = 7)
  summarize_vec <- function(v) {
    qs <- qfun(v)
    tibble::tibble(mean = mean(v), q025 = qs[1], q500 = qs[2], q975 = qs[3])
  }

  df_true <- purrr::map_dfr(seq_len(nrow(col_map)), function(k) {
    summarize_vec(true_rep[, k]) |>
      dplyr::mutate(s_idx = col_map$s[k], y_idx = col_map$y[k], component = "True reported")
  })
  df_miss <- purrr::map_dfr(seq_len(nrow(col_map)), function(k) {
    summarize_vec(missed[, k]) |>
      dplyr::mutate(s_idx = col_map$s[k], y_idx = col_map$y[k], component = "Missed reports")
  })
  df_fp <- purrr::map_dfr(seq_len(nrow(col_map)), function(k) {
    summarize_vec(fp_tot[, k]) |>
      dplyr::mutate(s_idx = col_map$s[k], y_idx = col_map$y[k], component = "False positives")
  })

  out <- dplyr::bind_rows(df_true, df_miss, df_fp)

  # Attach year labels if available
  if (exists("calendar", inherits = TRUE)) {
    year_map <- calendar |>
      dplyr::distinct(y_idx, year) |>
      dplyr::arrange(y_idx)
    out <- dplyr::left_join(out, year_map, by = "y_idx")
  }

  # Province labels if available
  if (exists("site_levels", inherits = TRUE)) {
    out <- out |>
      dplyr::mutate(province = factor(s_idx, levels = seq_len(S), labels = site_levels))
  } else {
    out <- out |>
      dplyr::mutate(province = factor(s_idx, levels = seq_len(S), labels = paste0("Province ", seq_len(S))))
  }

  out
}

plot_province_report_decomposition_model2 <- function(decomp_df) {
  x_var <- if ("year" %in% names(decomp_df)) "year" else "y_idx"
  ggplot(decomp_df, aes(x = .data[[x_var]], y = q500, fill = component)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = c(
      "True reported" = "#1f78b4",
      "Missed reports" = "#a6cee3",
      "False positives" = "#e31a1c"
    )) +
    labs(x = "Year",
         y = "Counts per year",
         title = "Decomposition of reported totals by year and province (Model2_simple)",
         fill = NULL) +
    theme_bw(base_size = 16)
}

# Example usage:
decomp2 <- summarize_province_report_decomposition_model2()
G2 <- plot_province_report_decomposition_model2(decomp2)
G2 + facet_wrap(province ~ ., scales = "free_y")

# png(here("Figures", "province_report_decomposition_model2_simple.png"),
#     width = 12, height = 8, units = "in", res = 300)
# G2 + facet_wrap(province ~ ., scales = "free_y")
# dev.off()
