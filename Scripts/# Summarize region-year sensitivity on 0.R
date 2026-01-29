# Summarize region-year sensitivity on 0–100% scale
summarize_region_sensitivity <- function(probs = c(0.025, 0.5, 0.975)) {
  stopifnot(exists("fit"), exists("stan_data"))
  Y <- stan_data$Y; R <- stan_data$R

  dr <- rstan::extract(fit, pars = c("alpha_Se_yr"), permuted = TRUE)
  alpha <- dr$alpha_Se_yr  # [draw, Y, R]
  n_draws <- dim(alpha)[1]

  # Se(y,r) = inv_logit(alpha_Se_yr[y,r]); summarize to 0–100%
  qfun <- function(x) stats::quantile(x, probs = probs, names = FALSE, type = 7)
  out <- purrr::map_dfr(seq_len(R), function(r) {
    purrr::map_dfr(seq_len(Y), function(y) {
      v <- plogis(alpha[, y, r]) * 100
      qs <- qfun(v)
      tibble::tibble(
        r_idx = r,
        y_idx = y,
        Se_pct_mean = mean(v),
        Se_pct_q025 = qs[1],
        Se_pct_q500 = qs[2],
        Se_pct_q975 = qs[3]
      )
    })
  })

  # Attach year labels if calendar is available
  if (exists("calendar", inherits = TRUE)) {
    year_map <- calendar |>
      dplyr::distinct(y_idx, year) |>
      dplyr::arrange(y_idx)
    out <- dplyr::left_join(out, year_map, by = "y_idx")
  }

  # Region labels if available; else fallback to "Region <r>"
 # if (exists("region_levels", inherits = TRUE)) {
 #   out <- out |>
 #     dplyr::mutate(region = factor(r_idx, levels = seq_len(R), labels = region_levels))
 # } else {
    out <- out |>
      dplyr::mutate(region = factor(r_idx, levels = seq_len(R), labels = paste0("Region ", seq_len(R))))
#  }

  out
}

# Plot sensitivity shift by year per region (median + 95% ribbon)
plot_region_sensitivity <- function(se_df) {
  x_var <- if ("year" %in% names(se_df)) "year" else "y_idx"
  ggplot(se_df, aes(x = .data[[x_var]], y = Se_pct_q500, color = region)) +
    geom_ribbon(aes(ymin = Se_pct_q025, ymax = Se_pct_q975, fill = region),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    labs(x = "Year", y = "Sensitivity (%)", title = "Region case detection sensitivity by year") +
    theme_bw(base_size = 16)
}

# Example usage
Se_region <- summarize_region_sensitivity()
G_se <- plot_region_sensitivity(Se_region)
G_se +
  facet_grid(rows = vars(region)) 
