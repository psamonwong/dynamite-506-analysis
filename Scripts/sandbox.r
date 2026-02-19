 
ppc_cases_by_province <- function(s_idx) {
  stopifnot(exists("fit"), exists("calendar"), exists("S"), exists("T"), exists("site_levels"))
  if (!(s_idx %in% seq_len(S))) stop("s_idx out of range 1..S")

  # Extract draws from Model2_simple: C_ts (expected observed cases)
  C_mat <- as.matrix(fit, pars = "C_ts")
  phi_vec <- as.numeric(as.matrix(fit, pars = "phi_case"))
  if (is.null(colnames(C_mat))) stop("C_ts columns have no names; re-run sampling with pars='C_ts'.")

  # Parse indices from column names [t,s]
  idx <- stringr::str_match(colnames(C_mat), "C_ts\\[(\\d+),(\\d+)\\]")
  if (any(is.na(idx))) stop("Failed to parse C_ts indices from column names.")
  t_idx_col <- as.integer(idx[, 2])
  s_idx_col <- as.integer(idx[, 3])

  param_map <- tibble::tibble(par = colnames(C_mat), t_idx = t_idx_col, s = s_idx_col) |>
    dplyr::filter(s == s_idx)

  C_long <- tibble::as_tibble(C_mat[, param_map$par, drop = FALSE]) |>
    dplyr::mutate(.draw = dplyr::row_number()) |>
    tidyr::pivot_longer(cols = - .draw, names_to = "par", values_to = "C_expected") |>
    dplyr::left_join(param_map, by = "par")

  phi_df <- tibble::tibble(.draw = seq_len(nrow(C_mat)), phi = phi_vec)

  # Summarize PPC for this province over time
  pred_summ <- C_long |>
    dplyr::left_join(phi_df, by = ".draw") |>
    dplyr::group_by(t_idx) |>
    dplyr::summarize(
      mu_mean   = mean(C_expected),
      mu_median = median(C_expected),
      mu_q025   = stats::quantile(C_expected, 0.025, names = FALSE),
      mu_q975   = stats::quantile(C_expected, 0.975, names = FALSE),
      yrep_q025 = stats::quantile(stats::rnbinom(n = dplyr::n(), size = phi, mu = C_expected), 0.025, names = FALSE),
      yrep_q500 = stats::quantile(stats::rnbinom(n = dplyr::n(), size = phi, mu = C_expected), 0.50, names = FALSE),
      yrep_q975 = stats::quantile(stats::rnbinom(n = dplyr::n(), size = phi, mu = C_expected), 0.975, names = FALSE),
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

######################################################################################

library(ggpubr)
library(fs)

# Create output directory if it doesn't exist
dir_create(here("Plots", "predictive_check", dataset))
S <- stan_data$S
T <- stan_data$T
Y <- stan_data$Y
# Collect plots in a list

plot_list <- vector("list", S)
for (s_idx in seq_len(S)) {
  plot_data <- ppc_cases_by_province(s_idx = s_idx)
  plot_list[[s_idx]] <- ggplot(plot_data, aes(x = date)) +
    geom_ribbon(aes(ymin = yrep_q025, ymax = yrep_q975), fill = "steelblue", alpha = 0.20) +
    geom_line(aes(y = yrep_q500), color = "steelblue", linewidth = 1) +
    geom_point(aes(y = observed_case), color = "black", size = 1.6, alpha = 0.6) +
    labs(
      x = "", y = "Reported cases",
      title = unique(plot_data$province)  # Province name as title
    ) +
    theme_bw(base_size = 18)
}

# Arrange plots in 3x4 grids and export each page as PNG
plots_per_page <- 9
n_pages <- ceiling(S / plots_per_page)
for (i in seq_len(n_pages)) {
  idx_start <- (i - 1) * plots_per_page + 1
  idx_end <- min(i * plots_per_page, S)
  plots_this_page <- plot_list[idx_start:idx_end]
  arranged <- ggarrange(plotlist = plots_this_page, ncol = 3, nrow = 3)
  png(here("Plots", "predictive_check", dataset, sprintf("predictive_check_page%02d.png", i)), width = 16, height = 12, units = "in", res = 300)
  print(arranged)
  dev.off()
}
