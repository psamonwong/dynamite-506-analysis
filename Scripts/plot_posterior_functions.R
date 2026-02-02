library(here)
library(tidyverse)
library(rstan)

load(here("Rout", "fit_03_test.RData"))  # loads 'fit', 'stan_data', 'site_levels', 'calendar'
load(here("Rout", "stan_data_03.RData"))  # loads 'stan_data', 'calendar', 'site_levels', etc.

s_idx <- 1  # Example province index for testing
S <- stan_data$S
T <- stan_data$T
Y <- stan_data$Y

 
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

library(ggpubr)
library(fs)

# Create output directory if it doesn't exist
dir_create("Plots/predictive_check")
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
plots_per_page <- 12
n_pages <- ceiling(S / plots_per_page)
for (i in seq_len(n_pages)) {
  idx_start <- (i - 1) * plots_per_page + 1
  idx_end <- min(i * plots_per_page, S)
  plots_this_page <- plot_list[idx_start:idx_end]
  arranged <- ggarrange(plotlist = plots_this_page, ncol = 4, nrow = 3)
  ggsave(
    filename = sprintf("Plots/predictive_check/predictive_check_page%02d.png", i),
    plot = arranged,
    width = 16, height = 12, units = "in", dpi = 300
  )
}





# ------------------------------------------------------------------------------
# Extract and summarize baseline FOI (lambda_base_ys) over years and provinces
# Returns posterior summaries with year and province labels
extract_lambda_base <- function(probs = c(0.025, 0.5, 0.975)) {
  stopifnot(exists("fit"), exists("stan_data"), exists("site_levels"), exists("calendar"))
  S <- stan_data$S
  Y <- stan_data$Y
  R <- stan_data$R

  # Extract lambda_base_ys draws [Y, S]
  lambda_mat <- as.matrix(fit, pars = "lambda_base_ys")
  if (is.null(colnames(lambda_mat))) stop("lambda_base_ys columns have no names; re-run sampling with pars='lambda_base_ys'.")

  # Parse indices from column names [y,s]
  idx <- stringr::str_match(colnames(lambda_mat), "lambda_base_ys\\[(\\d+),(\\d+)\\]")
  if (any(is.na(idx))) stop("Failed to parse lambda_base_ys indices from column names.")
  y_idx_col <- as.integer(idx[, 2])
  s_idx_col <- as.integer(idx[, 3])

  # Create parameter map
  param_map <- tibble::tibble(
    par = colnames(lambda_mat),
    y_idx = y_idx_col,
    s_idx = s_idx_col
  )

  # Convert to long format
  lambda_long <- tibble::as_tibble(lambda_mat) |>
    dplyr::mutate(.draw = dplyr::row_number()) |>
    tidyr::pivot_longer(cols = -.draw, names_to = "par", values_to = "lambda_base") |>
    dplyr::left_join(param_map, by = "par")

  # Summarize by year and province
  lambda_summ <- lambda_long |>
    dplyr::group_by(y_idx, s_idx) |>
    dplyr::summarize(
      mean = mean(lambda_base),
      q025 = stats::quantile(lambda_base, probs[1], names = FALSE),
      q500 = stats::quantile(lambda_base, probs[2], names = FALSE),
      q975 = stats::quantile(lambda_base, probs[3], names = FALSE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      province = site_levels[s_idx],
      year = min(calendar$year) + y_idx - 1  # Map y_idx to actual year
    )

  # Extract region info from stan_data
  region_id <- stan_data$region_id  # [S] province to region mapping
  region_map <- tibble::tibble(
    s_idx = seq_len(S),
    r_idx = region_id[s_idx],
    region = paste0("Zone", sprintf("%02d", region_id[s_idx]))
  )
  
  lambda_summ <- lambda_summ |>
    dplyr::left_join(region_map, by = "s_idx")

  lambda_summ
}

# Plot baseline FOI over time, faceted by region, lines by province
plot_lambda_base_by_region <- function(lambda_data = NULL) {
  if (is.null(lambda_data)) {
    lambda_data <- extract_lambda_base()
  }

  lambda_data |>
    dplyr::mutate(
      yearly_mean  = (1 - exp(-12 * mean)) * 1e5,
      yearly_q025  = (1 - exp(-12 * q025)) * 1e5,
      yearly_q975  = (1 - exp(-12 * q975)) * 1e5
    ) |>
    ggplot(aes(x = year, y = yearly_mean, group = province)) +
    geom_ribbon(aes(ymin = yearly_q025, ymax = yearly_q975), alpha = 0.15, fill = "gray40") +
    geom_line(linewidth = 0.8, color = "gray20", alpha = 0.7) +
    labs(
      x = "Year",
      y = "Yearly incidence per 100,000",
      title = "Baseline Force of Infection → Yearly Incidence",
      subtitle = "Faceted by region, each line represents a province"
    ) +
    theme_bw(base_size = 14) +
    theme(legend.position = "none") +
    facet_wrap(~ region, scales = "free_y", ncol = 4)
}

plot_lambda_base_by_region()

pname <- here("Plots", "lambda_base_by_region.png")
png(filename = pname, width = 12, height = 9, units = "in", res = 300)
plot_lambda_base_by_region()
dev.off()


################################################################

# Extract and plot Se_y over time
plot_Se_over_time <- function(probs = c(0.025, 0.5, 0.975)) {
  stopifnot(exists("fit"), exists("stan_data"), exists("calendar"))
  Y <- stan_data$Y
  
  # Extract Se_y draws [Y]
  Se_mat <- as.matrix(fit, pars = "Se_y")
  if (is.null(colnames(Se_mat))) stop("Se_y columns have no names; re-run sampling with pars='Se_y'.")
  
  # Parse indices from column names [y]
  idx <- stringr::str_match(colnames(Se_mat), "Se_y\\[(\\d+)\\]")
  if (any(is.na(idx))) stop("Failed to parse Se_y indices from column names.")
  y_idx_col <- as.integer(idx[, 2])
  
  # Create parameter map
  param_map <- tibble::tibble(
    par = colnames(Se_mat), 
    y_idx = y_idx_col
  )
  
  # Convert to long format
  Se_long <- tibble::as_tibble(Se_mat) |>
    dplyr::mutate(.draw = dplyr::row_number()) |>
    tidyr::pivot_longer(cols = -.draw, names_to = "par", values_to = "Se") |>
    dplyr::left_join(param_map, by = "par")
  
  # Summarize by year
  Se_summ <- Se_long |>
    dplyr::group_by(y_idx) |>
    dplyr::summarize(
      mean = mean(Se),
      q025 = stats::quantile(Se, probs[1], names = FALSE),
      q500 = stats::quantile(Se, probs[2], names = FALSE),
      q975 = stats::quantile(Se, probs[3], names = FALSE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      year = min(calendar$year) + y_idx - 1  # Map y_idx to actual year
    )
  
  # Plot
  ggplot(Se_summ, aes(x = year)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), fill = "darkgreen", alpha = 0.3) +
    geom_line(aes(y = mean), color = "darkgreen", linewidth = 1.2) +
    geom_point(aes(y = q500), color = "darkgreen", size = 2) +
    labs(
      x = "Year",
      y = "Case Detection Sensitivity (Se)",
      title = "Temporal Evolution of Case Detection Sensitivity",
      subtitle = "Posterior mean with 95% credible intervals"
    ) +
    scale_y_continuous(limits = c(0, 0.3), labels = scales::percent) +
    theme_bw(base_size = 14) +
    theme(panel.grid.minor = element_blank())
}

pname <- here("Plots", "Se_over_time.png")
png(filename = pname, width = 6, height = 6, units = "in", res = 300)
plot_Se_over_time()
dev.off()

################################################################


















# ------------------------------------------------------------------------------
# Monthly totals of true incident cases (I_ts[t,s])
# In Model2_simple, I_ts is directly available in transformed parameters.
# Returns posterior summaries by month, optionally joined to calendar.
derive_I_clinical_totals <- function(s_idx, probs = c(0.025, 0.5, 0.975)) {
  stopifnot(exists("fit"), exists("stan_data"))
  S <- stan_data$S; T <- stan_data$T
  if (!(s_idx %in% seq_len(S))) stop("s_idx out of range 1..S")

  # Extract I_ts draws from Model2_simple
  I_mat <- as.matrix(fit, pars = "I_ts")
  if (is.null(colnames(I_mat))) stop("I_ts columns have no names; re-run sampling with pars='I_ts'.")

  # Parse indices from column names [t,s]
  idx <- stringr::str_match(colnames(I_mat), "I_ts\\[(\\d+),(\\d+)\\]")
  if (any(is.na(idx))) stop("Failed to parse I_ts indices from column names.")
  t_idx_col <- as.integer(idx[, 2])
  s_idx_col <- as.integer(idx[, 3])

  # Filter for this province
  cols_for_s <- which(s_idx_col == s_idx)
  I_draws <- I_mat[, cols_for_s, drop = FALSE]
  t_indices <- t_idx_col[cols_for_s]

  summarize_t <- function(v) {
    qs <- stats::quantile(v, probs = probs, names = FALSE, type = 7)
    tibble::tibble(I_mean = mean(v), I_q025 = qs[1], I_q500 = qs[2], I_q975 = qs[3])
  }

  out <- purrr::map_dfr(seq_len(T), function(t) {
    col_idx <- which(t_indices == t)
    if (length(col_idx) == 0) return(NULL)
    summarize_t(I_draws[, col_idx]) |> dplyr::mutate(t_idx = t, .before = 1)
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


s_toplot <- 1
I_monthly <- derive_I_clinical_totals(s_idx = s_toplot)
plot_data <- ppc_cases_by_province(s_idx = s_toplot)

ggplot(plot_data, aes(x = date)) +
  geom_ribbon(aes(ymin = yrep_q025, ymax = yrep_q975), fill = "steelblue", alpha = 0.20) +
  geom_line(aes(y = yrep_q500), color = "steelblue", linewidth = 1) +
  geom_point(aes(y = observed_case), color = "black", size = 1.6, alpha = 0.6) +
  #Overlay expected new clinical infections (monthly totals)
  geom_ribbon(data = I_monthly,
              aes(x = date, ymin = I_q025, ymax = I_q975), fill = "firebrick",
                  alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = I_monthly,
            aes(x = date, y = I_q500), color = "firebrick", linewidth = 1,
                inherit.aes = FALSE) +
  labs(x = "", y = "Reported cases",
       title = "Prediction vs observed reported cases",
       subtitle = unique(plot_data$province)) +
  theme_bw(base_size = 18)

