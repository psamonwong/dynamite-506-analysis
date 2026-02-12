library(here)
library(tidyverse)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#############################
# MODULAR DATA PREPARATION FOR STAN MODELS
# Age-stratified model with age-specific susceptibility and case-fatality rate
#############################

# =============================================================================
# CONFIGURATION
# =============================================================================

config <- list(
  # Province selection (NULL = all provinces)
  areas = c("Chiang Rai", "Lampang", "Nakhon Ratchasima", "Surin", 
            "Songkhla", "Phatthalung", "Uttaradit", "Tak"),
  

  # Time period
  start_year = 2003,
  end_year = 2024,
  n_periods = 4,
  
  # Seasonality
  n_harmonics = 3,
  
  # Serology survey years
  sero_years = c(2007, 2008, 2012),
  
  # Fixed epidemiological parameters
  omega = 0.1882,           # immune waning rate
  nu = 0.15,                # clinical sign probability
  sigma_lambda_0 = 0.1,     # baseline FOI variation
  sigma_lambda_y = 0.122,   # year-to-year FOI variation
  sigma_lambda_init = 0.75, # initial FOI variation
  

  # Reference age group for susceptibility normalization
 ref_age_label = "25 to 35",
  
  # Age group definitions
  age_levels = c("0 to 5", "5 to 15", "15 to 25", "25 to 35",
                 "35 to 45", "45 to 55", "55 to 65", "above 65"),
  
  # Population age group mapping (raw -> standardized)
  pop_age_levels = c("0 to 4", "5 to 14", "15 to 24", "25 to 34",
                     "35 to 44", "45 to 54", "55 to 64", "above 65")
)

# =============================================================================
# DATA LOADING FUNCTIONS
# =============================================================================

load_case_data <- function(path = here("Analysis_data", "Scrub_typhus", "analysis_data.csv")) {
  read_csv(path, show_col_types = FALSE)
}

load_case_age_data <- function(path = here("Analysis_data", "Scrub_typhus", "case_death_by_age.csv")) {
  read_csv(path, show_col_types = FALSE)
}

load_population_data <- function(path = here("Analysis_data", "Population", "population_data_summarised.csv")) {
  read_csv(path, show_col_types = FALSE)
}

load_sero_data <- function(path = here("Analysis_data", "Gonwong_sero_data.csv")) {
  read_csv(path, show_col_types = FALSE)
}

# =============================================================================
# PERIOD ALLOCATION
# =============================================================================

#' Allocate years to period blocks
#' @param year_selected Vector of years
#' @param n_periods Number of period blocks
#' @return List with period_sizes and block_of_year mapping
allocate_periods <- function(year_selected, n_periods) {
  n_years <- length(year_selected)
  base_size <- n_years %/% n_periods
  remainder <- n_years %% n_periods
  
  period_sizes <- rep(base_size, n_periods)
  if (remainder > 0) {
    period_sizes[(n_periods - remainder + 1):n_periods] <- 
      period_sizes[(n_periods - remainder + 1):n_periods] + 1
  }
  
  block_of_year <- rep(1:n_periods, times = period_sizes)
  
  list(
    period_sizes = period_sizes,
    block_of_year = block_of_year,
    periods = split(year_selected, block_of_year)
  )
}

# =============================================================================
# DATA FILTERING AND PROCESSING
# =============================================================================

#' Filter and process case data
#' @param case_data Raw case data
#' @param areas Vector of province names to include
#' @param year_selected Vector of years to include
#' @return Filtered and arranged case data
filter_case_data <- function(case_data, areas, year_selected) {
  case_data |>
    filter(Area %in% areas, year %in% year_selected) |>
    arrange(Area, year, Month)
}

#' Filter and process age-stratified case data
#' @param case_age_data Raw age-stratified case data
#' @param areas Vector of province names
#' @param year_selected Vector of years
#' @param age_levels Ordered age group levels
#' @return Processed age-stratified data
filter_case_age_data <- function(case_age_data, areas, year_selected, age_levels) {
  case_age_data |>
    filter(Area %in% areas, year %in% year_selected, age_group_new != "unknown") |>
    group_by(Area, year, age_group_new) |>
    summarize(
      case = sum(case, na.rm = TRUE), 
      death = sum(death, na.rm = TRUE), 
      .groups = "drop"
    ) |>
    mutate(
      age_group_new = factor(age_group_new, levels = age_levels, ordered = TRUE)
    ) |>
    arrange(Area, year, age_group_new)
}

#' Filter and process population data
#' @param population_data Raw population data
#' @param areas Vector of province names
#' @param year_selected Vector of years
#' @param pop_age_levels Raw age group levels
#' @param age_levels Standardized age group labels
#' @return Processed population data
filter_population_data <- function(population_data, areas, year_selected, 
                                   pop_age_levels, age_levels) {
  population_data |>
    filter(Province %in% areas, year %in% year_selected) |>
    mutate(
      age_group_new = factor(
        age_group,
        levels = pop_age_levels,
        labels = age_levels,
        ordered = TRUE
      )
    )
}

# =============================================================================
# DIMENSION AND MAPPING EXTRACTION
# =============================================================================

#' Extract model dimensions and mappings
#' @param case_data_area Filtered case data
#' @param site_levels Province names
#' @param year_levels Year values
#' @param period_info Period allocation info
#' @return List of dimensions and mappings
extract_dimensions <- function(case_data_area, site_levels, year_levels, period_info) {
  S <- length(site_levels)
  Y <- length(year_levels)
  M <- 12L
  T <- Y * M
  P <- max(period_info$block_of_year)
  
  # Region mapping
  region_pairs <- case_data_area |>
    distinct(Area, Zone) |>
    arrange(match(Area, site_levels))
  region_levels <- region_pairs |> pull(Zone)
  R <- n_distinct(region_levels)
  region_id <- match(region_pairs$Zone, unique(region_levels))
  
  # Time mappings
  year_id <- rep(seq_len(Y), each = M)
  month_id <- rep(seq_len(M), times = Y)
  period_id <- period_info$block_of_year
  
  list(
    S = S, R = R, T = T, M = M, Y = Y, P = P,
    site_levels = site_levels,
    region_levels = region_levels,
    year_levels = year_levels,
    year_id = year_id,
    month_id = month_id,
    period_id = period_id,
    region_id = region_id,
    region_pairs = region_pairs
  )
}

# =============================================================================
# ARRAY CONSTRUCTION
# =============================================================================

#' Build calendar and monthly arrays
#' @param case_data_area Filtered case data
#' @param dims Dimensions object
#' @return List with calendar, R_case_ts, R_death_ts
build_monthly_arrays <- function(case_data_area, dims) {
  case_data_area <- case_data_area |>
    mutate(Area = factor(Area, levels = dims$site_levels))
  
  calendar <- expand_grid(
    Area = dims$site_levels, 
    year = dims$year_levels, 
    Month = 1:dims$M
  ) |>
    left_join(
      case_data_area |> select(Area, year, Month, case, death),
      by = c("Area", "year", "Month")
    ) |>
    mutate(
      case = coalesce(case, 0L), 
      death = coalesce(death, 0L),
      s = as.integer(factor(Area, levels = dims$site_levels)),
      y_idx = match(year, dims$year_levels),
      m_idx = Month,
      t_idx = (y_idx - 1L) * dims$M + m_idx
    )
  
  # Cases array [T, S]
  R_case_ts <- array(0L, dim = c(dims$S, dims$T))
  R_case_ts[cbind(calendar$s, calendar$t_idx)] <- as.integer(calendar$case)
  R_case_ts <- t(R_case_ts)
  
  # Deaths array [T, S]
  R_death_ts <- array(0L, dim = c(dims$S, dims$T))
  R_death_ts[cbind(calendar$s, calendar$t_idx)] <- as.integer(calendar$death)
  R_death_ts <- t(R_death_ts)
  
  list(calendar = calendar, R_case_ts = R_case_ts, R_death_ts = R_death_ts)
}

#' Build age-stratified arrays
#' @param case_age_data_area Filtered age-stratified data
#' @param dims Dimensions object
#' @param age_levels Age group levels
#' @return List with R_case_yas, R_death_yas, has_age_data_ys
build_age_arrays <- function(case_age_data_area, dims, age_levels) {
  A <- length(age_levels)
  
  case_age_data_area <- case_age_data_area |>
    mutate(Area = factor(Area, levels = dims$site_levels))
  
  age_comp <- case_age_data_area |>
    ungroup() |>
    mutate(year = as.integer(year)) |>
    complete(
      Area = dims$site_levels, 
      year = dims$year_levels, 
      age_group_new = age_levels,
      fill = list(case = 0L, death = 0L)
    ) |>
    mutate(
      a = as.integer(factor(age_group_new, levels = age_levels)),
      y_idx = match(year, dims$year_levels),
      s = as.integer(factor(Area, levels = dims$site_levels))
    )
  
  # Cases by age [Y, A, S]
  R_case_yas <- array(0L, dim = c(A, dims$Y, dims$S))
  R_case_yas[cbind(age_comp$a, age_comp$y_idx, age_comp$s)] <- as.integer(age_comp$case)
  R_case_yas <- aperm(R_case_yas, c(2, 1, 3))
  
  # Deaths by age [Y, A, S]
  R_death_yas <- array(0L, dim = c(A, dims$Y, dims$S))
  R_death_yas[cbind(age_comp$a, age_comp$y_idx, age_comp$s)] <- as.integer(age_comp$death)
  R_death_yas <- aperm(R_death_yas, c(2, 1, 3))
  
  # Indicator: has age data [Y, S]
  has_age_data_ys <- array(0L, dim = c(dims$Y, dims$S))
  age_data_summary <- case_age_data_area |>
    filter(case > 0) |>
    distinct(Area, year) |>
    mutate(
      y_idx = match(year, dims$year_levels),
      s_idx = match(Area, dims$site_levels)
    )
  
  for (i in seq_len(nrow(age_data_summary))) {
    has_age_data_ys[age_data_summary$y_idx[i], age_data_summary$s_idx[i]] <- 1L
  }
  
  list(
    A = A,
    age_levels = age_levels,
    R_case_yas = R_case_yas,
    R_death_yas = R_death_yas,
    has_age_data_ys = has_age_data_ys
  )
}

#' Build death data availability indicator
#' @param case_age_data Raw case age data (before filtering)
#' @param dims Dimensions object
#' @return has_death_data_ys array [Y, S]
build_death_data_indicator <- function(case_age_data, dims) {
  has_death_data_ys <- array(0L, dim = c(dims$Y, dims$S))
  
  death_data_summary <- case_age_data |>
    filter(Area %in% dims$site_levels, year %in% dims$year_levels) |>
    distinct(Area, year) |>
    mutate(
      y_idx = match(year, dims$year_levels),
      s_idx = match(Area, dims$site_levels)
    )
  
  for (i in seq_len(nrow(death_data_summary))) {
    has_death_data_ys[death_data_summary$y_idx[i], death_data_summary$s_idx[i]] <- 1L
  }
  
  has_death_data_ys
}

#' Build population array
#' @param population_data_area Filtered population data
#' @param dims Dimensions object
#' @param age_levels Age group levels
#' @return N_tas array [T, A, S]
build_population_array <- function(population_data_area, dims, age_levels) {
  A <- length(age_levels)
  
  # Find population column
  pop_candidates <- intersect(
    names(population_data_area), 
    c("pop", "population", "Population", "N", "n", "count", "Total")
  )
  if (length(pop_candidates) == 0) stop("No population column found.")
  pop_col <- pop_candidates[1]
  
  pop_tbl <- population_data_area |>
    rename(Area = Province, pop = all_of(pop_col)) |>
    mutate(
      Area = factor(Area, levels = dims$site_levels), 
      pop = as.integer(round(pop))
    ) |>
    filter(!is.na(pop)) |>
    complete(
      Area = dims$site_levels, 
      year = dims$year_levels, 
      age_group_new = age_levels,
      fill = list(pop = 0L)
    ) |>
    mutate(
      a = as.integer(factor(age_group_new, levels = age_levels)),
      y_idx = match(year, dims$year_levels),
      s = as.integer(factor(Area, levels = dims$site_levels))
    )
  
  N_ay <- array(0L, dim = c(A, dims$Y, dims$S))
  N_ay[cbind(pop_tbl$a, pop_tbl$y_idx, pop_tbl$s)] <- as.integer(pop_tbl$pop)
  N_ay <- aperm(N_ay, c(2, 1, 3))
  
  # Expand to monthly [T, A, S]
  N_tas <- N_ay[dims$year_id, , , drop = FALSE]
  
  N_tas
}

# =============================================================================
# FOURIER BASIS FUNCTIONS
# =============================================================================

#' Build Fourier basis matrices for seasonality
#' @param T Total number of time points
#' @param M Number of months per year
#' @param H Number of harmonics
#' @param month_id Month index for each time point
#' @return List with cos_t_h and sin_t_h matrices [T, H]
build_fourier_basis <- function(T, M, H, month_id) {
  cos_t_h <- matrix(0.0, nrow = T, ncol = H)
  sin_t_h <- matrix(0.0, nrow = T, ncol = H)
  
  for (t_ind in seq_len(T)) {
    theta <- 2 * pi * month_id[t_ind] / M
    for (h in seq_len(H)) {
      cos_t_h[t_ind, h] <- cos(h * theta)
      sin_t_h[t_ind, h] <- sin(h * theta)
    }
  }
  
  list(cos_t_h = cos_t_h, sin_t_h = sin_t_h)
}

# =============================================================================
# SEROLOGY DATA PROCESSING
# =============================================================================

#' Process serology data
#' @param sero_data Raw serology data
#' @param site_levels Province names
#' @param year_levels Year values
#' @param sero_years Years with serology surveys
#' @return List with sero_lower_s, sero_upper_s, sero_year_id
process_sero_data <- function(sero_data, site_levels, year_levels, sero_years) {
  sero_data_area <- sero_data |>
    filter(Province %in% site_levels) |>
    arrange(match(Province, site_levels))
  
  list(
    sero_lower_s = sero_data_area$sero_lower,
    sero_upper_s = sero_data_area$sero_upper,
    sero_year_id = match(sero_years, year_levels),
    K_sero = length(sero_years)
  )
}

# =============================================================================
# MAIN DATA PREPARATION FUNCTION
# =============================================================================

#' Prepare complete Stan data list
#' @param config Configuration list
#' @return List with stan_data and metadata
prepare_stan_data <- function(config) {
  # Load raw data
  case_data <- load_case_data()
  case_age_data <- load_case_age_data()
  population_data <- load_population_data()
  sero_data <- load_sero_data()
  
  # Determine areas
  if (is.null(config$areas)) {
    areas <- case_data |> distinct(Area) |> pull(Area)
  } else {
    areas <- config$areas
  }
  
  # Year selection
  year_selected <- config$start_year:config$end_year
  
  # Period allocation
  period_info <- allocate_periods(config$start_year, config$end_year, config$n_periods)
  
  # Filter data
  case_data_area <- filter_case_data(case_data, areas, year_selected)
  case_age_data_area <- filter_case_age_data(
    case_age_data, areas, year_selected, config$age_levels
  )
  population_data_area <- filter_population_data(
    population_data, areas, year_selected, 
    config$pop_age_levels, config$age_levels
  )
  
  # Extract dimensions
  dims <- extract_dimensions(case_data_area, areas, year_selected, period_info)
  
  # Build arrays
  monthly <- build_monthly_arrays(case_data_area, dims)
  age_arrays <- build_age_arrays(case_age_data_area, dims, config$age_levels)
  has_death_data_ys <- build_death_data_indicator(case_age_data, dims)
  N_tas <- build_population_array(population_data_area, dims, config$age_levels)
  fourier <- build_fourier_basis(dims$T, dims$M, config$n_harmonics, dims$month_id)
  sero <- process_sero_data(sero_data, dims$site_levels, dims$year_levels, config$sero_years)
  
  # Reference age index
  ref_age <- which(config$age_levels == config$ref_age_label)
  
  # Assemble Stan data
  stan_data <- list(
    # Dimensions
    S = dims$S,
    R = dims$R,
    T = dims$T,
    M = dims$M,
    Y = dims$Y,
    A = age_arrays$A,
    P = dims$P,
    K_sero = sero$K_sero,
    
    # Mappings
    year_id = dims$year_id,
    month_id = dims$month_id,
    period_id = dims$period_id,
    region_id = dims$region_id,
    
    # Case data
    R_case_ts = monthly$R_case_ts,
    R_case_yas = age_arrays$R_case_yas,
    has_age_data_ys = age_arrays$has_age_data_ys,
    
    # Death data
    R_death_ts = monthly$R_death_ts,
    R_death_yas = age_arrays$R_death_yas,
    has_death_data_ys = has_death_data_ys,
    
    # Population
    N_tas = N_tas,
    
    # Seasonality
    H = config$n_harmonics,
    cos_t_h = fourier$cos_t_h,
    sin_t_h = fourier$sin_t_h,
    
    # Reference age group
    ref_age = ref_age,
    
    # Serology data
    sero_lower_s = sero$sero_lower_s,
    sero_upper_s = sero$sero_upper_s,
    sero_year_id = sero$sero_year_id,
    
    # Fixed parameters
    omega = config$omega,
    nu = config$nu,
    sigma_lambda_0 = config$sigma_lambda_0,
    sigma_lambda_y = config$sigma_lambda_y,
    sigma_lambda_init = config$sigma_lambda_init
  )
  
  # Return data and metadata
  list(
    stan_data = stan_data,
    calendar = monthly$calendar,
    site_levels = dims$site_levels,
    region_levels = dims$region_levels,
    age_levels = age_arrays$age_levels,
    year_levels = dims$year_levels,
    period_info = period_info,
    config = config
  )
}

# =============================================================================
# EXECUTE DATA PREPARATION
# =============================================================================

# Prepare data using configuration
result <- prepare_stan_data(config)

# Extract for convenience
stan_data <- result$stan_data
calendar <- result$calendar
site_levels <- result$site_levels
region_levels <- result$region_levels
age_levels <- result$age_levels

# Save prepared data
save(
  stan_data, calendar, site_levels, region_levels, age_levels,
  file = here("Rout", "stan_data_09.RData")
)

# Print summary
cat("Data preparation complete!\n")
cat("Dimensions:\n")
cat(sprintf("  Provinces (S): %d\n", stan_data$S))
cat(sprintf("  Regions (R): %d\n", stan_data$R))
cat(sprintf("  Years (Y): %d\n", stan_data$Y))
cat(sprintf("  Months (T): %d\n", stan_data$T))
cat(sprintf("  Age groups (A): %d\n", stan_data$A))
cat(sprintf("  Periods (P): %d\n", stan_data$P))
cat(sprintf("  Harmonics (H): %d\n", stan_data$H))
