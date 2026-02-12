library(here)
library(tidyverse)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#############################
# DATA PREPARATION FOR MODEL 9
# Age-stratified model with age-specific susceptibility and case-fatality rate
#############################

# Load datasets
case_data <- read_csv(here("Analysis_data", "Scrub_typhus", "analysis_data.csv"))
case_age_data <- read_csv(here("Analysis_data", "Scrub_typhus", "case_death_by_age.csv"))
population_data <- read_csv(here("Analysis_data", "Population", "population_data_summarised.csv"))
sero_data <- read_csv(here("Analysis_data", "Gonwong_sero_data.csv"))

# Setup
area <- #case_data %>% distinct(Area) %>% pull(Area)
c( "Chiang Rai", "Lampang",  "Nakhon Ratchasima", "Surin", "Songkhla", "Phatthalung", "Uttaradit", "Tak")

start_year <- 2003
end_year   <- 2024
n_periods  <- 4

year_selected <- start_year:end_year
n_years <- length(year_selected)

# Base size and remainder
base_size <- n_years %/% n_periods
remainder <- n_years %% n_periods

# Distribute remainder to later periods
period_sizes <- rep(base_size, n_periods)
period_sizes[(n_periods - remainder + 1):n_periods] <- 
  period_sizes[(n_periods - remainder + 1):n_periods] + 1

periods <- split(
  year_selected,
  rep(1:n_periods, times = period_sizes)
)

# Filter case data  
case_data_area <- case_data %>%
  filter(Area %in% area, year %in% year_selected) %>%
  arrange(Area, year, Month)

# Filter unknown ages early and aggregate by age group
# Include both case and death counts
case_age_data_area <- case_age_data %>%
  filter(Area %in% area, year %in% year_selected, age_group_new != "unknown") %>%
  group_by(Area, year, age_group_new) %>%
  summarize(case = sum(case, na.rm = TRUE), death = sum(death, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    age_group_new = factor(
      age_group_new,
      levels = c("0 to 5", "5 to 15", "15 to 25", "25 to 35",
                 "35 to 45", "45 to 55", "55 to 65", "above 65"),
      ordered = TRUE
    )
  ) %>%
  arrange(Area, year, age_group_new)

population_data_area <- population_data %>%
  filter(Province %in% area, year %in% year_selected) %>%
  mutate(
    age_group_new = factor(
      age_group,
      levels = c("0 to 4", "5 to 14", "15 to 24", "25 to 34",
                 "35 to 44", "45 to 54", "55 to 64", "above 65"),
      labels = c("0 to 5", "5 to 15", "15 to 25", "25 to 35",
                 "35 to 45", "45 to 55", "55 to 65", "above 65"),
      ordered = TRUE
    )
  )

# Dimensions
site_levels <- area
S <- length(site_levels)

region_pairs <- case_data_area |>
  distinct(Area, Zone) |>
  arrange(match(Area, site_levels))
region_levels <- region_pairs |> pull(Zone)
R <- n_distinct(region_levels)

year_levels <- sort(unique(year_selected))
Y <- length(year_levels)
M <- 12L
T <- Y * M

block_of_year <- rep(1:n_periods, times = period_sizes)
P <- max(block_of_year)

# Mapping vectors
year_id <- rep(seq_len(Y), each = M)          # y(t)
month_id <- rep(seq_len(M), times = Y)         # m(t)
period_id <- block_of_year
region_id <- match(region_pairs$Zone, unique(region_levels))  # r(s)

# Monthly arrays
case_data_area <- case_data_area %>% mutate(Area = factor(Area, levels = site_levels))
calendar <- tidyr::expand_grid(Area = site_levels, year = year_levels, Month = 1:M) %>%
  left_join(case_data_area %>% select(Area, year, Month, case, death),
            by = c("Area", "year", "Month")) %>%
  mutate(case = coalesce(case, 0L), death = coalesce(death, 0L),
         s = as.integer(factor(Area, levels = site_levels)),
         y_idx = match(year, year_levels),
         m_idx = Month,
         t_idx = (y_idx - 1L) * M + m_idx)

# Cases
R_case_ts <- array(0L, dim = c(S, T))
R_case_ts[cbind(calendar$s, calendar$t_idx)] <- as.integer(calendar$case)
R_case_ts <- t(R_case_ts)

# Deaths (NEW for Model 9)
R_death_ts <- array(0L, dim = c(S, T))
R_death_ts[cbind(calendar$s, calendar$t_idx)] <- as.integer(calendar$death)
R_death_ts <- t(R_death_ts)

# Age arrays
age_levels <- levels(case_age_data_area$age_group_new)
A <- length(age_levels)
case_age_data_area <- case_age_data_area %>% mutate(Area = factor(Area, levels = site_levels))

age_comp <- case_age_data_area %>% ungroup() %>% mutate(year = as.integer(year)) %>%
  tidyr::complete(Area = site_levels, year = year_levels, age_group_new = age_levels,
                  fill = list(case = 0L, death = 0L)) %>%
  mutate(a = as.integer(factor(age_group_new, levels = age_levels)),
         y_idx = match(year, year_levels),
         s = as.integer(factor(Area, levels = site_levels)))

# Cases by age
R_case_yas <- array(0L, dim = c(A, Y, S))
R_case_yas[cbind(age_comp$a, age_comp$y_idx, age_comp$s)] <- as.integer(age_comp$case)
R_case_yas <- aperm(R_case_yas, c(2, 1, 3))

# Deaths by age (NEW for Model 9)
R_death_yas <- array(0L, dim = c(A, Y, S))
R_death_yas[cbind(age_comp$a, age_comp$y_idx, age_comp$s)] <- as.integer(age_comp$death)
R_death_yas <- aperm(R_death_yas, c(2, 1, 3))

# Create has_age_data_ys indicator: 1 if year-province has observed age-stratified case data
has_age_data_ys <- array(0L, dim = c(Y, S))
for (y_idx in seq_len(Y)) {
  for (s_idx in seq_len(S)) {
    year_val <- year_levels[y_idx]
    area_val <- site_levels[s_idx]
    # Check if this year-province has any observed age-stratified cases
    has_data <- case_age_data_area %>%
      filter(Area == area_val, year == year_val, case > 0) %>%
      nrow() > 0
    has_age_data_ys[y_idx, s_idx] <- as.integer(has_data)
  }
}

# Create has_death_data_ys indicator: 1 if year-province has death data (observed, even if 0 counts)
# Flag 0 only when data is completely missing (check raw data before "unknown" filtering)
has_death_data_ys <- array(0L, dim = c(Y, S))
for (y_idx in seq_len(Y)) {
  for (s_idx in seq_len(S)) {
    year_val <- year_levels[y_idx]
    area_val <- site_levels[s_idx]
    # Check if this year-province has any death data rows in raw data (before filtering "unknown")
    has_data <- case_age_data %>%
      filter(Area == area_val, year == year_val) %>%
      nrow() > 0
    has_death_data_ys[y_idx, s_idx] <- as.integer(has_data)
  }
}

# Population array
pop_candidates <- intersect(names(population_data_area), c("pop", "population", "Population", "N", "n", "count", "Total"))
if (length(pop_candidates) == 0) stop("No population column found.")
pop_col <- pop_candidates[1]

pop_tbl <- population_data_area %>% rename(Area = Province, pop = all_of(pop_col)) %>%
  mutate(Area = factor(Area, levels = site_levels), pop = as.integer(round(pop))) %>%
  filter(!is.na(pop)) %>%
  tidyr::complete(Area = site_levels, year = year_levels, age_group_new = age_levels,
                  fill = list(pop = 0L)) %>%
  mutate(a = as.integer(factor(age_group_new, levels = age_levels)),
         y_idx = match(year, year_levels),
         s = as.integer(factor(Area, levels = site_levels)))

N_ay <- array(0L, dim = c(A, Y, S))
N_ay[cbind(pop_tbl$a, pop_tbl$y_idx, pop_tbl$s)] <- as.integer(pop_tbl$pop)
N_ay <- aperm(N_ay, c(2, 1, 3))

N_tas <- N_ay[year_id, , , drop = FALSE]   # dims: [T, A, S]

# Fourier seasonality basis functions
H <- 3
cos_hm <- outer(1:H, 1:M, function(h, m) cos(2 * pi * m / M * h))  # [H, M]
sin_hm <- outer(1:H, 1:M, function(h, m) sin(2 * pi * m / M * h))  # [H, M]

# Pre-compute full Fourier matrices for all time points [T, H]
cos_t_h <- matrix(0.0, nrow = T, ncol = H)
sin_t_h <- matrix(0.0, nrow = T, ncol = H)
for(t_ind in 1:T) {
  m <- month_id[t_ind]
  theta <- 2 * pi * m / M
  for(h in 1:H) {
    cos_t_h[t_ind, h] <- cos(h * theta)
    sin_t_h[t_ind, h] <- sin(h * theta)
  }
}

# Sero data processing
sero_data_area <- sero_data %>% 
  filter(Province %in% area) %>%
  arrange(match(Province, site_levels))

omega <- 0.1882  # immune waning rate
nu <- 0.15       # clinical sign probability

# Reference age group: choose which age group has rho = 1
ref_age_label <- "25 to 35"
ref_age <- which(age_levels == ref_age_label)

# Prepare data for Stan (Model 9)
stan_data <- list(
  # Dimensions
  S = S,
  R = R,
  T = T,
  M = M,
  Y = Y,
  A = A,
  P = P,
  K_sero = length(c(2007, 2008, 2012)),
  
  # Mappings
  year_id = year_id,
  month_id = month_id,
  period_id = period_id,
  region_id = region_id,
  
  # Case data
  R_case_ts = R_case_ts,
  R_case_yas = R_case_yas,
  has_age_data_ys = has_age_data_ys,
  
  # Death data (NEW for Model 9)
  R_death_ts = R_death_ts,
  R_death_yas = R_death_yas,
  has_death_data_ys = has_death_data_ys,
  
  # Population
  N_tas = N_tas,
  
  # Seasonality (precomputed Fourier bases)
  H = H,
  cos_t_h = cos_t_h,
  sin_t_h = sin_t_h,
  
  # Reference age group (for normalization of rho_pa and CFR_pa)
  ref_age = ref_age,
  
  # Serology data
  sero_lower_s = sero_data_area$sero_lower,
  sero_upper_s = sero_data_area$sero_upper,
  sero_year_id = match(c(2007, 2008, 2012), year_levels),
  
  # Fixed parameters
  omega = omega,
  nu = nu,
  sigma_lambda_0 = 0.1,
  sigma_lambda_y = 0.122,
  sigma_lambda_init = 0.75
)

save(stan_data, calendar, site_levels, region_levels, age_levels, 
     file = here("Rout", "stan_data_09.RData"))
