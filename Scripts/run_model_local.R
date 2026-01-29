library(here)
library(tidyverse)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#############################
# DATA PREPARATION (same as sandbox3.R except p_y removed)
#############################

# Load datasets
case_data <- read_csv(here("Analysis_data", "Scrub_typhus", "analysis_data.csv"))
case_age_data <- read_csv(here("Analysis_data", "Scrub_typhus", "case_death_by_age.csv"))
population_data <- read_csv(here("Analysis_data", "Population", "population_data_summarised.csv"))

# Setup
area <- 
  #case_data %>% distinct(Area) %>% pull(Area)
   c( "Chiang Rai", "Lampang",  "Lamphun", "Songkhla", "Phatthalung", "Narathiwat",
    "Uttaradit", "Tak", "Phitsanulok")
period_size <- 5
n_periods <- 4
start_year <- 2004
year_selected <- seq(from = start_year, to = start_year + period_size * n_periods - 1)

# Filter case data
case_data_area <- case_data %>%
  filter(Area %in% area, year %in% year_selected) %>%
  arrange(Area, year, Month)

# Modal age per Area-year (exclude unknown)
modal_age <- case_age_data %>%
  filter(Area %in% area, year %in% year_selected, age_group_new != "unknown") %>%
  group_by(Area, year, age_group_new) %>%
  summarize(total_case = sum(case, na.rm = TRUE), .groups = "drop") %>%
  group_by(Area, year) %>%
  slice_max(total_case, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Area, year, modal_age = age_group_new)

# Impute unknown ages + aggregate
case_age_data_area <- case_age_data %>%
  filter(Area %in% area, year %in% year_selected) %>%
  left_join(modal_age, by = c("Area", "year")) %>%
  mutate(age_group_new = if_else(age_group_new == "unknown" & !is.na(modal_age), modal_age, age_group_new)) %>%
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

block_of_year <- 1L + (year_levels - min(year_levels)) %/% period_size
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

R_case_ts <- array(0L, dim = c(S, T))
R_death_ts <- array(0L, dim = c(S, T))
R_case_ts[cbind(calendar$s, calendar$t_idx)] <- as.integer(calendar$case)
R_death_ts[cbind(calendar$s, calendar$t_idx)] <- as.integer(calendar$death)

R_case_ts <- t(R_case_ts)
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

R_case_yas <- array(0L, dim = c(A, Y, S))
R_death_yas <- array(0L, dim = c(A, Y, S))
R_case_yas[cbind(age_comp$a, age_comp$y_idx, age_comp$s)] <- as.integer(age_comp$case)
R_death_yas[cbind(age_comp$a, age_comp$y_idx, age_comp$s)] <- as.integer(age_comp$death)

R_case_yas <- aperm(R_case_yas, c(2, 1, 3))
R_death_yas <- aperm(R_death_yas, c(2, 1, 3))

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

# Susceptibility
S_tas <- array(1.0, dim = c(T, A, S))

H <- 3



# Prepare data for Stan
stan_data <- list(
  S = S,
  R = R,
  T = T,
  M = M,
  Y = Y,
  A = A,
  P = P,
  R_case_ts = R_case_ts,
  R_death_ts = R_death_ts,
  R_case_yas = R_case_yas,
  R_death_yas = R_death_yas,
  N_tas = N_tas,
  S_tas = S_tas,
  year_id = year_id,
  month_id = month_id,
  period_id = period_id,
  region_id = region_id,
  H = H
)


save(stan_data, calendar, site_levels, region_levels, year_levels, age_levels,
     file = here("Rout", "stan_data_01.RData"))
###############################
system(paste('Rscript --vanilla Scripts/run_cmd.R'))
##############################


# Sampling settings
iter <- getOption("stan.iter", 2000)
warmup <- getOption("stan.warmup", 1000)
chains <- getOption("stan.chains", 4)
thin <- getOption("stan.thin", 4)
seed <- getOption("stan.seed", 124L)
message(sprintf("Sampling: iter=%d warmup=%d chains=%d thin=%d", iter, warmup, chains, thin))

mod <- stan_model(here("Stan_models", "Model1_fourier.stan"), verbose = T)
fit <- rstan::sampling(object = mod, data = stan_data, iter = iter, warmup = warmup, cores = 4,
                       chains = chains, thin = thin, seed = seed, 
                       control = list(adapt_delta = 0.90, max_treedepth = 12),
                       pars = c(
                         # Hyperparameters
                         "phi_c", "phi_d",
                         "sigma_lambda_s", "sigma_epsilon_lambda",
                         "sigma_alpha_Se", "sigma_delta_Se",
                         "sigma_alpha_FP", "sigma_delta_FP",
                         "sigma_gamma_kappa", "sigma_delta_kappa",
                         # "sigma_gamma_CFR",  # fixed in transformed data
                         
                         # Key parameters (high-level, interpretable)
                         "nu",
                         "alpha_lambda_s",
                         "beta_kappa_a", "h_kappa_a", "gamma_kappa_pr_base", "gamma_kappa_pr",
                         "beta_CFR_a", "h_CFR_a", "gamma_CFR_pr_base", "gamma_CFR_pr",
                         "alpha_Se_yr_base", "alpha_Se_yr", "delta_Se_s",
                         "alpha_FP_yr_base", "alpha_FP_yr", "delta_FP_s",
                         "a_cos", "a_sin",
                         
                         # Transformed parameters (derived quantities)
                         "kappa_age_yas", "CFR_yas",
                         "mu_case_ts", "mu_death_ts",
                         "q_case_ys", "q_death_ys",
                         
                         # Generated quantities
                         "log_lik_case", "log_lik_death", "log_lik_case_age", "log_lik_death_age"
                       ),
                       include = TRUE)

save(fit, stan_data,
     file = here("Rout", "fit_01.RData"))


load('Rout/fit_01.RData')


# Kappa (age-transmission modifier)
rstan::traceplot(fit, pars = 'beta_kappa_a')
rstan::traceplot(fit, pars = 'h_kappa_a')
rstan::traceplot(fit, pars = 'gamma_kappa_pr_base')
rstan::traceplot(fit, pars = 'sigma_gamma_kappa')
rstan::traceplot(fit, pars = 'sigma_delta_kappa')
rstan::traceplot(fit, pars = 'delta_kappa_ys')


# CFR (case fatality rate)
rstan::traceplot(fit, pars = 'beta_CFR_a')
rstan::traceplot(fit, pars = 'h_CFR_a')
rstan::traceplot(fit, pars = 'gamma_CFR_pr_base')
#rstan::traceplot(fit, pars = 'sigma_gamma_CFR')

# Case detection sensitivity (Se)
rstan::traceplot(fit, pars = 'alpha_Se_yr_base')
rstan::traceplot(fit, pars = 'sigma_alpha_Se')
rstan::traceplot(fit, pars = 'sigma_delta_Se')
rstan::traceplot(fit, pars = 'delta_Se_s')


# False positives (FP)
rstan::traceplot(fit, pars = 'alpha_FP_yr_base')
rstan::traceplot(fit, pars = 'sigma_alpha_FP')
rstan::traceplot(fit, pars = 'sigma_delta_FP')
rstan::traceplot(fit, pars = 'delta_FP_s')

# FOI (force of infection)
rstan::traceplot(fit, pars = 'sigma_lambda_s') 
rstan::traceplot(fit, pars = 'sigma_epsilon_lambda')
#rstan::traceplot(fit, pars = 'rho_lambda')  # AR(1) coefficient

# Clinical fraction
rstan::traceplot(fit, pars = 'nu')

# NB2 dispersion
rstan::traceplot(fit, pars = 'phi_c')
rstan::traceplot(fit, pars = 'phi_d')


names(fit)
