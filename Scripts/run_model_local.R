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
sero_data <- read_csv(here("Analysis_data", "Gonwong_sero_data.csv"))

# Setup
area <- #case_data %>% distinct(Area) %>% pull(Area)
c( "Chiang Rai", "Lampang",  "Nakhon Ratchasima", "Surin", "Songkhla", "Phatthalung", "Uttaradit", "Tak")
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
cos_hm <- outer(1:H, 1:M, function(h, m) cos(2 * pi * m / M * h))  # [H, M]
sin_hm <- outer(1:H, 1:M, function(h, m) sin(2 * pi * m / M * h))  # [H, M]

# Sero data processing (not used in model but saved for reference)
sero_data_area <- sero_data %>% 
  filter(Province %in% area) %>%
  arrange(Province)

omega <- 1/2  # assumed immune waning rate

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
  H = H, cos_hm = cos_hm, sin_hm = sin_hm,
  sero_lower_s = sero_data_area$sero_lower,
  sero_upper_s = sero_data_area$sero_upper,
  omega = omega,
  Y_sero = nrow(sero_data_area),
  sero_year_id = match(c(2007, 2008, 2012), year_levels)
)


save(stan_data,
     file = here("Rout", "stan_data_01.RData"))


#####################################
# PARAMETER ESTIMATION VIA OPTIMIZATION
#####################################
library(here)
library(tidyverse)
library(rstan)
load(here("Rout", "stan_data_01.RData"))

mod <- stan_model(here("Stan_models", "Model1_fourier.stan"), verbose = FALSE)

fits <- replicate(20,
 optimizing(mod, data = stan_data, iter = 10000, verbose = TRUE, 
                      hessian = FALSE, as_vector = FALSE),
  simplify = FALSE
)

values <- sapply(fits, function(f) f$par$sigma_delta_Se)
values %>% hist()



opt_fit$par[names(opt_fit$par) %in% c("sigma_delta_Se", "sigma_delta_Se_y", "alpha_Se_r", "beta_FP",
                                        "sigma_delta_FP", "sigma_delta_FP_y", "alpha_FP_r")]

save(opt_fit, stan_data, file = here("Rout", "opt_fit_01.RData"))
message(sprintf("Optimization complete. Value: %.2f", opt_fit$value))


#####################################
system(paste('Rscript --vanilla Scripts/run_cmd.R'))
##############################

load('Rout/fit_01.RData')




# "lambda_0", "alpha_0_tau", "gamma", 'beta_0, nu, B0, beta_a_raw

# sigma_tau_pr sigma_tau_ys  sigma_lambda_pr sigma_lambda_ys
# sigma_seasonality_pr sigma_age_pr


library(rstan)

# Model 1 traceplots - High-level interpretable parameters only
plots <- list(
  dispersions = rstan::traceplot(
    fit,
    pars = c("phi_c", "phi_d"),
    inc_warmup = FALSE
  ),
  
  foi_hierarchy = rstan::traceplot(
    fit,
    pars = c("alpha_lambda_r", "sigma_lambda_r", "sigma_lambda_s"),
    inc_warmup = FALSE
  ),
  
  seasonality_baselines = rstan::traceplot(
    fit,
    pars = c("a_cos"),
    inc_warmup = FALSE
  ),

    seasonality_baselines2 = rstan::traceplot(
    fit,
    pars = c("a_sin"),
    inc_warmup = FALSE
  ),

      seasonality_baselines2 = rstan::traceplot(
    fit,
    pars = c("a_sin"),
    inc_warmup = FALSE
  ),
  
  age_kappa_structure1 = rstan::traceplot(
    fit,
    pars = c("alpha_kappa_a"),
    inc_warmup = FALSE
  ),
    age_kappa_structure2 = rstan::traceplot(
    fit,
    pars = c("alpha_kappa_age"),
    inc_warmup = FALSE
  ),
    age_kappa_structure3 = rstan::traceplot(
    fit,
    pars = c("gamma_kappa_pr", "sigma_gamma_kappa"),
    inc_warmup = FALSE
  ),
  
  se_fp_hyperparameters = rstan::traceplot(
    fit,
    pars = c("sigma_delta_Se", "sigma_delta_Se_y", "sigma_delta_FP", "sigma_delta_FP_y"),
    inc_warmup = FALSE
  ),
  
  cfr_structure = rstan::traceplot(
    fit,
    pars = c("alpha_CFR_a"),
    inc_warmup = FALSE
  ),

    cfr_structure2 = rstan::traceplot(
    fit,
    pars = c("alpha_CFR_age"),
    inc_warmup = FALSE
  ),
      cfr_structure2 = rstan::traceplot(
    fit,
    pars = c("gamma_CFR_pr"),
    inc_warmup = FALSE
  )
)

plots


# plots <- list(
#   fixed = rstan::traceplot(
#     fit,
#     pars = c("lambda_0", "alpha_0_tau", "alpha_0_se", "pi_fp"),
#     inc_warmup = FALSE
#   ),
#   seasonality_baselines = rstan::traceplot(
#     fit,
#     pars = c("A_pr", "B_pr"),
#     inc_warmup = FALSE
#   ),
#   age_effects = rstan::traceplot(
#     fit,
#     pars = c("beta_a_raw", "beta_a", "sigma_age"),
#     inc_warmup = FALSE
#   ),
#   scales = rstan::traceplot(
#     fit,
#     pars = c(
#       "sigma_tau_pr", "sigma_tau_ys",
#       "sigma_lambda_pr", "sigma_lambda_ys",
#       "sigma_A_ys", "sigma_B_ys",
#       "sigma_age_pr"
#     ),
#     inc_warmup = FALSE
#   ),
#   overdispersion = rstan::traceplot(
#     fit,
#     pars = "phi_c",
#     inc_warmup = FALSE
#   ),
#   random_effects_foi = rstan::traceplot(
#     fit,
#     pars = c("epsilon_lambda_pr", "epsilon_lambda_ys"),
#     inc_warmup = FALSE
#   ),
#   random_effects_seasonality = rstan::traceplot(
#     fit,
#     pars = c("epsilon_A_ys", "epsilon_B_ys"),
#     inc_warmup = FALSE
#   ),
#   random_effects_age = rstan::traceplot(
#     fit,
#     pars = "epsilon_age_pr",
#     inc_warmup = FALSE
#   ),
#   derived_quantities = rstan::traceplot(
#     fit,
#     pars = c("log_lambda_base", "kappa_month_yms", "kappa_age_yas"),
#     inc_warmup = FALSE
#   )
# )

# plots
