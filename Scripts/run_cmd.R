library(here)
library(tidyverse)
library(rstan)
library(tictoc)

## RUN MODEL
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

tic()
##############################
print(sprintf("Started at %s", Sys.time()))
##############################
# LOAD PREPARED DATA
load(here("Rout", "stan_data_01.RData"))
##############################
# MODEL FITTING
iter <- getOption("stan.iter", 400)
warmup <- getOption("stan.warmup", 100)
chains <- getOption("stan.chains", 4)
thin <- getOption("stan.thin", 4)
seed <- getOption("stan.seed", 124L)
message(sprintf("Sampling: iter=%d warmup=%d chains=%d thin=%d", iter, warmup, chains, thin))
mod <- stan_model(here("Stan_models", "Model_additive_simple.stan"), verbose = F)
print("Stan model compiled.")
fit <- rstan::sampling(object = mod, data = stan_data, iter = iter, warmup = warmup, cores = 4,
                       chains = chains, thin = thin, seed = seed, 
                       refresh = max(1, iter / 20),  # print progress every 5%
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
                         "alpha_lambda_s", "delta_kappa_ys",
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
# SAVE OUTPUT
save(fit, stan_data, calendar, site_levels, region_levels, year_levels, age_levels,
     file = here("Rout", "fit_02.RData"))
###############################
message(sprintf("Finished at %s", Sys.time()))
## clear workspace
rm(list = ls())
toc()
