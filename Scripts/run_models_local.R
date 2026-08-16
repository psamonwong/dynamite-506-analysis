#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
################################
library(here)
library(tidyverse)
library(rstan)

source('Scripts/init_function.R')
################################
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
##############################
print(sprintf("Started at %s", Sys.time()))
##############################
# LOAD PREPARED DATA
dataset <- args[1]
seed <- args[2]

load(here("Rout", paste0("stan_data_", dataset, ".RData")))
load(here("Rout", paste0("model_settings_", dataset, ".RData")))

message(sprintf("Model: %s", model_file))
message(sprintf("Dataset: %s", dataset))
message(sprintf("Sampling: iter=%d warmup=%d chains=%d thin=%d", 
                model_settings$iter,
                model_settings$warmup,
                model_settings$chains,
                model_settings$thin))

mod <- stan_model(here("Stan_models", model_file), verbose = FALSE)

fit <- rstan::sampling(
  object = mod, 
  data = stan_data, 
  iter = model_settings$iter, 
  warmup = model_settings$warmup, 
  cores = 4,
  chains = model_settings$chains, 
  thin = model_settings$thin, 
  seed = seed, 
  refresh = max(1, model_settings$iter * 0.01), 
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  init = function() init_fun(stan_data),
  pars = c(
    # === Dispersion parameters === //
    "phi_case",
    "phi_death",
    
    # === Detection sensitivity (essential for interpretation) === //
    "Se_y",
    "sigma_Se",
    
    # === Force of infection === //
    "mu_log_lambda_0",
    "log_lambda_0_s",
 #   "sigma_lambda_init",
    "epsilon_lambda_ys",
    "lambda_base_ys",
    "lambda_year_ys",
    
    # === Seasonality === //
    "a_cos",
    "a_sin",
    "tau_cos_base",
    "tau_sin_base",
    "sigma_delta_cos",
    "sigma_delta_sin",
    
    # === Age-specific susceptibility === //
    "rho_pa",
    "sigma_rho",
    
    # === Case-Fatality Rate === //
    "mu_logit_CFR",
    "sigma_CFR",
    "CFR_pa",
    
    # === Posterior predictive checks === //
    "C_ts",
    "D_ts",
    
    # === Seroprevalence === //
    "p_sero_s",
    "H0_s",
    
    # === Log-likelihood (for LOO-IC/WAIC) === //
    "log_lik_case",
    "log_lik_death",
    "log_lik_case_age",
    "log_lik_death_age",
    
    # === Posterior predictive replicates === //
    "R_case_ts_rep",
    "R_death_ts_rep",
    
    # === Seasonality multipliers === //
    "seasonal_multiplier",
    
    # === Diagnostics === //
    "sum_log_lik"
  ),
  include = TRUE
)

# SAVE OUTPUT
save(fit, stan_data, calendar, site_levels, region_levels, age_levels,
     file = here("Rout", paste0("fit_", dataset, ".RData")))

###############################
message(sprintf("Finished at %s", Sys.time()))

rm(list = ls())
