library(here)
library(tidyverse)
library(rstan)

source('Scripts/init_function.R')

## RUN MODEL 9
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

##############################
print(sprintf("Started at %s", Sys.time()))
##############################
# LOAD PREPARED DATA
load(here("Rout", "stan_data_09.RData"))
stan_data$S_t <- 1

##############################
# MODEL FITTING
iter <- getOption("stan.iter", 200)
warmup <- getOption("stan.warmup", 100)
chains <- getOption("stan.chains", 4)
thin <- getOption("stan.thin", 1)
seed <- getOption("stan.seed", 124L)

message(sprintf("Sampling: iter=%d warmup=%d chains=%d thin=%d", iter, warmup, chains, thin))

mod <- stan_model(here("Stan_models", "Model9_CFR.stan"), verbose = FALSE)

fit <- rstan::sampling(
  object = mod, 
  data = stan_data, 
  iter = iter, 
  warmup = warmup, 
  cores = 4,
  chains = chains, 
  thin = thin, 
  seed = seed, 
  refresh = max(1, iter %/% 20),  # print progress every 1%
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  init = function() init_fun(stan_data),
  pars = c(
    # Dispersion
    "phi_case",
    "phi_death",
    
    # Detection sensitivity
    "Se_y",
    "Sp_y",
    "sigma_Se",
    "logit_Se_init",
    "logit_Se_innov",
    
    # Force of infection - baseline
    "mu_log_lambda_0",
    "log_lambda_0_s",
    "sigma_lambda_init",
    "epsilon_lambda_ys",
    "lambda_base_ys",
    
    # Seasonality
    "a_cos",
    "a_sin",
    "a_cos_base",
    "a_sin_base",
    "tau_cos_base",
    "tau_sin_base",
    "sigma_delta_cos",
    "sigma_delta_sin",
    
    # Age-specific susceptibility
    "log_rho_pa",
    "rho_pa",
    "sigma_rho",
    
    # Case-Fatality Rate (NEW for Model 9)
    "mu_logit_CFR",
    "logit_CFR_pa_raw",
    "sigma_CFR",
    "CFR_pa",
    
    # Expected cases
    "C_ts",
    "C_tas",
    "I_tas",
    
    # Expected deaths (NEW for Model 9)
    "D_ts",
    "D_tas",
    "D_yas",
    
    # Seroprevalence
    "p_sero_s",
    "H0_s",
    "lambda_year_ys"
  ),
  include = TRUE
)

# SAVE OUTPUT
save(fit, stan_data, calendar, site_levels, region_levels, age_levels,
     file = here("Rout", "fit_09.RData"))

###############################
message(sprintf("Finished at %s", Sys.time()))

## Summary statistics
#print(fit, probs = c(0.025, 0.5, 0.975), digits = 3)

## clear workspace
rm(list = ls())
