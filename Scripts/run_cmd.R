library(here)
library(tidyverse)
library(rstan)

## RUN MODEL
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

##############################
print(sprintf("Started at %s", Sys.time()))
##############################
# LOAD PREPARED DATA
load(here("Rout", "stan_data_02.RData"))
##############################
# MODEL FITTING
iter <- getOption("stan.iter", 400)
warmup <- getOption("stan.warmup", 100)
chains <- getOption("stan.chains", 4)
thin <- getOption("stan.thin", 1)
seed <- getOption("stan.seed", 124L)
message(sprintf("Sampling: iter=%d warmup=%d chains=%d thin=%d", iter, warmup, chains, thin))
mod <- stan_model(here("Stan_models", "Model2_simple.stan"), verbose = F)
fit <- rstan::sampling(object = mod, data = stan_data, iter = iter, warmup = warmup, cores = 4,
                       chains = chains, thin = thin, seed = seed, 
                       refresh = max(1, iter * 0.05),  # print progress every 5%
                       control = list(adapt_delta = 0.95, max_treedepth = 12),
                       pars = c(
                         # Fixed effects
                         "lambda_0", "alpha_0_tau", "alpha_0_se",
                         "beta_0", "pi_fp",
                         
                         # Seasonality baselines
                          "A_pr", "B_pr",
                         
                         # Age effects
                         "beta_a_raw", "sigma_age",
                         
                         # Scale parameters (hierarchical SDs)
                         "sigma_tau_pr", "sigma_tau_ys",
                         "sigma_lambda_pr", "sigma_lambda_ys",
                        "sigma_A_ys", "sigma_B_ys",
                         "sigma_age_pr", #"sigma_age_ys",
                         
                         # Overdispersion
                         "phi_c",
                         
                         # Transformed parameters to keep
                         "tau_ys", "log_lambda_base",
                         "mu_case_ts", "C_yas", "q_case_ys"
                       ),
                       include = TRUE)
# SAVE OUTPUT
save(fit, stan_data,
     file = here("Rout", "fit_02.RData"))
###############################
message(sprintf("Finished at %s", Sys.time()))
## clear workspace
rm(list = ls())

