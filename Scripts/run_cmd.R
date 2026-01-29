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
load(here("Rout", "stan_data_01.RData"))
##############################
# MODEL FITTING
iter <- getOption("stan.iter", 250)
warmup <- getOption("stan.warmup", 100)
chains <- getOption("stan.chains", 4)
thin <- getOption("stan.thin", 1)
seed <- getOption("stan.seed", 124L)
message(sprintf("Sampling: iter=%d warmup=%d chains=%d thin=%d", iter, warmup, chains, thin))
mod <- stan_model(here("Stan_models", "Model1_fourier.stan"), verbose = F)
fit <- rstan::sampling(object = mod, data = stan_data, iter = iter, warmup = warmup, cores = 4,
                       chains = chains, thin = thin, seed = seed, 
                       refresh = max(1, iter * 0.05),  # print progress every 5%
                       control = list(adapt_delta = 0.95, max_treedepth = 12),
                       pars = c(
                         # Overdispersion
                         "phi_c", "phi_d",
                         
                         # Baseline FOI (with regional hierarchy)
                         "alpha_lambda_r", "alpha_lambda_s",
                         "sigma_lambda_r", "sigma_lambda_s",
                         
                         # Seasonality (Fourier)
                         "a_cos", "a_sin", "sigma_seasonality_pr",
                         
                         # Age transmission (sum-to-zero age effects)
                         "alpha_kappa_a", "alpha_kappa_age",
                         "gamma_kappa_pr_base", "gamma_kappa_pr",
                         "sigma_gamma_kappa",
                         
                         # Case detection sensitivity (Se) - additive decomposition
                         "alpha_Se_r", "delta_Se_s", "delta_Se_y",
                         "sigma_delta_Se", "sigma_delta_Se_y",
                         
                         # False positives (FP) - additive decomposition
                         "alpha_FP_r", "delta_FP_s", "delta_FP_y",
                         "sigma_delta_FP", "sigma_delta_FP_y",
                         
                         # Case fatality rate (sum-to-zero age effects)
                         "alpha_CFR_a", "alpha_CFR_age",
                         "gamma_CFR_pr_base", "gamma_CFR_pr",
                         "sigma_gamma_CFR",
                         
                         # Transformed parameters
                         "log_lambda_base", "kappa_month_yms", "kappa_age_yas",
                         "CFR_yas", "mu_case_ts", "mu_death_ts",
                         "q_case_ys", "q_death_ys",
                         
                         # Log-likelihoods
                         "log_lik_case", "log_lik_death",
                         "log_lik_case_age", "log_lik_death_age"
                         
                         # ===== MODEL 2 PARAMETERS (COMMENTED OUT) =====
                         # # Fixed effects
                         # "lambda_0", "alpha_0_tau", "alpha_0_se", "pi_fp",
                         # 
                         # # Seasonality baselines (Fourier)
                         # "A_pr", "B_pr",
                         # 
                         # # Age effects
                         # "beta_a_raw", "beta_a", "sigma_age",
                         # 
                         # # Hierarchical SDs
                         # "sigma_tau_pr", "sigma_tau_ys",
                         # "sigma_lambda_pr", "sigma_lambda_ys",
                         # "sigma_A_ys", "sigma_B_ys",
                         # "sigma_age_pr",
                         # 
                         # # Overdispersion
                         # "phi_c",
                         # 
                         # # Scaled random effects
                         # "epsilon_tau_pr", "epsilon_tau_ys",
                         # "epsilon_lambda_pr", "epsilon_lambda_ys",
                         # "epsilon_A_ys", "epsilon_B_ys",
                         # "epsilon_age_pr",
                         # 
                         # # Transformed parameters
                         # "tau_ys", "log_lambda_base",
                         # "kappa_month_yms", "kappa_age_yas",
                         # "mu_case_ts", "C_yas", "q_case_ys"
                       ),
                       include = TRUE)
# SAVE OUTPUT
save(fit, stan_data,
     file = here("Rout", "fit_01.RData"))
###############################
message(sprintf("Finished at %s", Sys.time()))
## clear workspace
rm(list = ls())

