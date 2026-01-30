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
stan_data$S_t <- 1
stan_data$N_tas <- stan_data$N_tas  %>% apply(c(1,3), sum)
##############################
# MODEL FITTING
iter <- getOption("stan.iter", 250)
warmup <- getOption("stan.warmup", 100)
chains <- getOption("stan.chains", 4)
thin <- getOption("stan.thin", 1)
seed <- getOption("stan.seed", 124L)
message(sprintf("Sampling: iter=%d warmup=%d chains=%d thin=%d", iter, warmup, chains, thin))
mod <- stan_model(here("Stan_models", "Model2_simple_sero.stan"), verbose = F)
fit <- rstan::sampling(object = mod, data = stan_data, iter = iter, warmup = warmup, cores = 4,
                       chains = chains, thin = thin, seed = seed, 
                       refresh = max(1, iter * 0.05),  # print progress every 5%
                                                        control = list(adapt_delta = 0.95, max_treedepth = 12),
                                                        pars = c(
                                                             "phi_case",
                                                             "Se_y",
                                                             "Sp_y",
                                                             "C_ts",
                                                             "log_lambda_0",
                                                             "epsilon_lambda_ys",
                                                             "lambda_base_ys",
                                                             "a_cos",
                                                             "a_sin",
                                                             "sigma_lambda_y",
                                                             "sigma_Se"
                                                        ),
                                                        include = TRUE)
# SAVE OUTPUT
save(fit, stan_data,
     file = here("Rout", "fit_02.RData"))
###############################
message(sprintf("Finished at %s", Sys.time()))
## clear workspace
rm(list = ls())

