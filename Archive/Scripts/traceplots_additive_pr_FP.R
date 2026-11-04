library(here)
library(rstan)

load(here("Rout", "fit_additive_pr_FP.RData"))

# Dispersion
rstan::traceplot(fit, pars = 'phi_c')
rstan::traceplot(fit, pars = 'phi_d')

# Clinical fraction
rstan::traceplot(fit, pars = 'nu')

# FOI additive baseline
rstan::traceplot(fit, pars = 'alpha_lambda')
rstan::traceplot(fit, pars = 'delta_lambda_s')
rstan::traceplot(fit, pars = 'delta_lambda_y')
rstan::traceplot(fit, pars = 'sigma_lambda_noise')

# Seasonality (Fourier)
rstan::traceplot(fit, pars = 'a_cos')
rstan::traceplot(fit, pars = 'a_sin')
rstan::traceplot(fit, pars = 'delta_cos')
rstan::traceplot(fit, pars = 'delta_sin')

# Age transmission (Lee–Carter)
rstan::traceplot(fit, pars = 'beta_kappa_a')
rstan::traceplot(fit, pars = 'h_kappa_a')
rstan::traceplot(fit, pars = 'gamma_kappa_pr_base')
rstan::traceplot(fit, pars = 'gamma_kappa_pr')
rstan::traceplot(fit, pars = 'sigma_gamma_kappa')
rstan::traceplot(fit, pars = 'sigma_delta_kappa')
rstan::traceplot(fit, pars = 'delta_kappa_ys')

# CFR
rstan::traceplot(fit, pars = 'beta_CFR_a')
rstan::traceplot(fit, pars = 'h_CFR_a')
rstan::traceplot(fit, pars = 'gamma_CFR_pr_base')
rstan::traceplot(fit, pars = 'gamma_CFR_pr')

# Se
rstan::traceplot(fit, pars = 'alpha_Se_yr_base')
rstan::traceplot(fit, pars = 'alpha_Se_yr')
rstan::traceplot(fit, pars = 'delta_Se_s')
rstan::traceplot(fit, pars = 'sigma_alpha_Se')
rstan::traceplot(fit, pars = 'sigma_delta_Se')

# FP (alpha_pr + delta_s + delta_y + noise)
rstan::traceplot(fit, pars = 'alpha_FP_pr')
rstan::traceplot(fit, pars = 'delta_FP_s')
rstan::traceplot(fit, pars = 'delta_FP_y')
rstan::traceplot(fit, pars = 'sigma_delta_FP')
rstan::traceplot(fit, pars = 'sigma_FP_noise')
