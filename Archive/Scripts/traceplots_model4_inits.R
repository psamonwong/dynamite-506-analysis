library(here)
library(rstan)

load(here("Rout", "fit_04.RData"))

# Dispersion
rstan::traceplot(fit, pars = 'phi_case')

# Baseline FOI
rstan::traceplot(fit, pars = 'log_lambda_0')
rstan::traceplot(fit, pars = 'sigma_lambda_init')
rstan::traceplot(fit, pars = 'sigma_lambda_y')

# Seasonality baseline (province-specific, time-invariant)
rstan::traceplot(fit, pars = 'a_cos_base')
rstan::traceplot(fit, pars = 'a_sin_base')
rstan::traceplot(fit, pars = 'tau_cos_base')
rstan::traceplot(fit, pars = 'tau_sin_base')

# Year-specific deviations from baseline
rstan::traceplot(fit, pars = 'sigma_delta_cos')
rstan::traceplot(fit, pars = 'sigma_delta_sin')

# Seasonality (year-province specific, derived)
rstan::traceplot(fit, pars = 'a_cos')
rstan::traceplot(fit, pars = 'a_sin')

# Sensitivity (global temporal trend)
rstan::traceplot(fit, pars = 'logit_Se_init')
rstan::traceplot(fit, pars = 'logit_Se_innov')
rstan::traceplot(fit, pars = 'sigma_Se')

# Seroprevalence fit
rstan::traceplot(fit, pars = 'p_sero_s')
