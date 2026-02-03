library(here)
library(rstan)

load(here("Rout", "fit_04.RData"))

# Dispersion
rstan::traceplot(fit, pars = 'phi_case')

# Baseline FOI
rstan::traceplot(fit, pars = 'log_lambda_0')
rstan::traceplot(fit, pars = 'sigma_lambda_init')
rstan::traceplot(fit, pars = 'sigma_lambda_y')

# Seasonality (Fourier coefficients - province-specific)
rstan::traceplot(fit, pars = 'a_cos')
rstan::traceplot(fit, pars = 'a_sin')

# Sensitivity (global temporal trend)
rstan::traceplot(fit, pars = 'logit_Se_init')
rstan::traceplot(fit, pars = 'logit_Se_innov')
rstan::traceplot(fit, pars = 'sigma_Se')

# Seroprevalence fit
rstan::traceplot(fit, pars = 'p_sero_s')
