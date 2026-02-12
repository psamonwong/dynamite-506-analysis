library(here)
library(rstan)

load(here("Rout", "fit_09.RData"))

# =============================================
# DISPERSION PARAMETERS
# =============================================
rstan::traceplot(fit, pars = "phi_case")
rstan::traceplot(fit, pars = "phi_death")

# =============================================
# BASELINE FORCE OF INFECTION
# =============================================
rstan::traceplot(fit, pars = "mu_log_lambda_0")
rstan::traceplot(fit, pars = "log_lambda_0_s")
rstan::traceplot(fit, pars = "sigma_lambda_init")

# Province-specific random walk innovations
rstan::traceplot(fit, pars = "epsilon_lambda_ys_init")

# =============================================
# SEASONALITY - BASELINE (province-specific, time-invariant)
# =============================================
rstan::traceplot(fit, pars = "a_cos_base")
rstan::traceplot(fit, pars = "a_sin_base")
rstan::traceplot(fit, pars = "tau_cos_base")
rstan::traceplot(fit, pars = "tau_sin_base")

# =============================================
# SEASONALITY - PERIOD DEVIATIONS
# =============================================
rstan::traceplot(fit, pars = "sigma_delta_cos")
rstan::traceplot(fit, pars = "sigma_delta_sin")

# Derived period-specific coefficients
rstan::traceplot(fit, pars = "a_cos")
rstan::traceplot(fit, pars = "a_sin")

# =============================================
# SENSITIVITY (region-year random walk)
# =============================================
rstan::traceplot(fit, pars = "logit_Se_init")
rstan::traceplot(fit, pars = "sigma_Se")

# Derived sensitivity
rstan::traceplot(fit, pars = "Se_y")

# =============================================
# AGE-SPECIFIC SUSCEPTIBILITY (Model 9)
# =============================================
rstan::traceplot(fit, pars = "log_rho_pa")
rstan::traceplot(fit, pars = "sigma_rho")

# Derived (normalized to ref_age = 1)
rstan::traceplot(fit, pars = "rho_pa")

# =============================================
# CASE-FATALITY RATE (Model 9 core)
# =============================================
rstan::traceplot(fit, pars = "mu_logit_CFR")
rstan::traceplot(fit, pars = "sigma_CFR")
rstan::traceplot(fit, pars = "logit_CFR_pa_raw")

# Derived CFR by period and age
rstan::traceplot(fit, pars = "CFR_pa")

# =============================================
# SEROPREVALENCE
# =============================================
rstan::traceplot(fit, pars = "H0_s")
rstan::traceplot(fit, pars = "p_sero_s")

# =============================================
# DIAGNOSTICS - KEY DERIVED QUANTITIES
# =============================================
# Baseline FOI by year-province
rstan::traceplot(fit, pars = "lambda_base_ys")

# Expected cases and deaths (subset for diagnostics)
# rstan::traceplot(fit, pars = "C_ts")  # Large - use sparingly
# rstan::traceplot(fit, pars = "D_ts")  # Large - use sparingly
