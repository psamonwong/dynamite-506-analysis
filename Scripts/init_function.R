init_fun <- function(stan_data) {
  S <- stan_data$S
  Y <- stan_data$Y
  R <- stan_data$R
  H <- stan_data$H
  P <- stan_data$P

  list(
    # Dispersion
    phi_case = abs(rnorm(1, 6, 1)),

    # Baseline FOI
    mu_log_lambda_0 = rnorm(1, -10, 2),
    sigma_lambda_0  = abs(rnorm(1, 0.1, 0.05)),
    log_lambda_0_s  = rnorm(S, -10, 0.1),

    # Province random walk
    epsilon_lambda_ys_init  = rnorm(S, 0, 0.02),
    epsilon_lambda_ys_innov = matrix(rnorm((Y-1)*S, 0, 0.05), Y-1, S),
    sigma_lambda_y     = abs(rnorm(1, 0.05, 0.01)),
    sigma_lambda_init  = rexp(1, 10),

    # Sensitivity random walk
    logit_Se_init  = rnorm(R, qlogis(0.4), 0.3),
    logit_Se_innov = matrix(rnorm(R*(Y-1), 0, 0.1), R, Y-1),
    sigma_Se = rexp(1, 1),

    # Seasonality
    a_cos_base = matrix(rnorm(H*S, 0, 1), H, S),
    a_sin_base = matrix(rnorm(H*S, 0, 1), H, S),

    tau_cos_base = abs(rnorm(S, 0, 0.5)),
    tau_sin_base = abs(rnorm(S, 0, 0.5)),

    sigma_delta_cos = abs(rnorm(P, 0, 0.05)),
    sigma_delta_sin = abs(rnorm(P, 0, 0.05)),

    delta_cos_raw = array(rnorm(H*P*S, 0, 1), c(H,P,S)),
    delta_sin_raw = array(rnorm(H*P*S, 0, 1), c(H,P,S)),

    # Background hazard
    H0_s = rexp(S, 1)
  )
}
