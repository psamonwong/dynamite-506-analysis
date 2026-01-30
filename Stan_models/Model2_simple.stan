data {
    // Dimensions
    int<lower=1> T;                 // Number of months = Y * M
    int<lower=1> Y;                 // number of years
    int<lower=1> M;                 // number of months per year
    int<lower=1> P;                 // number of period blocks
    int<lower=1> S;                 // number of provinces
    int<lower=1> R;                 // number of regions

  // Mappings
  array[T] int<lower=1, upper=Y> year_id;         // y(t): year corresponding to month t
  array[T] int<lower=1, upper=M> month_id;        // m(t): month corresponding to month t
  array[Y] int<lower=1, upper=P> period_id;       // p(y(t)): period block corresponding to month t
  array[S] int<lower=1, upper=R> region_id;       // r(s): region corresponding to province s

  // Population size
  array[T,S] int<lower=0> N_tas;          // N_t,s population size at time t, province s
  // Susceptible proportions
  int S_t;                            //  S_t susceptible proportion at time t (assumed constant across provinces and time)
  // Reported cases and deaths
  array[T, S] int<lower=0> R_case_ts;          // R_case_t,s reported cases at time t, province s

  // Seasonality configurations
  int<lower=1> H;                 // number of Fourier harmonics
}

parameters{
  // === Dispersion parameters === //
  real<lower=0> phi_case;               // Dispersion parameter for reported cases
  // === Test sensitivity and specificity parameters (yearly constants) === //
  array[Y] real<lower=0, upper=1> Se_y;        // Test sensitivity for year y
  array[Y] real<lower=0, upper=1> Sp_y;        // Test specificity for year y

  // === Force of infection parameters === //

  // Baseline force of infection parameters
  real log_lambda_0;                      // overall baseline log force of infection
  array[S] real epsilon_lambda_s;         // province-specific random effect on baseline log force of infection
  array[Y] real epsilon_lambda_y;         // year-specific random effect on baseline log force of infection

  // Seasonal variation parameters
  array[H, S] real a_cos;   // Fourier cosine coefficients for seasonal variation for harmonic h, province s
  array[H, S] real a_sin;   // Fourier sine coefficients for seasonal variation for harmonic h, province s

  // === Hyperparameters === //
  real<lower=0> sigma_lambda_s;          // standard deviation of province-specific random effect on baseline log force of infection
  real<lower=0> sigma_lambda_y;          // standard deviation of year-specific random effect on baseline log force of infection

  // === Testing intensity (for false positives) === //
  real alpha_test;                        // intercept on test-intensity logit
  real beta_test;                         // slope for log(lambda_ts) effect on testing intensity
  real<lower=0> tau_test_base;            // baseline tested fraction multiplier
  array[S] real epsilon_test_s;           // province random effect on testing intensity
  array[Y] real epsilon_test_y;           // year random effect on testing intensity
  real<lower=0> sigma_test_s;             // sd for province testing RE
  real<lower=0> sigma_test_y;             // sd for year testing RE

}

transformed parameters {
  // === Parameters declarations === //
  array[T, S] real<lower=0> C_ts;        // C_t,s expected observed cases at time t, province s
  array[T, S] real<lower=0> I_ts;        // I_t,s true incident cases at time t, province s
  array[T, S] real lambda_ts;   // log force of infection at time t, province s
  array[Y, S] real lambda_base_ys;   // baseline log force of infection at year y, province s
  array[Y, M, S] real<lower=0> kappa_seasonal_yms;   // seasonal variation component of log force of infection at year y, province s, month m
  array[T, S] real<lower=0> N_eff_ts;   // effective number tested driving false positives


  // === Baseline force of infection calculation === //
  for(y_ind in 1:Y)
    for(s_ind in 1:S)
      lambda_base_ys[y_ind, s_ind] = exp(log_lambda_0 + epsilon_lambda_s[s_ind] + epsilon_lambda_y[y_ind]);     // Baseline log force of infection at year y, province s 

  // === Seasonal variation calculation with Fourier with H harmonics === //
  for(y_ind in 1:Y){
    for(s_ind in 1:S){
      vector[M] log_kappa_yms; // log seasonal modifier for month m, year y, province s
        for(m_ind in 1:M){
            real theta = 2 * pi() * m_ind / M; // angle for month m
            real f = 0;
            for(h in 1:H){
              // Fourier series expansion for month m with harmonics h
              f += a_cos[h, s_ind] * cos(h * theta) +
                   a_sin[h, s_ind] * sin(h * theta);
            }
            log_kappa_yms[m_ind] = f;
        }
        // Normalize across months to mean 1 (identifiability)
        real sum_exp_m = 0; 
        for(m_ind in 1:M) sum_exp_m += exp(log_kappa_yms[m_ind]); // sum over months 
        for(m_ind in 1:M) kappa_seasonal_yms[y_ind, m_ind, s_ind] = M * exp(log_kappa_yms[m_ind])/ sum_exp_m; // normalize to mean 1
      }
  }

  // == Force of infection calculation === //
  for(t_ind in 1:T)
    for(s_ind in 1:S)
      lambda_ts[t_ind, s_ind] =  lambda_base_ys[year_id[t_ind] , s_ind] *       // Baseline log force of infection
                                 kappa_seasonal_yms[year_id[t_ind],  month_id[t_ind], s_ind];     // Seasonal variation component of log force of infection

  // === True incident cases calculation === //
  for(t_ind in 1:T)
    for(s_ind in 1:S)
      I_ts[t_ind, s_ind] = N_tas[t_ind, s_ind] * S_t * (1 - exp(-lambda_ts[t_ind, s_ind]));      // True incident cases at time t, province s

  // === Effective testing volume (drives false positives) === //
  for(t_ind in 1:T){
    int y = year_id[t_ind];
    for(s_ind in 1:S){
      real logit_f = alpha_test + epsilon_test_s[s_ind] + epsilon_test_y[y]
                      + beta_test * log(lambda_ts[t_ind, s_ind] + 1e-9);
      real f = inv_logit(logit_f); // fraction in (0,1)
      N_eff_ts[t_ind, s_ind] = N_tas[t_ind, s_ind] * tau_test_base * f;
    }
  }

  // === Observed cases calculation === //
  for(t_ind in 1:T)
    for(s_ind in 1:S)
      C_ts[t_ind, s_ind] = (I_ts[t_ind, s_ind] * Se_y[year_id[t_ind]]) + ((1 - Sp_y[year_id[t_ind]]) * N_eff_ts[t_ind, s_ind]);      // Expected observed cases at time t, province s 

}

model {
    // === Priors === //
    // Dispersion parameters
    phi_case ~ normal(0, 5);               // Prior for dispersion parameter for reported cases

    // Baseline force of infection parameters
    log_lambda_0 ~ normal(-20, 5);               // Prior for overall baseline log force of infection
    epsilon_lambda_s ~ normal(0, sigma_lambda_s);         // Prior for province-specific random effect on baseline log force of infection
    epsilon_lambda_y ~ normal(0, sigma_lambda_y);         // Prior for year-specific random effect on baseline log

    // Seasonal variation parameters
    for(h in 1:H){
      a_cos[h] ~ normal(0, 1);        // Prior for Fourier cosine coefficients for seasonal variation
      a_sin[h] ~ normal(0, 1);        // Prior for Fourier sine coefficients for seasonal variation
    }

    // Sensitivity and Specificity parameters (yearly)
    for(y_ind in 1:Y){
      Se_y[y_ind] ~ beta(8, 2);        // Prior for test sensitivity at year y (mean 0.8)
      Sp_y[y_ind] ~ beta(9, 1);        // Prior for test specificity at year y (mean 0.9)
    }

    // Hyperparameters
    sigma_lambda_s ~ exponential(1);          // Prior for standard deviation of province-specific random effect on baseline log force of infection
    sigma_lambda_y ~ exponential(1);          // Prior for standard deviation of year-specific random effect on baseline log force of infection

    // Testing intensity priors
    alpha_test ~ normal(-2.94, 1);            // centers f around ~0.05 before other effects
    beta_test ~ normal(0, 1);
    tau_test_base ~ lognormal(log(0.01), 0.5); // small baseline tested fraction
    sigma_test_s ~ exponential(1);
    sigma_test_y ~ exponential(1);
    epsilon_test_s ~ normal(0, sigma_test_s);
    epsilon_test_y ~ normal(0, sigma_test_y);


    // === Likelihood  === // 
    for(t_ind in 1:T)
      for(s_ind in 1:S)
      R_case_ts[t_ind, s_ind] ~ neg_binomial_2(C_ts[t_ind, s_ind] + 1e-9, phi_case);      // Reported cases likelihood (ensure strictly positive mean)
}
