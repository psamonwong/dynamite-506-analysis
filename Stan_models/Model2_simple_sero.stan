data {
    // Dimensions
    int<lower=1> T;                 // Number of months = Y * M
    int<lower=1> Y;                 // number of years
    int<lower=1> M;                 // number of months per year
    int<lower=1> P;                 // number of period blocks
    int<lower=1> S;                 // number of provinces
    int<lower=1> R;                 // number of regions
    int<lower=1> Y_sero;            // number of years with serology data

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

  // Serology data
  array[S] real<lower=0, upper=1> sero_lower_s;   // lower bounds of seroprevalence in each province
  array[S] real<lower=0, upper=1> sero_upper_s;   // upper bounds of seroprevalence in each province
  array[Y_sero] int sero_year_id;               // year indices for serology data
  // Immune waning rate
  real<lower=0> omega;               // immune waning rate
}

parameters{
  // === Dispersion parameters === //
  real<lower=0> phi_case;               // Dispersion parameter for reported cases
  
  // === Test sensitivity and specificity parameters (yearly constants) === //
  real logit_Se_init;                   // Initial sensitivity (year 1) on logit scale
  array[Y-1] real logit_Se_innov;       // Year-to-year innovations for sensitivity (random walk)

  // === Force of infection parameters === //

  // Baseline force of infection parameters
  real log_lambda_0;                      // overall baseline log force of infection
  array[S] real epsilon_lambda_s;         // site-level effect (constant across years)

  // Seasonal variation parameters
  array[H, S] real a_cos;   // Fourier cosine coefficients for seasonal variation for harmonic h, province s
  array[H, S] real a_sin;   // Fourier sine coefficients for seasonal variation for harmonic h, province s

  // === Hyperparameters === //
  real<lower=0> sigma_lambda_s;          // standard deviation of site-level effects
  real<lower=0> sigma_Se;                // standard deviation of sensitivity innovations (random walk)

}

transformed parameters {
  // === Parameters declarations === //
  array[T, S] real<lower=0> C_ts;        // C_t,s expected observed cases at time t, province s
  array[T, S] real<lower=0> I_ts;        // I_t,s true incident cases at time t, province s
  array[T, S] real lambda_ts;   // force of infection at time t, province s
  array[S] real lambda_base_s;   // baseline force of infection at site s (constant across years)
  array[Y, M, S] real<lower=0> kappa_seasonal_yms;   // seasonal variation component of log force of infection at year y, province s, month m
  
  array[Y] real<lower=0, upper=1> Sp_y;        // Test specificity for year y
  array[Y] real<lower=0, upper=1> Se_y;        // Test sensitivity for year y (random walk)
  // removed province-specific year effects; lambda varies only seasonally

  for(i in 1:Y) Sp_y[i] = 1; // initialize specificity

   // === Seroprevalence calculation === //
  array[S] real<lower=0, upper=1> p_sero_s;   // seroprevalence in each province
  array[S] real lambda_star_s;              // average force of infection for seroprevalence calculation during serology year
  
  // Build sensitivity random walk on logit scale
  {
    array[Y] real logit_Se_y;
    logit_Se_y[1] = logit_Se_init;
    for (y in 2:Y) {
      logit_Se_y[y] = logit_Se_y[y-1] + logit_Se_innov[y-1];
    }
    for (y in 1:Y) {
      Se_y[y] = inv_logit(logit_Se_y[y]);
    }
  }
  
  // === Baseline force of infection calculation (site-constant across years) === //
  for(s_ind in 1:S)
    lambda_base_s[s_ind] = exp(log_lambda_0 + epsilon_lambda_s[s_ind]);

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
      lambda_ts[t_ind, s_ind] =  lambda_base_s[s_ind] *       // Baseline force of infection (constant across years)
                                 kappa_seasonal_yms[year_id[t_ind],  month_id[t_ind], s_ind];     // Seasonal variation component of log force of infection

  // === True incident cases calculation === //
  for(t_ind in 1:T)
    for(s_ind in 1:S)
      I_ts[t_ind, s_ind] = N_tas[t_ind, s_ind] * S_t * (1 - exp(-lambda_ts[t_ind, s_ind]));      // True incident cases at time t, province s

  // === Observed cases calculation === //
  for(t_ind in 1:T)
    for(s_ind in 1:S)
      C_ts[t_ind, s_ind] = (I_ts[t_ind, s_ind] * Se_y[year_id[t_ind]]) + ((1 - Sp_y[year_id[t_ind]]) * N_tas[t_ind, s_ind]);      // Expected observed cases at time t, province s 


  // === Seroprevalence calculation === //
  for(s_ind in 1:S){
    //real accum_lambda = 0;
    //accum_lambda = 12 * Y_sero * lambda_base_s[s_ind]; // sum baseline over serology years (constant across years)
    lambda_star_s[s_ind] = 12 * lambda_base_s[s_ind]; // average baseline force of infection (equals 12 * lambda_base_s)
    p_sero_s[s_ind] = lambda_star_s[s_ind] / (lambda_star_s[s_ind] + omega); // seroprevalence calculation for province s with immune waning rate omega
  }
}

model {
    // === Priors === //
    // Dispersion parameters
    phi_case ~ normal(6, 1);               // Prior for dispersion parameter for reported cases) 
    // Baseline force of infection parameters
    log_lambda_0 ~ normal(-15, 2);               // Prior for overall baseline log force of infection
    epsilon_lambda_s ~ normal(0, sigma_lambda_s);   // Prior for site-level effects (each province)

    // Seasonal variation parameters
    for(h in 1:H){
      a_cos[h] ~ normal(0, 1);        // Prior for Fourier cosine coefficients for seasonal variation
      a_sin[h] ~ normal(0, 1);        // Prior for Fourier sine coefficients for seasonal variation
    }

    // Sensitivity and Specificity parameters (yearly)
    logit_Se_init ~ normal(logit(0.8), 0.2);           // Prior for initial sensitivity (centered at 0.8)
    logit_Se_innov ~ normal(0, sigma_Se);            // Prior for year-to-year innovations (allows gradual changes and shocks)
    sigma_Se ~ exponential(1);                       // Prior for sensitivity innovation SD (allows flexibility for shocks)}

    // Hyperparameters
    sigma_lambda_s ~ normal(0.1, 0.05); // Prior for standard deviation of site-level effects

    // === Likelihood  for reported cases === //
    for(t_ind in 1:T)
      for(s_ind in 1:S)
      R_case_ts[t_ind, s_ind] ~ neg_binomial_2(C_ts[t_ind, s_ind] + 1e-9, phi_case);      // Reported cases likelihood (ensure strictly positive mean)

    // === Likelihood for seroprevalence data (interval-censored)=== //
    for(s_ind in 1:S){
      real sigma_sero = 0.02;  // measurement uncertainty
      target += log_diff_exp(
        std_normal_lcdf((sero_upper_s[s_ind] - p_sero_s[s_ind]) / sigma_sero),
        std_normal_lcdf((sero_lower_s[s_ind] - p_sero_s[s_ind]) / sigma_sero)
      );
      }
}
