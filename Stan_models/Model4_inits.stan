data {
    // Dimensions
    int<lower=1> T;                 // Number of months = Y * M
    int<lower=1> Y;                 // number of years
    int<lower=1> M;                 // number of months per year
    int<lower=1> P;                 // number of period blocks
    int<lower=1> S;                 // number of provinces
    int<lower=1> R;                 // number of regions
    int<lower=1> K_sero;            // number of serology measurement years (2007, 2008, and 2012 = 3)
    

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
  array[K_sero] int sero_year_id;               // year indices for serology data

  // Immune waning rate
  real<lower=0> omega;               // immune waning rate

  // Clinical fraction
  real<lower=0, upper=1> nu;        // clinical fraction
}

parameters{
  // === Dispersion parameters === //
  real<lower=0> phi_case;               // Dispersion parameter for reported cases
  
  // === Test sensitivity and specificity parameters (yearly constants) === //
  real logit_Se_init;                   // Initial sensitivity (year 1) on logit scale (global)
  array[Y-1] real logit_Se_innov;       // Year-to-year innovations for sensitivity (random walk)

  // === Force of infection parameters === //

  // Baseline force of infection parameters
  real log_lambda_0;                      // overall baseline log force of infection
  array[S] real epsilon_lambda_ys_init;   // initial year effect for each province (year 1)
  array[Y-1, S] real epsilon_lambda_ys_innov; // year-to-year innovations for each province (random walk)

  // Seasonal variation parameters (baseline + year deviations)
  array[H, S] real a_cos_base;            // province baseline seasonality (time-invariant)
  array[H, S] real a_sin_base;            // province baseline seasonality (time-invariant)
  array[H, Y, S] real delta_cos_raw;      // year-specific deviations from baseline (raw)
  array[H, Y, S] real delta_sin_raw;      // year-specific deviations from baseline (raw)
  
  // Hierarchical hyperparameters for seasonality
  array[S] real<lower=0> tau_cos_base;    // amplitude of baseline seasonality for each province
  array[S] real<lower=0> tau_sin_base;    // amplitude of baseline seasonality for each province
  array[Y] real<lower=0> sigma_delta_cos; // how much year y deviates from baseline
  array[Y] real<lower=0> sigma_delta_sin; // how much year y deviates from baseline

  // === Hyperparameters === //
  real<lower=0> sigma_lambda_y;          // standard deviation of year-to-year innovations (random walk)
  real<lower=0> sigma_lambda_init;      // standard deviation of initial year effects
  real<lower=0> sigma_Se;                // standard deviation of sensitivity innovations (random walk)
  
}

transformed parameters {
  // === Parameters declarations === //
  array[T, S] real<lower=0> C_ts;        // C_t,s expected observed cases at time t, province s
  array[T, S] real<lower=0> I_ts;        // I_t,s true incident cases at time t, province s
  array[T, S] real lambda_ts;   // force of infection at time t, province s
  array[Y, S] real lambda_base_ys;   // baseline log force of infection at year y, province s
  array[Y, M, S] real<lower=0> kappa_seasonal_yms;   // seasonal variation component of log force of infection at year y, province s, month m
  
  array[Y] real<lower=0, upper=1> Sp_y;        // Test specificity for year y
  array[Y] real<lower=0, upper=1> Se_y;        // Test sensitivity for year y (random walk)
  array[Y, S] real epsilon_lambda_ys;          // province-specific year effects (random walk)
  
  // Year-specific seasonal coefficients (baseline + deviation)
  array[H, Y, S] real a_cos;
  array[H, Y, S] real a_sin;

  for(i in 1:Y) Sp_y[i] = 1; // initialize specificity

   // === Seroprevalence calculation === //
  array[S] real<lower=0, upper=1> p_sero_s;   // seroprevalence in each province
  array[S] real lambda_star_s;              // average force of infection for seroprevalence calculation during serology year
  
  // Build year-specific seasonal coefficients: baseline + year deviation
  for(h in 1:H){
    for(y in 1:Y){
      for(s in 1:S){
        a_cos[h, y, s] = tau_cos_base[s] * a_cos_base[h, s] + 
                         sigma_delta_cos[y] * delta_cos_raw[h, y, s];
        a_sin[h, y, s] = tau_sin_base[s] * a_sin_base[h, s] + 
                         sigma_delta_sin[y] * delta_sin_raw[h, y, s];
      }
    }
  }
  
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
  
  // Build province-specific random walks
  for (s_ind in 1:S) {
    epsilon_lambda_ys[1, s_ind] = epsilon_lambda_ys_init[s_ind];
    for (y in 2:Y) {
      epsilon_lambda_ys[y, s_ind] = epsilon_lambda_ys[y-1, s_ind] + epsilon_lambda_ys_innov[y-1, s_ind];
    }
  }

  // === Baseline force of infection calculation === //
  for(y_ind in 1:Y)
    for(s_ind in 1:S)
      lambda_base_ys[y_ind, s_ind] = exp(log_lambda_0 + epsilon_lambda_ys[y_ind, s_ind]);     // Baseline log force of infection at year y, province s 

  // === Seasonal variation calculation with Fourier with H harmonics (year-province) === //
  for(y_ind in 1:Y){
    for(s_ind in 1:S){
      vector[M] log_kappa_yms; // log seasonal modifier for month m, year y, province s
        for(m_ind in 1:M){
            real theta = 2 * pi() * m_ind / M; // angle for month m
            real f = 0;
            for(h in 1:H){
              // Fourier series expansion for month m with harmonics h
              f += a_cos[h, y_ind, s_ind] * cos(h * theta) +
                   a_sin[h, y_ind, s_ind] * sin(h * theta);
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
      I_ts[t_ind, s_ind] = nu * N_tas[t_ind, s_ind] * S_t * (1 - exp(-lambda_ts[t_ind, s_ind]));      // True incident cases at time t, province s

  // === Observed cases calculation === //
  for(t_ind in 1:T)
    for(s_ind in 1:S)
      C_ts[t_ind, s_ind] = (I_ts[t_ind, s_ind] * Se_y[year_id[t_ind]]) + ((1 - Sp_y[year_id[t_ind]]) * N_tas[t_ind, s_ind]);      // Expected observed cases at time t, province s 


  // === Seroprevalence calculation === //
  for(s_ind in 1:S){
    real p_sum = 0;
    for (k in 1:K_sero) {
      int y0 = sero_year_id[k];                 // serology measurement year
      real Hazards = 0;                                // cumulative annual hazard up to y0
      for (y_ind in 1:y0) {
        real years_since = y0 - y_ind;
        real lambda_year = 12 * lambda_base_ys[y_ind, s_ind];   // year^-1
        Hazards += lambda_year * exp(-omega * years_since);       // waning discount
      }
      real p_y = 1 - exp(-Hazards);                    // prevalence for year y0
      p_sum += p_y;
    }
    p_sero_s[s_ind] = p_sum / K_sero;            // average across 2007–2012
  }
}

model {
    // === Priors === //
    // Dispersion parameters
    phi_case ~ normal(6, 1);               // Prior for dispersion parameter for reported cases) 
    // Baseline force of infection parameters
    
    log_lambda_0 ~ normal(-6, 2);    // Prior for overall baseline log force of infection
    epsilon_lambda_ys_init ~ normal(0, sigma_lambda_init);   // Prior for initial year effects (each province)
    sigma_lambda_init ~ exponential(10); // Prior for initial year effects SD (mean 0.02)

    for (y in 1:(Y-1)) {
      epsilon_lambda_ys_innov[y] ~ normal(0, sigma_lambda_y);  // Prior for year-to-year innovations (each province)
    }

    // Seasonal variation parameters (baseline + year-specific deviations)
    // Province baseline seasonality (time-invariant)
    for(h in 1:H){
      for(s in 1:S){
        a_cos_base[h, s] ~ normal(0, 1);
        a_sin_base[h, s] ~ normal(0, 1);
      }
    }
    
    // Baseline amplitude (province-specific)
    for(s in 1:S){
      tau_cos_base[s] ~ normal(0, 0.5);
      tau_sin_base[s] ~ normal(0, 0.5);
    }
    
    // Year-specific deviation scales with temporal smoothing
    sigma_delta_cos[1] ~ normal(0, 0.05);  // first year: small deviation from baseline
    sigma_delta_sin[1] ~ normal(0, 0.05);
    for(y in 2:Y){
      sigma_delta_cos[y] ~ normal(sigma_delta_cos[y-1], 0.05);  // smooth evolution
      sigma_delta_sin[y] ~ normal(sigma_delta_sin[y-1], 0.05);
    }
    
    // Raw year deviations (standard normal)
    for(h in 1:H){
      for(y in 1:Y){
        for(s in 1:S){
          delta_cos_raw[h, y, s] ~ std_normal();
          delta_sin_raw[h, y, s] ~ std_normal();
        }
      }
    }

    // Sensitivity and Specificity parameters (yearly)
    logit_Se_init ~ normal(logit(0.5), 0.3);           // Prior for initial sensitivity (centered at 0.5)
    logit_Se_innov ~ normal(0, sigma_Se);              // Prior for year-to-year innovations (allows gradual changes and shocks)
    sigma_Se ~ exponential(1);                         // Prior for sensitivity innovation SD (allows flexibility for shocks)

    // Hyperparameters
    sigma_lambda_y ~ normal(0.05, 0.01); // Prior for standard deviation of year-to-year innovations (constrained for gradual changes)

    // Immune waning adjustment parameter
   // log_delta ~ normal(0, 1); // Prior for immune waning adjustment parameter

    // === Likelihood  for reported cases === //
    for(t_ind in 1:T)
      for(s_ind in 1:S)
      R_case_ts[t_ind, s_ind] ~ neg_binomial_2(C_ts[t_ind, s_ind] + 1e-9, phi_case);      // Reported cases likelihood (ensure strictly positive mean)

    // === Likelihood for seroprevalence data (interval-censored)=== //
    for(s_ind in 1:S){
      real sigma_sero = 0.01;  // measurement uncertainty
      target += log_diff_exp(
        std_normal_lcdf((sero_upper_s[s_ind] - p_sero_s[s_ind]) / sigma_sero),
        std_normal_lcdf((sero_lower_s[s_ind] - p_sero_s[s_ind]) / sigma_sero)
      );
      }
}
