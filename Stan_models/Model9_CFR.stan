data {
    // Dimensions
    int<lower=1> T;                 // Number of months = Y * M
    int<lower=1> Y;                 // number of years
    int<lower=1> M;                 // number of months per year
    int<lower=1> P;                 // number of period blocks
    int<lower=1> S;                 // number of provinces
    int<lower=1> R;                 // number of regions
    int<lower=1> A;                 // number of age groups
    int<lower=1> K_sero;            // number of serology measurement years (2007, 2008, and 2012 = 3)
    
  // Mappings
  array[T] int<lower=1, upper=Y> year_id;         // y(t): year corresponding to month t
  array[T] int<lower=1, upper=M> month_id;        // m(t): month corresponding to month t
  array[Y] int<lower=1, upper=P> period_id;       // p(y(t)): period block corresponding to month t
  array[S] int<lower=1, upper=R> region_id;       // r(s): region corresponding to province s

  // Population size
  array[T, A, S] int<lower=0> N_tas;          // Population size at time t, age a, province s
  
  // Reported cases
  array[T, S] int<lower=0> R_case_ts;              // Monthly total cases (aggregated over age)
  array[Y, A, S] int<lower=0> R_case_yas;         // Annual age-stratified cases
  array[Y, S] int<lower=0, upper=1> has_age_data_ys;  // Indicator: 1 if age-stratified data available for year y, province s
  
  // Reported deaths
  array[T, S] int<lower=0> R_death_ts;             // Monthly total deaths (aggregated over age)
  array[Y, A, S] int<lower=0> R_death_yas;        // Annual age-stratified deaths
  array[Y, S] int<lower=0, upper=1> has_death_data_ys;  // Indicator: 1 if age-stratified death data available
  
  // Susceptible proportions
  int S_t;                            // Susceptible proportion (assumed constant)

  // Seasonality configurations
  int<lower=1> H;                 // number of Fourier harmonics
  matrix[T, H] cos_t_h;           // Precomputed cosine basis [T x H] for each time point
  matrix[T, H] sin_t_h;           // Precomputed sine basis [T x H] for each time point

  // Reference age group for susceptibility and CFR
  int<lower=1, upper=A> ref_age;  // Reference age group index (1 to A)

  // Serology data
  array[S] real<lower=0, upper=1> sero_lower_s;
  array[S] real<lower=0, upper=1> sero_upper_s;
  array[K_sero] int sero_year_id;

  // Immune waning rate
  real<lower=0> omega;

  // Clinical fraction
  real<lower=0, upper=1> nu;
  
  // Fixed hyperparameter for baseline FOI variation across provinces
  real<lower=0> sigma_lambda_0;
  
  // Fixed year-to-year FOI variation
  real<lower=0> sigma_lambda_y;
}

parameters{
  // === Dispersion parameters === //
  real<lower=0> phi_case;               // Dispersion for monthly case totals
  real<lower=0> phi_death;              // Dispersion for monthly death totals
  
  // === Test sensitivity and specificity parameters === //
  array[R] real logit_Se_init;
  array[R, Y-1] real logit_Se_innov;
  real<lower=0> sigma_Se;

  // === Force of infection parameters === //
  real mu_log_lambda_0;
  array[S] real log_lambda_0_s;
  array[S] real epsilon_lambda_ys_init;
  array[Y-1, S] real epsilon_lambda_ys_innov;

  // Seasonal variation parameters (baseline + period deviations)
  array[H, S] real a_cos_base;
  array[H, S] real a_sin_base;
  array[H, P, S] real delta_cos_raw;
  array[H, P, S] real delta_sin_raw;
  
  array[S] real<lower=0> tau_cos_base;
  array[S] real<lower=0> tau_sin_base;
  array[P] real<lower=0> sigma_delta_cos;
  array[P] real<lower=0> sigma_delta_sin;

  // === Age-specific susceptibility === //
  array[P, A] real log_rho_pa;            // log relative susceptibility by age group and period
  real<lower=0> sigma_rho;                // SD for age effects across periods

  // === Case-Fatality Rate (CFR) by period-age === //
  array[P, A] real logit_CFR_pa_raw;      // logit CFR deviations by period and age
  real<lower=0> sigma_CFR;                // SD for CFR effects across periods
  real mu_logit_CFR;                      // Mean logit-CFR (hyperprior)

  array[S] real<lower=0> H0_s;

  // === Hyperparameters === //
  real<lower=0> sigma_lambda_init;
}

transformed parameters {
  // === Force of infection components === //
  matrix[T, S] lambda_ts;
  matrix[Y, S] lambda_base_ys;
  matrix[Y, S] lambda_year_ys;
  
  // === Sensitivity and specificity === //
  matrix[R, Y] Se_y;
  vector[Y] Sp_y;
  
  // === Expected latent cases and reported cases === //
  array[T, A, S] real<lower=0> I_tas;
  matrix[T, S] I_ts;
  array[T, A, S] real<lower=0> C_tas;
  matrix[T, S] C_ts;
  
  // === Annual age-stratified cases === //
  array[Y, A, S] real<lower=0> I_yas;
  array[Y, A, S] real<lower=0> C_yas;
  
  // === Expected deaths === //
  array[T, A, S] real<lower=0> D_tas;
  matrix[T, S] D_ts;
  array[Y, A, S] real<lower=0> D_yas;
  
  // === Seroprevalence === //
  vector[S] p_sero_s;
  
  // === Age and period effects === //
  array[H, P, S] real a_cos;
  array[H, P, S] real a_sin;
  matrix[P, A] rho_pa;      // Age-specific relative susceptibility
  matrix[P, A] CFR_pa;      // Case-fatality rate by period and age
  
  // === Random walk processes === //
  matrix[Y, S] epsilon_lambda_ys;

  // Initialize specificity
  Sp_y = rep_vector(1.0, Y);

  // === Build period-specific seasonal coefficients === //
  for(h_ind in 1:H){
    for(p_ind in 1:P){
      for(s_ind in 1:S){
        a_cos[h_ind, p_ind, s_ind] = tau_cos_base[s_ind] * a_cos_base[h_ind, s_ind] + 
                         sigma_delta_cos[p_ind] * delta_cos_raw[h_ind, p_ind, s_ind];
        a_sin[h_ind, p_ind, s_ind] = tau_sin_base[s_ind] * a_sin_base[h_ind, s_ind] + 
                         sigma_delta_sin[p_ind] * delta_sin_raw[h_ind, p_ind, s_ind];
      }
    }
  }
  
  // === Build sensitivity random walk === //
  {
    matrix[R, Y] logit_Se_y;
    for(r_ind in 1:R){
      logit_Se_y[r_ind, 1] = logit_Se_init[r_ind];
      for (y_ind in 2:Y) {
        logit_Se_y[r_ind, y_ind] = logit_Se_y[r_ind, y_ind-1] + logit_Se_innov[r_ind, y_ind-1];
      }
    }
    Se_y = inv_logit(logit_Se_y);
  }
  
  // === Build province-specific random walks === //
  {
    epsilon_lambda_ys[1, ] = to_row_vector(epsilon_lambda_ys_init);
    for (y_ind in 2:Y) {
      epsilon_lambda_ys[y_ind, ] = epsilon_lambda_ys[y_ind-1, ] + to_row_vector(epsilon_lambda_ys_innov[y_ind-1, ]);
    }
  }

  // === Baseline force of infection === //
  for(y_ind in 1:Y){
    lambda_base_ys[y_ind, ] = exp(to_row_vector(log_lambda_0_s) + epsilon_lambda_ys[y_ind, ]);
  }

  // === Age-specific susceptibility (constraint: reference age = 1) === //
  for(p_ind in 1:P){
    for(a_ind in 1:A){
      rho_pa[p_ind, a_ind] = exp(log_rho_pa[p_ind, a_ind]);
    }
    // Normalize: set reference age group to rho = 1
    real ref_rho = rho_pa[p_ind, ref_age];
    for(a_ind in 1:A){
      rho_pa[p_ind, a_ind] = rho_pa[p_ind, a_ind] / ref_rho;
    }
  }

  // === Case-Fatality Rate (CFR) by period-age (constraint: reference age = 1) === //
  for(p_ind in 1:P){
    for(a_ind in 1:A){
      CFR_pa[p_ind, a_ind] = inv_logit(mu_logit_CFR + sigma_CFR * logit_CFR_pa_raw[p_ind, a_ind]);
    }
    // Normalize: set reference age group to CFR = 1 (proportional adjustment)
    real ref_CFR = CFR_pa[p_ind, ref_age];
    for(a_ind in 1:A){
      CFR_pa[p_ind, a_ind] = CFR_pa[p_ind, a_ind] / ref_CFR;
      CFR_pa[p_ind, a_ind] = fmin(CFR_pa[p_ind, a_ind], 0.99);  // Cap at 99%
    }
  }

  // === Force of infection by time and site (using precomputed Fourier basis) === //
  for(t_ind in 1:T){
    int y_ind = year_id[t_ind];
    int p_ind = period_id[y_ind];
    for(s_ind in 1:S){
      real f_seasonal = 0;
      for(h_ind in 1:H){
        f_seasonal += a_cos[h_ind, p_ind, s_ind] * cos_t_h[t_ind, h_ind] +
                      a_sin[h_ind, p_ind, s_ind] * sin_t_h[t_ind, h_ind];
      }
      lambda_ts[t_ind, s_ind] = lambda_base_ys[y_ind, s_ind] * exp(f_seasonal);
    }
  }

  // === Latent true clinical cases and deaths (age-stratified, monthly) === //
  // I_{t,a,s} = nu * N_{t,a,s} * (1 - exp(-lambda_{t,a,s}))
  // D_{t,a,s} = C_{t,a,s} * CFR_{p(t),a}
  for(t_ind in 1:T){
    int y_ind = year_id[t_ind];
    int p_ind = period_id[y_ind];
    for(s_ind in 1:S){
      int r_ind = region_id[s_ind];
      vector[A] lambda_ta = lambda_ts[t_ind, s_ind] * rho_pa[p_ind, ]';
      vector[A] I_a = nu * to_vector(N_tas[t_ind, , s_ind]) .* (1 - exp(-lambda_ta));
      vector[A] C_a = Se_y[r_ind, y_ind] * I_a;
      vector[A] D_a = C_a .* CFR_pa[p_ind, ]';
      
      for(a_ind in 1:A){
        I_tas[t_ind, a_ind, s_ind] = I_a[a_ind];
        C_tas[t_ind, a_ind, s_ind] = C_a[a_ind];
        D_tas[t_ind, a_ind, s_ind] = D_a[a_ind];
      }
      I_ts[t_ind, s_ind] = sum(I_a);
      C_ts[t_ind, s_ind] = sum(C_a);
      D_ts[t_ind, s_ind] = sum(D_a);
    }
  }

  // === Annual aggregations === //
  for(y_ind in 1:Y){
    for(a_ind in 1:A){
      for(s_ind in 1:S){
        int t_start = (y_ind - 1) * M + 1;
        int t_end = y_ind * M;
        I_yas[y_ind, a_ind, s_ind] = sum(I_tas[t_start:t_end, a_ind, s_ind]);
        C_yas[y_ind, a_ind, s_ind] = sum(C_tas[t_start:t_end, a_ind, s_ind]);
        D_yas[y_ind, a_ind, s_ind] = sum(D_tas[t_start:t_end, a_ind, s_ind]);
      }
    }
  }

  // === Annual aggregated force of infection (for seroprevalence) === //
  for(y_ind in 1:Y){
    int t_start = (y_ind - 1) * M + 1;
    int t_end = y_ind * M;
    for(s_ind in 1:S){
      lambda_year_ys[y_ind, s_ind] = mean(lambda_ts[t_start:t_end, s_ind]);
    }
  }

  // === Seroprevalence calculation with waning immunity === //
  for(s_ind in 1:S){
    real p_sum = 0;
    for (k_ind in 1:K_sero) {
      int y0 = sero_year_id[k_ind];
      real cumulative_hazard = H0_s[s_ind] * exp(-omega * y0);
      // Accumulate past infections with exponential waning
      for (y_ind in 1:y0) {
        real years_since = y0 - y_ind;
        cumulative_hazard += lambda_year_ys[y_ind, s_ind] * exp(-omega * years_since);
      }
      real p_y = 1 - exp(-cumulative_hazard);
      p_sum += p_y;
    }
    p_sero_s[s_ind] = p_sum / K_sero;  // Average across survey years
  }
}

model {
    // === Priors === //
    phi_case ~ normal(6, 1);
    phi_death ~ normal(3, 1);

    mu_log_lambda_0 ~ normal(-10, 2);
    log_lambda_0_s ~ normal(mu_log_lambda_0, sigma_lambda_0);
    H0_s ~ exponential(1);
    
    epsilon_lambda_ys_init ~ normal(0, sigma_lambda_init);
    sigma_lambda_init ~ exponential(10);

    for (y_ind in 1:(Y-1)) {
      epsilon_lambda_ys_innov[y_ind, ] ~ normal(0, sigma_lambda_y);
    }

    // Seasonal variation priors
    for(h_ind in 1:H){
      a_cos_base[h_ind, ] ~ normal(0, 1);
      a_sin_base[h_ind, ] ~ normal(0, 1);
    }
    
    tau_cos_base ~ normal(0, 0.5);
    tau_sin_base ~ normal(0, 0.5);
    
    sigma_delta_cos[1] ~ normal(0, 0.05);
    sigma_delta_sin[1] ~ normal(0, 0.05);
    for(p_ind in 2:P){
      sigma_delta_cos[p_ind] ~ normal(sigma_delta_cos[p_ind-1], 0.05);
      sigma_delta_sin[p_ind] ~ normal(sigma_delta_sin[p_ind-1], 0.05);
    }
    
    for(h_ind in 1:H){
      for(p_ind in 1:P){
        delta_cos_raw[h_ind, p_ind, ] ~ std_normal();
        delta_sin_raw[h_ind, p_ind, ] ~ std_normal();
      }
    }

    // Sensitivity priors
    logit_Se_init ~ normal(logit(0.4), 0.3);
    for(r_ind in 1:R){
      logit_Se_innov[r_ind, ] ~ normal(0, sigma_Se);
    }
    sigma_Se ~ exponential(1);

    // Age effect priors
    sigma_rho ~ exponential(1);
    for(p_ind in 1:P){
      log_rho_pa[p_ind, ] ~ normal(0, sigma_rho);
    }

    // CFR priors (similar to age susceptibility)
    mu_logit_CFR ~ normal(logit(0.05), 1);  // Prior around 5% baseline CFR
    sigma_CFR ~ exponential(1);
    for(p_ind in 1:P){
      logit_CFR_pa_raw[p_ind, ] ~ std_normal();
    }

    // === Likelihoods === //
    // Monthly case totals: negative binomial (vectorized per row)
    for(t_ind in 1:T)
      R_case_ts[t_ind, ] ~ neg_binomial_2(C_ts[t_ind, ], phi_case);

    // Monthly death totals: negative binomial (vectorized per row)
    for(t_ind in 1:T)
      R_death_ts[t_ind, ] ~ neg_binomial_2(D_ts[t_ind, ], phi_death);

    // Annual age-stratified cases: multinomial (only when data available)
    for(y_ind in 1:Y){
      for(s_ind in 1:S){
        if(has_age_data_ys[y_ind, s_ind]){
          vector[A] q_case_yas;
          int N_case_ys = 0;
          for(a_ind in 1:A){
            q_case_yas[a_ind] = C_yas[y_ind, a_ind, s_ind];
            N_case_ys += R_case_yas[y_ind, a_ind, s_ind];
          }
          if(N_case_ys > 0){
            q_case_yas = q_case_yas / sum(q_case_yas);  // Normalize to proportions
            R_case_yas[y_ind, , s_ind] ~ multinomial(q_case_yas);
          }
        }
      }
    }

    // Annual age-stratified deaths: multinomial (only when data available)
    for(y_ind in 1:Y){
      for(s_ind in 1:S){
        if(has_death_data_ys[y_ind, s_ind]){
          vector[A] q_death_yas;
          int N_death_ys = 0;
          for(a_ind in 1:A){
            q_death_yas[a_ind] = D_yas[y_ind, a_ind, s_ind];
            N_death_ys += R_death_yas[y_ind, a_ind, s_ind];
          }
          if(N_death_ys > 0){
            q_death_yas = q_death_yas / sum(q_death_yas);  // Normalize to proportions
            R_death_yas[y_ind, , s_ind] ~ multinomial(q_death_yas);
          }
        }
      }
    }

    // Seroprevalence likelihood: interval-censored normal
    for(s_ind in 1:S){
      real sigma_sero = 0.01;
      target += log_diff_exp(
        std_normal_lcdf((sero_upper_s[s_ind] - p_sero_s[s_ind]) / sigma_sero),
        std_normal_lcdf((sero_lower_s[s_ind] - p_sero_s[s_ind]) / sigma_sero)
      );
    }
}

generated quantities {
  // === Reconstruct seasonal multipliers for posterior inference === //
  array[M, P, S] real seasonal_multiplier;
  for(m_ind in 1:M){
    for(p_ind in 1:P){
      for(s_ind in 1:S){
        real f_seasonal = 0;
        real theta_m = 2.0 * pi() * m_ind / M;
        for(h_ind in 1:H){
          f_seasonal += a_cos[h_ind, p_ind, s_ind] * cos(h_ind * theta_m) +
                        a_sin[h_ind, p_ind, s_ind] * sin(h_ind * theta_m);
        }
        seasonal_multiplier[m_ind, p_ind, s_ind] = exp(f_seasonal);
      }
    }
  }
}
