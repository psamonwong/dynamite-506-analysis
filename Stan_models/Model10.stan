data {
    // === DIMENSIONS ===
    int<lower=1> T;                 // Total months (264 = 22 years × 12 months)
    int<lower=1> Y;                 // Years (22: 2003-2024)
    int<lower=1> M;                 // Months per year (12)
    int<lower=1> P;                 // Period blocks (4)
    int<lower=1> S;                 // Provinces (77)
    int<lower=1> R;                 // Regions (13)
    int<lower=1> A;                 // Age groups (8)
    int<lower=1> K_sero;            // Serology years (3: 2007, 2008, 2012)
    
  // === MAPPINGS ===
  array[T] int<lower=1, upper=Y> year_id;         // Year for each month
  array[T] int<lower=1, upper=M> month_id;        // Month (1-12) for each time
  array[Y] int<lower=1, upper=P> period_id;       // Period block (1-4) per year
  array[S] int<lower=1, upper=R> region_id;       // Region per province

  // === POPULATION ===
  array[T, A, S] int<lower=0> N_tas;              // Population by time, age, province
  
  // === OBSERVED CASES ===
  array[T, S] int<lower=0> R_case_ts;             // Monthly total cases
  array[Y, A, S] int<lower=0> R_case_yas;         // Annual age-stratified cases
  array[Y, S] int<lower=0, upper=1> has_age_data_ys;  // 1 if age case data available
  
  // === OBSERVED DEATHS (Model 9 core) ===
  array[T, S] int<lower=0> R_death_ts;            // Monthly total deaths
  array[Y, A, S] int<lower=0> R_death_yas;        // Annual age-stratified deaths
  array[Y, S] int<lower=0, upper=1> has_death_data_ys;  // 1 if age death data available

  // === SEASONALITY ===
  int<lower=1> H;                 // Fourier harmonics (3)
  matrix[T, H] cos_t_h;           // Cosine basis [T x H]
  matrix[T, H] sin_t_h;           // Sine basis [T x H]

  // === REFERENCE ===
  int<lower=1, upper=A> ref_age;  // Reference age for relative susceptibility

  // === SEROLOGY DATA ===
  array[S] real<lower=0, upper=1> sero_lower_s;  // Seroprevalence lower bound
  array[S] real<lower=0, upper=1> sero_upper_s;  // Seroprevalence upper bound
  array[K_sero] int sero_year_id;                 // Years with serology surveys

  // === FIXED PARAMETERS ===
  real<lower=0> omega;             // Immune waning rate
  real<lower=0, upper=1> nu;       // Clinical manifestation probability
  real<lower=0> sigma_lambda_0;    // Baseline FOI variation (0.1 fixed)
  real<lower=0> sigma_lambda_y;    // Year-to-year FOI variation (0.122 fixed)
  real<lower=0> sigma_lambda_init; // Initial FOI variation (0.75 fixed)

  // === SEROLOGY AGE GROUP ===
  int<lower=1, upper=A> sero_age;  // Age group for serology (e.g., 5-14 years)
}

parameters{
  // === DISPERSION (NegBinom k) ===
  real<lower=0> phi_case;          // Case overdispersion
  real<lower=0> phi_death;         // Death overdispersion
  
  // === DETECTION SENSITIVITY (RW by region-year) ===
  array[R] real logit_Se_init;     // Initial logit(Se)
  array[R, Y-1] real logit_Se_innov;  // Year innovations
  real<lower=0> sigma_Se;          // Innovation SD

  // === FORCE OF INFECTION (RW by province-year) ===
  real mu_log_lambda_0;            // Global mean log-FOI
  array[S] real log_lambda_0_s;    // Province baseline
  array[S] real epsilon_lambda_ys_init;  // Initial innovations
  array[Y-1, S] real epsilon_lambda_ys_innov;  // Year innovations

  // === SEASONALITY (Fourier + period effects) ===
  array[H, S] real a_cos_base;     // Baseline cosine
  array[H, S] real a_sin_base;     // Baseline sine
  array[H, P, S] real delta_cos_raw;  // Period deviations (cos)
  array[H, P, S] real delta_sin_raw;  // Period deviations (sin)
  array[S] real<lower=0> tau_cos_base;  // Baseline variance (cos)
  array[S] real<lower=0> tau_sin_base;  // Baseline variance (sin)
  array[P] real<lower=0> sigma_delta_cos;  // Period SD (cos)
  array[P] real<lower=0> sigma_delta_sin;  // Period SD (sin)

  // === AGE-SPECIFIC SUSCEPTIBILITY ===
  array[P, A] real log_rho_pa;     // Log susceptibility (relative to ref_age)
  real<lower=0> sigma_rho;         // Between-period variation

  // === CASE-FATALITY RATE (Model 9 core) ===
  array[P, A] real logit_CFR_pa_raw;  // CFR deviations
  real<lower=0> sigma_CFR;            // CFR variation
  real mu_logit_CFR;                  // Mean logit-CFR prior (~5%)

  // === SEROPREVALENCE ===
  array[S] real<lower=0> H0_s;     // Initial cumulative hazard
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

  // === Age-specific susceptibility: normalize to ref_age = 1 ===
  for(p_ind in 1:P){
    for(a_ind in 1:A){
      rho_pa[p_ind, a_ind] = exp(log_rho_pa[p_ind, a_ind]);
    }
    real ref_rho = rho_pa[p_ind, ref_age];
    for(a_ind in 1:A){
      rho_pa[p_ind, a_ind] = rho_pa[p_ind, a_ind] / ref_rho;  // Constraint
    }
  }

  // === CFR: absolute probability [0,1] by period-age ===
  // inv_logit naturally constrains; no reference normalization
  for(p_ind in 1:P){
    for(a_ind in 1:A){
      CFR_pa[p_ind, a_ind] = inv_logit(mu_logit_CFR + sigma_CFR * logit_CFR_pa_raw[p_ind, a_ind]);
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

  // === Latent cases: I -> C -> D chain ===
  // I = clinical cases (age-adjusted FOI)
  // C = detected cases (I × sensitivity)
  // D = deaths (C × period-age CFR)
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
      lambda_year_ys[y_ind, s_ind] = sum(lambda_ts[t_start:t_end, s_ind]);
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
        int p_here = period_id[y_ind];
        cumulative_hazard += rho_pa[p_here, sero_age] * 
                             lambda_year_ys[y_ind, s_ind] * 
                             exp(-omega * years_since);
      }
      real p_y = 1 - exp(-cumulative_hazard);
      p_sum += p_y;
    }
    p_sero_s[s_ind] = p_sum / K_sero;  // Average across survey years
  }
}

model {
    // === PRIORS ===
    phi_case ~ normal(5, 2);
    phi_death ~ normal(3, 2);

    mu_log_lambda_0 ~ normal(-10, 0.5);
    log_lambda_0_s ~ normal(mu_log_lambda_0, sigma_lambda_0);
    H0_s ~ exponential(1);
    
    epsilon_lambda_ys_init ~ normal(0, sigma_lambda_init);

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
    sigma_Se ~ normal(0.6, 0.01);

    // Age effect priors
    sigma_rho ~ normal(0.3, 0.01);
    for(p_ind in 1:P){
      log_rho_pa[p_ind, ] ~ normal(0, sigma_rho);
    }

    // CFR priors
    mu_logit_CFR ~ normal(logit(0.05), 1);  // Prior ~5% baseline CFR
    sigma_CFR ~ exponential(1);
    for(p_ind in 1:P){
      logit_CFR_pa_raw[p_ind, ] ~ std_normal();
    }

    // === LIKELIHOODS ===
    // Monthly case totals: NegBinom
    for(t_ind in 1:T)
      R_case_ts[t_ind, ] ~ neg_binomial_2(C_ts[t_ind, ] + 1e-10, phi_case);

    // Monthly death totals: NegBinom
    for(t_ind in 1:T)
      R_death_ts[t_ind, ] ~ neg_binomial_2(D_ts[t_ind, ] + 1e-10, phi_death);

    // Annual age-stratified cases: Multinomial (when age data available)
    for(y_ind in 1:Y){
      for(s_ind in 1:S){
        if(has_age_data_ys[y_ind, s_ind]){ // Check if age data available
          vector[A] q_case_yas;
          int N_case_ys = 0;
          for(a_ind in 1:A){
            q_case_yas[a_ind] = C_yas[y_ind, a_ind, s_ind];
            N_case_ys += R_case_yas[y_ind, a_ind, s_ind];
          }
          if(N_case_ys > 0 && sum(q_case_yas) > 0){ // Only model if there are observed cases and expected > 0
            q_case_yas = q_case_yas / sum(q_case_yas);  // Normalize to proportions
            R_case_yas[y_ind, , s_ind] ~ multinomial(q_case_yas);
          }
        }
      }
    }

    // Annual age-stratified deaths: Multinomial (when age-death data available)
    for(y_ind in 1:Y){
      for(s_ind in 1:S){
        if(has_death_data_ys[y_ind, s_ind]){ // Check if age-death data available
          vector[A] q_death_yas;
          int N_death_ys = 0;
          for(a_ind in 1:A){
            q_death_yas[a_ind] = D_yas[y_ind, a_ind, s_ind];
            N_death_ys += R_death_yas[y_ind, a_ind, s_ind];
          }
          if(N_death_ys > 0 && sum(q_death_yas) > 0){ // Only model if there are observed deaths and expected > 0
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
  // === LOG-LIKELIHOOD (for LOO-IC/WAIC model comparison) ===
  vector[T] log_lik_case;           // Per-month log-lik (cases)
  vector[T] log_lik_death;          // Per-month log-lik (deaths) [Model 9]
  vector[Y * S] log_lik_case_age;   // Annual age-stratified (cases)
  vector[Y * S] log_lik_death_age;  // Annual age-stratified (deaths) [Model 9]
  
  // === POSTERIOR PREDICTIVE CHECKS ===
  array[T, S] int R_case_ts_rep;    // Simulated cases (for PPC)
  array[T, S] int R_death_ts_rep;   // Simulated deaths (for PPC) [Model 9]
  
  // === SEASONALITY EFFECTS ===
  array[M, P, S] real seasonal_multiplier;  // exp(Fourier terms)
  
  // === DIAGNOSTICS ===
  real sum_log_lik;                 // Total log-likelihood
  
  // Initialize
  sum_log_lik = 0;
  
  // === Compute log-likelihoods === //
  for(t_ind in 1:T){
    // Case log-likelihood
    log_lik_case[t_ind] = neg_binomial_2_lpmf(R_case_ts[t_ind, 1] | C_ts[t_ind, 1] + 1e-10, phi_case);
    for(s_ind in 2:S){
      log_lik_case[t_ind] += neg_binomial_2_lpmf(R_case_ts[t_ind, s_ind] | C_ts[t_ind, s_ind] + 1e-10, phi_case);
    }
    
    // Death log-likelihood
    log_lik_death[t_ind] = neg_binomial_2_lpmf(R_death_ts[t_ind, 1] | D_ts[t_ind, 1] + 1e-10, phi_death);
    for(s_ind in 2:S){
      log_lik_death[t_ind] += neg_binomial_2_lpmf(R_death_ts[t_ind, s_ind] | D_ts[t_ind, s_ind] + 1e-10, phi_death);
    }
    
    sum_log_lik += log_lik_case[t_ind] + log_lik_death[t_ind];
  }
  
  // Age-stratified case log-likelihood
  {
    int idx = 1;
    for(y_ind in 1:Y){
      for(s_ind in 1:S){
        if(has_age_data_ys[y_ind, s_ind]){
          vector[A] q_case_yas;
          real N_case_ys = 0;
          for(a_ind in 1:A){
            q_case_yas[a_ind] = C_yas[y_ind, a_ind, s_ind];
            N_case_ys += R_case_yas[y_ind, a_ind, s_ind];
          }
          if(N_case_ys > 0 && sum(q_case_yas) > 0){
            q_case_yas = q_case_yas / sum(q_case_yas);
            log_lik_case_age[idx] = multinomial_lpmf(R_case_yas[y_ind, , s_ind] | q_case_yas);
            sum_log_lik += log_lik_case_age[idx];
          } else {
            log_lik_case_age[idx] = 0;
          }
        } else {
          log_lik_case_age[idx] = 0;
        }
        idx += 1;
      }
    }
  }
  
  // Age-stratified death log-likelihood
  {
    int idx = 1;
    for(y_ind in 1:Y){
      for(s_ind in 1:S){
        if(has_death_data_ys[y_ind, s_ind]){
          vector[A] q_death_yas;
          real N_death_ys = 0;
          for(a_ind in 1:A){
            q_death_yas[a_ind] = D_yas[y_ind, a_ind, s_ind];
            N_death_ys += R_death_yas[y_ind, a_ind, s_ind];
          }
          if(N_death_ys > 0 && sum(q_death_yas) > 0){
            q_death_yas = q_death_yas / sum(q_death_yas);
            log_lik_death_age[idx] = multinomial_lpmf(R_death_yas[y_ind, , s_ind] | q_death_yas);
            sum_log_lik += log_lik_death_age[idx];
          } else {
            log_lik_death_age[idx] = 0;
          }
        } else {
          log_lik_death_age[idx] = 0;
        }
        idx += 1;
      }
    }
  }
  
  // === Posterior predictive replicates === //
  for(t_ind in 1:T){
    for(s_ind in 1:S){
      R_case_ts_rep[t_ind, s_ind] = neg_binomial_2_rng(C_ts[t_ind, s_ind] + 1e-10, phi_case);
      R_death_ts_rep[t_ind, s_ind] = neg_binomial_2_rng(D_ts[t_ind, s_ind] + 1e-10, phi_death);
    }
  }
  
  // === Seasonal multipliers === //
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
