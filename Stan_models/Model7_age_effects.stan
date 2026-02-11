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
  
  // Susceptible proportions
  int S_t;                            // Susceptible proportion (assumed constant)

  // Seasonality configurations
  int<lower=1> H;                 // number of Fourier harmonics
  matrix[T, H] cos_t_h;           // Precomputed cosine basis [T x H] for each time point
  matrix[T, H] sin_t_h;           // Precomputed sine basis [T x H] for each time point

  // Reference age group for susceptibility
  int<lower=1, upper=A> ref_age;  // Reference age group index (1 to A)

  // Serology data
  array[S] real<lower=0, upper=1> sero_lower_s;
  array[S] real<lower=0, upper=1> sero_upper_s;
  array[K_sero] int sero_year_id;

  // Immune waning rate
  real<lower=0> omega;

  // Clinical fraction
  real<lower=0, upper=1> nu;
}

parameters{
  // === Dispersion parameters === //
  real<lower=0> phi_case;               // Dispersion for monthly totals
  
  // === Test sensitivity and specificity parameters === //
  array[R] real logit_Se_init;
  array[R, Y-1] real logit_Se_innov;

  // === Force of infection parameters === //
  real mu_log_lambda_0;
  array[S] real log_lambda_0_s;
  real<lower=0> sigma_lambda_0;
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

  array[S] real<lower=0> H0_s;

  // === Hyperparameters === //
  real<lower=0> sigma_lambda_y;
  real<lower=0> sigma_lambda_init;
  real<lower=0> sigma_Se;
}

transformed parameters {
  // === Force of infection components === //
  array[T, S] real lambda_ts;            // Force of infection at time t, province s (baseline * seasonal)
  array[Y, S] real lambda_base_ys;       // Baseline FOI at year y, province s
  array[Y, S] real lambda_year_ys;       // Annual aggregated FOI (for seroprevalence)
  
  // === Sensitivity and specificity === //
  array[R, Y] real<lower=0, upper=1> Se_y;  // Case detection sensitivity by region-year
  array[Y] real<lower=0, upper=1> Sp_y;     // Specificity (set to 1)
  
  // === Expected latent cases and reported cases === //
  array[T, A, S] real<lower=0> I_tas;    // Monthly latent true clinical cases (before detection)
  array[T, A, S] real<lower=0> C_tas;    // Monthly reported cases (after detection sensitivity)
  array[T, S] real<lower=0> I_ts;        // Monthly latent total cases (sum over age)
  array[T, S] real<lower=0> C_ts;        // Monthly reported total cases (sum over age)
  
  // === Annual age-stratified cases === //
  array[Y, A, S] real<lower=0> I_yas;    // Annual latent cases by age
  array[Y, A, S] real<lower=0> C_yas;    // Annual reported cases by age
  
  // === Seroprevalence === //
  array[S] real<lower=0, upper=1> p_sero_s;  // Province-level model-predicted seroprevalence
  
  // === Age and period effects === //
  array[H, P, S] real a_cos;
  array[H, P, S] real a_sin;
  array[P, A] real<lower=0> rho_pa;      // Age-specific relative susceptibility
  
  // === Random walk processes === //
  array[Y, S] real epsilon_lambda_ys;

  // Initialize specificity
  for(i in 1:Y) Sp_y[i] = 1;

  // === Build period-specific seasonal coefficients === //
  for(h in 1:H){
    for(p in 1:P){
      for(s in 1:S){
        a_cos[h, p, s] = tau_cos_base[s] * a_cos_base[h, s] + 
                         sigma_delta_cos[p] * delta_cos_raw[h, p, s];
        a_sin[h, p, s] = tau_sin_base[s] * a_sin_base[h, s] + 
                         sigma_delta_sin[p] * delta_sin_raw[h, p, s];
      }
    }
  }
  
  // === Build sensitivity random walk === //
  {
    array[R, Y] real logit_Se_y;
    for(r_ind in 1:R){
      logit_Se_y[r_ind, 1] = logit_Se_init[r_ind];
      for (y in 2:Y) {
        logit_Se_y[r_ind, y] = logit_Se_y[r_ind, y-1] + logit_Se_innov[r_ind, y-1];
      }
      for (y in 1:Y) {
        Se_y[r_ind, y] = inv_logit(logit_Se_y[r_ind, y]);
      }
    }
  }
  
  // === Build province-specific random walks === //
  for (s_ind in 1:S) {
    epsilon_lambda_ys[1, s_ind] = epsilon_lambda_ys_init[s_ind];
    for (y in 2:Y) {
      epsilon_lambda_ys[y, s_ind] = epsilon_lambda_ys[y-1, s_ind] + epsilon_lambda_ys_innov[y-1, s_ind];
    }
  }

  // === Baseline force of infection === //
  for(y_ind in 1:Y)
    for(s_ind in 1:S)
      lambda_base_ys[y_ind, s_ind] = exp(log_lambda_0_s[s_ind] + epsilon_lambda_ys[y_ind, s_ind]);

  // === Age-specific susceptibility (constraint: reference age = 1) === //
  for(p in 1:P){
    for(a in 1:A){
      rho_pa[p, a] = exp(log_rho_pa[p, a]);
    }
    // Normalize: set reference age group to rho = 1
    real ref_rho = rho_pa[p, ref_age];
    for(a in 1:A){
      rho_pa[p, a] = rho_pa[p, a] / ref_rho;
    }
  }

  // === Force of infection by time and site (using precomputed Fourier basis) === //
  for(t_ind in 1:T){
    int y = year_id[t_ind];
    int p = period_id[y];
    for(s_ind in 1:S){
      real f_seasonal = 0;
      for(h in 1:H){
        f_seasonal += a_cos[h, p, s_ind] * cos_t_h[t_ind, h] +
                      a_sin[h, p, s_ind] * sin_t_h[t_ind, h];
      }
      lambda_ts[t_ind, s_ind] = lambda_base_ys[y, s_ind] * exp(f_seasonal);
    }
  }

  // === Latent true clinical cases (age-stratified, monthly) === //
  // I_{t,a,s} = nu * N_{t,a,s} * (1 - exp(-lambda_{t,a,s}))
  for(t_ind in 1:T){
    int y = year_id[t_ind];
    int p = period_id[y];
    for(s_ind in 1:S){
      real I_sum = 0;
      for(a_ind in 1:A){
        real lambda_ta = lambda_ts[t_ind, s_ind] * rho_pa[p, a_ind];
        I_tas[t_ind, a_ind, s_ind] = nu * N_tas[t_ind, a_ind, s_ind] * (1 - exp(-lambda_ta));
        I_sum += I_tas[t_ind, a_ind, s_ind];
      }
      I_ts[t_ind, s_ind] = I_sum;
    }
  }

  // === Reported cases after detection sensitivity === //
  // C_{t,a,s} = Se_{y(t),r(s)} * I_{t,a,s}
  for(t_ind in 1:T){
    int y = year_id[t_ind];
    for(s_ind in 1:S){
      int r = region_id[s_ind];
      real C_sum = 0;
      for(a_ind in 1:A){
        C_tas[t_ind, a_ind, s_ind] = Se_y[r, y] * I_tas[t_ind, a_ind, s_ind];
        C_sum += C_tas[t_ind, a_ind, s_ind];
      }
      C_ts[t_ind, s_ind] = C_sum;
    }
  }

  // === Annual aggregations === //
  for(y_ind in 1:Y){
    for(a_ind in 1:A){
      for(s_ind in 1:S){
        real I_sum = 0;
        real C_sum = 0;
        for(m_ind in 1:M){
          int t = (y_ind - 1) * M + m_ind;
          I_sum += I_tas[t, a_ind, s_ind];
          C_sum += C_tas[t, a_ind, s_ind];
        }
        I_yas[y_ind, a_ind, s_ind] = I_sum;
        C_yas[y_ind, a_ind, s_ind] = C_sum;
      }
    }
  }

  // === Annual aggregated force of infection (for seroprevalence) === //
  for(y_ind in 1:Y){
    for(s_ind in 1:S){
      real lambda_sum = 0;
      for(m_ind in 1:M){
        int t = (y_ind - 1) * M + m_ind;
        lambda_sum += lambda_ts[t, s_ind];
      }
      lambda_year_ys[y_ind, s_ind] = lambda_sum / 12;  // Monthly average
    }
  }

  // === Seroprevalence calculation with waning immunity === //
  for(s_ind in 1:S){
    real p_sum = 0;
    for (k in 1:K_sero) {
      int y0 = sero_year_id[k];
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

    mu_log_lambda_0 ~ normal(-10, 2);
    sigma_lambda_0 ~ normal(0.1, 0.05);
    for(s in 1:S){
      log_lambda_0_s[s] ~ normal(mu_log_lambda_0, sigma_lambda_0);
      H0_s[s] ~ exponential(1);
    }
    epsilon_lambda_ys_init ~ normal(0, sigma_lambda_init);
    sigma_lambda_init ~ exponential(10);

    for (y in 1:(Y-1)) {
      epsilon_lambda_ys_innov[y] ~ normal(0, sigma_lambda_y);
    }

    // Seasonal variation priors
    for(h in 1:H){
      for(s in 1:S){
        a_cos_base[h, s] ~ normal(0, 1);
        a_sin_base[h, s] ~ normal(0, 1);
      }
    }
    
    for(s in 1:S){
      tau_cos_base[s] ~ normal(0, 0.5);
      tau_sin_base[s] ~ normal(0, 0.5);
    }
    
    sigma_delta_cos[1] ~ normal(0, 0.05);
    sigma_delta_sin[1] ~ normal(0, 0.05);
    for(p in 2:P){
      sigma_delta_cos[p] ~ normal(sigma_delta_cos[p-1], 0.05);
      sigma_delta_sin[p] ~ normal(sigma_delta_sin[p-1], 0.05);
    }
    
    for(h in 1:H){
      for(p in 1:P){
        for(s in 1:S){
          delta_cos_raw[h, p, s] ~ std_normal();
          delta_sin_raw[h, p, s] ~ std_normal();
        }
      }
    }

    // Sensitivity priors
    for(r in 1:R)
      for(y in 1:(Y-1))
        logit_Se_innov[r, y] ~ normal(0, sigma_Se);

    logit_Se_init ~ normal(logit(0.4), 0.3);
    sigma_Se ~ exponential(1);

    // Age effect priors
    sigma_rho ~ exponential(1);
    for(p in 1:P){
      for(a in 1:A){
        log_rho_pa[p, a] ~ normal(0, sigma_rho);
      }
    }

    // Hyperparameter priors
    sigma_lambda_y ~ normal(0.05, 0.01);

    // === Likelihoods === //
    // Monthly case totals: negative binomial
    for(t_ind in 1:T)
      for(s_ind in 1:S)
        R_case_ts[t_ind, s_ind] ~ neg_binomial_2(C_ts[t_ind, s_ind] + 1e-9, phi_case);

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

    // Seroprevalence likelihood: interval-censored normal
    for(s_ind in 1:S){
      real sigma_sero = 0.01;
      target += log_diff_exp(
        std_normal_lcdf((sero_upper_s[s_ind] - p_sero_s[s_ind]) / sigma_sero),
        std_normal_lcdf((sero_lower_s[s_ind] - p_sero_s[s_ind]) / sigma_sero)
      );
    }
}
