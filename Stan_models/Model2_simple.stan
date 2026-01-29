data {
  // ===== Dimensions =====
  int<lower=1> T;  // Number of months (= Y * M)
  int<lower=1> Y;  // Number of years
  int<lower=1> A;  // Number of age groups
  int<lower=1> S;  // Number of provinces
  int<lower=1> R;  // Number of regions
  int<lower=1> P;  // Number of period blocks
  int<lower=1> M;  // Months per year (=12)
  int<lower=1> H;  // Number of Fourier harmonics for seasonality
  
  // ===== Precomputed seasonality basis (H x M) =====
  matrix[H, M] cos_hm;  // cosine basis for each harmonic and month
  matrix[H, M] sin_hm;  // sine basis for each harmonic and month
  
  // ===== Index mappings =====
  array[T] int<lower=1, upper=Y> year_id;    // Map month t to year y
  array[T] int<lower=1, upper=M> month_id;   // Map month t to month-of-year m
  array[Y] int<lower=1, upper=P> period_id;  // Map year y to period block p
  array[S] int<lower=1, upper=R> region_id;  // Map province s to region r
  
  // ===== Population and susceptibility =====
  array[T,A,S] int<lower=0> N_tas;                   // Population size
  array[T,A,S] real<lower=0, upper=1> S_tas;        // Proportion susceptible
  
  // ===== Observed data =====
  array[T,S] int<lower=0> R_case_ts;      // Monthly reported cases
  array[Y,A,S] int<lower=0> R_case_yas;   // Yearly reported cases by age
}

parameters {
  // ===== Fixed effects =====
  real lambda_0;                         // Baseline log force of infection
  real alpha_0_tau;                      // Baseline log testing intensity
  real<lower=0, upper=1> alpha_0_se;     // Baseline sensitivity (probability scale)
  real<lower=0, upper=1> pi_fp;         // False positive rate
  
  // ===== Fourier seasonality baselines (period–region) =====
  array[H, P, R] real A_pr;  // Cosine coefficients by harmonic/period/region
  array[H, P, R] real B_pr;  // Sine coefficients by harmonic/period/region
  
  // ===== Age effect baselines =====
  vector[A-1] beta_a_raw;     // Age-specific deviations (centered, first group = 0)
  
  // ===== Hierarchical SDs =====
  real<lower=0> sigma_tau_pr;            // Testing: period-region
  real<lower=0> sigma_tau_ys;            // Testing: year-province
  real<lower=0> sigma_lambda_pr;         // FOI: period-region
  real<lower=0> sigma_lambda_ys;         // FOI: year-province
  real<lower=0> sigma_A_ys;              // Seasonality cos: year-province noise
  real<lower=0> sigma_B_ys;              // Seasonality sin: year-province noise
  real<lower=0> sigma_age_pr;            // Age: period-region
  real<lower=0> sigma_age;               // Age: between age groups
  
  // ===== Overdispersion =====
  real<lower=0> phi_c;  // Cases
  
  // ===== Random effects (non-centered: z-scores) =====
  array[P,R] real z_tau_pr;        // Testing: period-region
  array[Y,S] real z_tau_ys;        // Testing: year-province
  array[P,R] real z_lambda_pr;     // FOI: period-region
  array[Y,S] real z_lambda_ys;     // FOI: year-province
  array[H,Y,S] real z_A_ys;        // Seasonality cos: year-province noise
  array[H,Y,S] real z_B_ys;        // Seasonality sin: year-province noise
  array[P,R] real z_age_pr;        // Age: period-region
}


transformed parameters {
  // ===== Scaled random effects (from non-centered z-scores) =====
  array[P,R] real epsilon_tau_pr;
  array[Y,S] real epsilon_tau_ys;
  array[P,R] real epsilon_lambda_pr;
  array[Y,S] real epsilon_lambda_ys;
  array[H,Y,S] real epsilon_A_ys;
  array[H,Y,S] real epsilon_B_ys;
  array[P,R] real epsilon_age_pr;
  
  
  // ===== Derived parameters =====
  array[Y,S] real<lower=0> tau_ys;              // Testing intensity
  array[Y,S] real log_lambda_base;              // Baseline log FOI
  array[Y,M,S] real kappa_month_yms;            // Monthly seasonality multipliers
  array[Y,A,S] real kappa_age_yas;              // Age-specific multipliers
  vector[A] beta_a;                             // Full age effects (incl. reference)
  
  // ===== Expected values =====
  array[T,S] real<lower=0> mu_case_ts;          // Monthly mean cases
  array[Y,A,S] real C_yas;                      // Yearly cases by age
  array[Y,S] vector[A] q_case_ys;               // Age composition proportions
  
  // ===== Scale random effects (non-centered parameterization) =====
  for (p_ind in 1:P) {
    for (r_ind in 1:R) {
      epsilon_tau_pr[p_ind, r_ind] = z_tau_pr[p_ind, r_ind] * sigma_tau_pr;
      epsilon_lambda_pr[p_ind, r_ind] = z_lambda_pr[p_ind, r_ind] * sigma_lambda_pr;
      epsilon_age_pr[p_ind, r_ind] = z_age_pr[p_ind, r_ind] * sigma_age_pr;
    }
  }
  
  for (h_ind in 1:H) {
    for (y_ind in 1:Y) {
      for (s_ind in 1:S) {
        epsilon_A_ys[h_ind, y_ind, s_ind] = z_A_ys[h_ind, y_ind, s_ind] * sigma_A_ys;
        epsilon_B_ys[h_ind, y_ind, s_ind] = z_B_ys[h_ind, y_ind, s_ind] * sigma_B_ys;
      }
    }
  }
  
  for (y_ind in 1:Y) {
    for (s_ind in 1:S) {
      epsilon_tau_ys[y_ind, s_ind] = z_tau_ys[y_ind, s_ind] * sigma_tau_ys;
      epsilon_lambda_ys[y_ind, s_ind] = z_lambda_ys[y_ind, s_ind] * sigma_lambda_ys;
    }
  }
  
  
  
  // ===== Compute baseline FOI =====
  for (y_ind in 1:Y) {
    for (s_ind in 1:S) {
      log_lambda_base[y_ind, s_ind] = lambda_0 + 
        epsilon_lambda_pr[period_id[y_ind], region_id[s_ind]] + 
        epsilon_lambda_ys[y_ind, s_ind];
    }
  }
  
  // ===== Compute seasonality multipliers using precomputed basis =====
  for (y_ind in 1:Y) {
    for (s_ind in 1:S) {
      vector[M] log_kappa_yms;
      // Build local seasonality coefficient vectors (length H)
      vector[H] Avec;
      vector[H] Bvec;
      for (h_ind in 1:H) {
        Avec[h_ind] = A_pr[h_ind, period_id[y_ind], region_id[s_ind]] +
          epsilon_A_ys[h_ind, y_ind, s_ind];
        Bvec[h_ind] = B_pr[h_ind, period_id[y_ind], region_id[s_ind]] +
          epsilon_B_ys[h_ind, y_ind, s_ind];
      }
      // Project onto month-specific basis columns
      for (m_ind in 1:M) {
        log_kappa_yms[m_ind] = dot_product(Avec, cos_hm[, m_ind]) +
                               dot_product(Bvec, sin_hm[, m_ind]);
      }
      // Normalize to mean 1 (identifiability constraint)
      real sum_exp_m = 0;
      for (m_ind in 1:M) sum_exp_m += exp(log_kappa_yms[m_ind]);
      for (m_ind in 1:M) {
        kappa_month_yms[y_ind, m_ind, s_ind] = M * exp(log_kappa_yms[m_ind]) / sum_exp_m;
      }
    }
  }
  
  // ===== Compute age-specific multipliers =====
  beta_a[1] = 0;  // Reference age group
  beta_a[2:A] = beta_a_raw;
  
  for (y_ind in 1:Y) {
    for (s_ind in 1:S) {
      vector[A] log_kappa_age;
      for (a_ind in 1:A) {
        log_kappa_age[a_ind] = beta_a[a_ind] + 
          epsilon_age_pr[period_id[y_ind], region_id[s_ind]];
        kappa_age_yas[y_ind, a_ind, s_ind] = exp(log_kappa_age[a_ind]);
      }
      // Normalize to mean 1 across ages (identifiability constraint)
      real sum_kappa_a = 0;
      for (a_ind in 1:A) sum_kappa_a += kappa_age_yas[y_ind, a_ind, s_ind];
      for (a_ind in 1:A) {
        kappa_age_yas[y_ind, a_ind, s_ind] = A * kappa_age_yas[y_ind, a_ind, s_ind] / sum_kappa_a;
      }
    }
  }
  
  // ===== Compute testing intensity =====
  for (y_ind in 1:Y) {
    for (s_ind in 1:S) {
      tau_ys[y_ind, s_ind] = exp(alpha_0_tau + 
        epsilon_tau_pr[period_id[y_ind], region_id[s_ind]] + 
        epsilon_tau_ys[y_ind, s_ind]);
    }
  }
  
  // ===== Initialize yearly aggregates =====
  for (y_ind in 1:Y) {
    for (a_ind in 1:A) {
      for (s_ind in 1:S) {
        C_yas[y_ind, a_ind, s_ind] = 0;
      }
    }
  }
  
  // ===== Compute monthly totals and accumulate yearly by age =====
  for (t_ind in 1:T) {
    int y_ind = year_id[t_ind];
    int m_ind = month_id[t_ind];
    
    for (s_ind in 1:S) {
      real mu_I_sum = 0;   // sum of clinical infections across ages
      real mu_FP_sum = 0;  // sum of false positives across ages
      
      for (a_ind in 1:A) {
        int N_local = N_tas[t_ind, a_ind, s_ind];
        real lambda_local = exp(log_lambda_base[y_ind, s_ind]) *
          kappa_month_yms[y_ind, m_ind, s_ind] *
          kappa_age_yas[y_ind, a_ind, s_ind];
        real I_local = S_tas[t_ind, a_ind, s_ind] * N_local * (1 - exp(-lambda_local));
        real FP_local = pi_fp * lambda_local * N_local;  // FP scales with FOI
        
        mu_I_sum += I_local;   // accumulate clinical infections
        mu_FP_sum += FP_local; // accumulate false positives
        
        // Update yearly-by-age
        C_yas[y_ind, a_ind, s_ind] += I_local * alpha_0_se + FP_local;
      }
      // Monthly total: infections * Se + FP
      mu_case_ts[t_ind, s_ind] = mu_I_sum * alpha_0_se + mu_FP_sum;
    }
  }
  
  // ===== Compute age composition =====
  for (y_ind in 1:Y) {
    for (s_ind in 1:S) {
      vector[A] C_age;
      for (a_ind in 1:A) {
        C_age[a_ind] = C_yas[y_ind, a_ind, s_ind];
      }
      if (sum(C_age) > 0) {
        q_case_ys[y_ind, s_ind] = C_age / sum(C_age);
      } else {
        q_case_ys[y_ind, s_ind] = rep_vector(1.0 / A, A);
      }
    }
  }
}


model {
  // ========================================
  // PRIORS
  // ========================================
  
  // ----- Fixed effects -----
  lambda_0 ~ normal(-3, 2);      // Baseline log FOI centered at reasonable value
  alpha_0_tau ~ normal(-2, 2);   // Baseline log testing intensity
  
  // ----- Bounded parameters -----
  alpha_0_se ~ beta(50, 10);     // Prior mean ~0.83 (reasonable sensitivity)
  pi_fp ~ beta(2, 98);           // Prior mean ~2%
  
  // ----- Seasonality baselines -----
  for (h in 1:H) {
    for (p_ind in 1:P) {
      for (r_ind in 1:R) {
        A_pr[h, p_ind, r_ind] ~ normal(0, 0.5);
        B_pr[h, p_ind, r_ind] ~ normal(0, 0.5);
      }
    }
  }
  
  // ----- Age effects -----
  beta_a_raw ~ normal(0, sigma_age);
  sigma_age ~ normal(0, 0.1);
  
  // ----- Hierarchical SDs -----
  sigma_tau_pr ~ exponential(1);
  sigma_tau_ys ~ exponential(10);    // mean = 0.1
  sigma_lambda_pr ~ exponential(1);
  sigma_lambda_ys ~ exponential(10); // mean = 0.1
  sigma_A_ys ~ exponential(5);       // mean = 0.2
  sigma_B_ys ~ exponential(5);       // mean = 0.2
  sigma_age_pr ~ exponential(1);
  
  // ----- Overdispersion -----
  phi_c ~ normal(5, 0.25);  // Better prior for positive scale parameter
  
  // ----- Non-centered random effects (z-scores) -----
  for (p_ind in 1:P) {
    for (r_ind in 1:R) {
      z_tau_pr[p_ind, r_ind] ~ std_normal();
      z_lambda_pr[p_ind, r_ind] ~ std_normal();
      z_age_pr[p_ind, r_ind] ~ std_normal();
    }
  }
  
  for (y_ind in 1:Y) {
    for (s_ind in 1:S) {
      z_tau_ys[y_ind, s_ind] ~ std_normal();
      z_lambda_ys[y_ind, s_ind] ~ std_normal();
      for (h_ind in 1:H) {
        z_A_ys[h_ind, y_ind, s_ind] ~ std_normal();
        z_B_ys[h_ind, y_ind, s_ind] ~ std_normal();
      }
    }
  }
  
  // ========================================
  // LIKELIHOODS
  // ========================================
  
  // ----- Monthly case counts (Negative Binomial - vectorized) -----
  {
    vector[T * S] mu_flat;
    array[T * S] int R_flat;
    int idx = 1;
    for (t_ind in 1:T) {
      for (s_ind in 1:S) {
        mu_flat[idx] = mu_case_ts[t_ind, s_ind];
        R_flat[idx] = R_case_ts[t_ind, s_ind];
        idx += 1;
      }
    }
    R_flat ~ neg_binomial_2(mu_flat, phi_c);
  }
  
  // ----- Yearly age composition (Multinomial) -----
  for (y_ind in 1:Y) {
    for (s_ind in 1:S) {
      array[A] int R_case_age;
      int N_case_ys = 0;
      
      for (a_ind in 1:A) {
        R_case_age[a_ind] = R_case_yas[y_ind, a_ind, s_ind];
        N_case_ys += R_case_age[a_ind];
      }
      
      if (N_case_ys > 0) {
        R_case_age ~ multinomial(q_case_ys[y_ind, s_ind]);
      }
    }
  }
}
