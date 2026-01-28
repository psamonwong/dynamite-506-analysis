data {
  // Dimensions
  int<lower=1> T;  // months = Y * M
  int<lower=1> Y;  // years
  int<lower=1> A;  // age groups
  int<lower=1> S;  // provinces
  int<lower=1> R;  // regions
  int<lower=1> P;  // period blocks
  int<lower=1> M;  // months per year
  int<lower=1> H;  // Fourier harmonics

  // Index mappings
  array[T] int<lower=1, upper=Y> year_id;    // t -> y
  array[T] int<lower=1, upper=M> month_id;   // t -> m
  array[Y] int<lower=1, upper=P> period_id;  // y -> p
  array[S] int<lower=1, upper=R> region_id;  // s -> r

  // Population and susceptibility
  array[T,A,S] int<lower=0> N_tas;
  array[T,A,S] real<lower=0, upper=1> S_tas;

  // Observations
  array[T,S] int<lower=0> R_case_ts;      // monthly cases
  array[Y,A,S] int<lower=0> R_case_yas;   // yearly cases by age
}

parameters {
  // Fixed effects
  real lambda_0;               // baseline log FOI
  real alpha_0_tau;            // baseline log testing intensity
  real alpha_0_se;             // baseline logit sensitivity
  real gamma;                  // effect of testing on sensitivity
  real<lower=0, upper=1> pi_fp; // false positive rate
  real<lower=0, upper=1> nu;    // clinical fraction

  // Seasonality baselines
  array[H] real A0;  // cosine
  array[H] real B0;  // sine

  // Age effects
  real beta_0;            // age intercept
  vector[A-1] beta_a_raw; // age deviations (centered)

  // SDs (hierarchical scales)
  real<lower=0> sigma_tau_pr;
  real<lower=0> sigma_tau_ys;
  real<lower=0> sigma_lambda_pr;
  real<lower=0> sigma_lambda_ys;
  real<lower=0> sigma_seasonality_pr;
  real<lower=0> sigma_age_pr;
  real<lower=0> sigma_age_ys;
  real<lower=0> sigma_age;

  // Overdispersion
  real<lower=0> phi_c;

  // Non-centered random effects (standard normal)
  array[P,R] real z_tau_pr;
  array[Y,S] real z_tau_ys;
  array[P,R] real z_lambda_pr;
  array[Y,S] real z_lambda_ys;
  array[H,P,R] real z_A_pr;
  array[H,P,R] real z_B_pr;
  array[P,R] real z_age_pr;
  array[Y,S] real z_age_ys;
}

transformed parameters {
  // Random effects (scaled)
  array[P,R] real epsilon_tau_pr;
  array[Y,S] real epsilon_tau_ys;
  array[P,R] real epsilon_lambda_pr;
  array[Y,S] real epsilon_lambda_ys;
  array[H,P,R] real epsilon_A_pr;
  array[H,P,R] real epsilon_B_pr;
  array[P,R] real epsilon_age_pr;
  array[Y,S] real epsilon_age_ys;

  // Small derived arrays we may want to save
  array[Y,S] real<lower=0> tau_ys;
  array[Y,S] real log_lambda_base;

  // Expected values
  array[T,S] real<lower=0> mu_case_ts;
  array[Y,A,S] real C_yas;                      // yearly totals by age (for composition)
  array[Y,S] vector[A] q_case_ys;               // age composition proportions

  // Scale random effects
  for (p in 1:P) for (r in 1:R) {
    epsilon_tau_pr[p,r] = sigma_tau_pr * z_tau_pr[p,r];
    epsilon_lambda_pr[p,r] = sigma_lambda_pr * z_lambda_pr[p,r];
    epsilon_age_pr[p,r] = sigma_age_pr * z_age_pr[p,r];
    for (h in 1:H) {
      epsilon_A_pr[h,p,r] = sigma_seasonality_pr * z_A_pr[h,p,r];
      epsilon_B_pr[h,p,r] = sigma_seasonality_pr * z_B_pr[h,p,r];
    }
  }
  for (y in 1:Y) for (s in 1:S) {
    epsilon_tau_ys[y,s] = sigma_tau_ys * z_tau_ys[y,s];
    epsilon_lambda_ys[y,s] = sigma_lambda_ys * z_lambda_ys[y,s];
    epsilon_age_ys[y,s] = sigma_age_ys * z_age_ys[y,s];
    tau_ys[y,s] = exp(alpha_0_tau +
                      epsilon_tau_pr[period_id[y], region_id[s]] +
                      epsilon_tau_ys[y,s]);
    log_lambda_base[y,s] = lambda_0 +
                           epsilon_lambda_pr[period_id[y], region_id[s]] +
                           epsilon_lambda_ys[y,s];
  }
  

  // Initialize accumulators
  for (y in 1:Y) for (a in 1:A) for (s in 1:S) C_yas[y,a,s] = 0;

  // Compute monthly totals and accumulate yearly-by-age
  for (t in 1:T) {
    int y = year_id[t];
    int m = month_id[t];
    for (s in 1:S) {
      // Sensitivity depends on testing
      real Se = inv_logit(alpha_0_se + gamma * log(tau_ys[y,s]));
      real mu_sum = 0;

      // Precompute seasonality for (y,m,s)
      real theta = 2 * pi() * m / M;
      real log_kappa_m = 0;
      for (h in 1:H) {
        real A_loc = A0[h] + epsilon_A_pr[h, period_id[y], region_id[s]];
        real B_loc = B0[h] + epsilon_B_pr[h, period_id[y], region_id[s]];
        log_kappa_m += A_loc * cos(h * theta) + B_loc * sin(h * theta);
      }

      for (a in 1:A) {
        // Age multiplier
        real beta_a = (a == 1) ? 0 : beta_a_raw[a-1];
        real log_kappa_a = beta_0 + beta_a +
                           epsilon_age_pr[period_id[y], region_id[s]] +
                           epsilon_age_ys[y,s];

        // False positives
        real FP = N_tas[t,a,s] * tau_ys[y,s] * pi_fp;

        // FOI and clinical infections
        real lambda = exp(log_lambda_base[y,s] + log_kappa_m + log_kappa_a);
        real I = nu * S_tas[t,a,s] * N_tas[t,a,s] * (1 - exp(-lambda));

        // Reported cases
        real C = I * Se + FP;
        mu_sum += C;
        C_yas[y,a,s] += C;
      }
      mu_case_ts[t,s] = mu_sum;
    }
  }

  // Age composition proportions
  for (y in 1:Y) for (s in 1:S) {
    vector[A] C_age;
    for (a in 1:A) C_age[a] = C_yas[y,a,s];
    if (sum(C_age) > 0) q_case_ys[y,s] = C_age / sum(C_age);
    else q_case_ys[y,s] = rep_vector(1.0 / A, A);
  }
}

model {
  // Priors - fixed effects
  lambda_0 ~ normal(0, 5);
  alpha_0_tau ~ normal(0, 5);
  alpha_0_se ~ normal(0, 5);
  gamma ~ normal(0, 2);
  beta_0 ~ normal(0, 5);
  pi_fp ~ beta(2, 98); // ~2%
  nu ~ beta(2, 2);     // ~50%
  for (h in 1:H) { A0[h] ~ normal(0, 2); B0[h] ~ normal(0, 2); }
  beta_a_raw ~ normal(0, sigma_age);
  sigma_age ~ normal(0, 5);

  // Priors - scales
  sigma_tau_pr ~ normal(0, 5);
  sigma_tau_ys ~ normal(0, 5);
  sigma_lambda_pr ~ normal(0, 5);
  sigma_lambda_ys ~ normal(0, 5);
  sigma_seasonality_pr ~ normal(0, 5);
  sigma_age_pr ~ normal(0, 5);
  sigma_age_ys ~ normal(0, 5);
  phi_c ~ normal(0, 5);

  // Non-centered random effects
  for(p in 1:P) for(r in 1:R) {
    z_tau_pr[p,r] ~ std_normal();
    z_lambda_pr[p,r] ~ std_normal();
    z_age_pr[p,r] ~ std_normal();
    for(h in 1:H) {
      z_A_pr[h,p,r] ~ std_normal();
      z_B_pr[h,p,r] ~ std_normal();
    }
  }
  
  for(y in 1:Y) for(s in 1:S) {
    z_tau_ys[y,s] ~ std_normal();
    z_lambda_ys[y,s] ~ std_normal();
    z_age_ys[y,s] ~ std_normal();
  }


//   to_vector(z_tau_pr) ~ std_normal();
//   to_vector(z_tau_ys) ~ std_normal();
//   to_vector(z_lambda_pr) ~ std_normal();
//   to_vector(z_lambda_ys) ~ std_normal();
//   to_vector(z_A_pr) ~ std_normal();
//   to_vector(z_A_ys) ~ std_normal();
//   to_vector(z_B_pr) ~ std_normal();
//   to_vector(z_B_ys) ~ std_normal();
//   to_vector(z_age_pr) ~ std_normal();
//   to_vector(z_age_ys) ~ std_normal();

  // Likelihoods
  for (t in 1:T) for (s in 1:S)
    R_case_ts[t,s] ~ neg_binomial_2(mu_case_ts[t,s], phi_c);

  for (y in 1:Y) for (s in 1:S) {
    array[A] int R_case_age;
    int N_case = 0;
    for (a in 1:A) { R_case_age[a] = R_case_yas[y,a,s]; N_case += R_case_age[a]; }
    if (N_case > 0) R_case_age ~ multinomial(q_case_ys[y,s]);
  }
}
