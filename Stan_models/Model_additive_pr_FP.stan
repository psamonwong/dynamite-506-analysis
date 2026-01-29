data {
  // Dimensions
  int<lower=1> T;                 // Number of months = Y * M
  int<lower=1> A;                 // number of age strata
  int<lower=1> Y;                 // number of years
  int<lower=1> M;                 // number of months per year
  int<lower=1> P;                 // number of period blocks
  int<lower=1> S;                 // number of provinces
  int<lower=1> R;                 // number of regions

  // Mappings
  array[T] int<lower=1, upper=Y> year_id;         // y(t): year corresponding to month t
  array[T] int<lower=1, upper=M> month_id;        // m(t): month corresponding to month t
  array[Y] int<lower=1, upper=P> period_id;       // p(y): period block
  array[S] int<lower=1, upper=R> region_id;       // r(s): region corresponding to province s

  // Population size and susceptibility
  array[T, A, S] int<lower=0> N_tas;          // N_t,a,s population size at time t, age a, province s
  array[T, A, S] real<lower=0, upper=1> S_tas;   // S_t,a,s susceptible proportion at time t, age a, province s

  // Reported monthly totals
  array[T, S] int<lower=0> R_case_ts;          // reported cases at time t, province s
  array[T, S] int<lower=0> R_death_ts;         // reported deaths at time t, province s

  // Yearly by-age counts
  array[Y, A, S] int<lower=0> R_case_yas;
  array[Y, A, S] int<lower=0> R_death_yas;

  // Seasonality configurations
  int<lower=1> H;                 // number of Fourier harmonics
}

transformed data {
  real sigma_gamma_CFR = 0.01;  // fixed small value — data don't inform this
}

parameters {
  // Overdispersion parameters for monthly totals (NB2)
  real<lower=1e-09> phi_c;
  real<lower=1e-09> phi_d;

  // Clinical fraction (probability an infection is clinical)
  real<lower=0, upper=1> nu;

  // FOI baseline on log scale: additive structure
  real alpha_lambda;                   // global baseline
  array[S] real delta_lambda_s;        // province deviation
  array[Y] real delta_lambda_y;        // year deviation
  real<lower=0> sigma_lambda_noise;    // SD for province-year noise
  array[Y, S] real epsilon_lambda_ys;  // province-year noise

  // Seasonality (Fourier): region–period baseline + year–province deviations (raw)
  array[H, P, R] real a_cos;        // cosine coefficients by harmonic/period/region
  array[H, P, R] real a_sin;        // sine coefficients by harmonic/period/region
  array[H, Y, S] real raw_delta_cos;    // year–province cosine deviations (raw, centered in TP)
  array[H, Y, S] real raw_delta_sin;    // year–province sine deviations (raw, centered in TP)

  // Age transmission modifier (Lee–Carter)
  array[A] real beta_kappa_a;                  // baseline age-specific RR
  array[R] real gamma_kappa_pr_base;           // period-1 anchor per region (>0 by constraint below)
  array[P-1, R] real gamma_kappa_pr_innov;     // RW innovations for periods 2..P
  array[Y, S] real raw_delta_kappa_ys;         // province–year deviation (raw, centered in TP)
  ordered[A] raw_h_kappa_a;                    // monotone-increasing raw age-loading

  // Case detection sensitivity (Se) on logit scale (non-centered RW)
  array[R] real alpha_Se_yr_base;             // year-1 anchor per region
  array[Y-1, R] real alpha_Se_yr_innov;       // RW innovations for years 2..Y
  array[S] real raw_delta_Se_s;               // province deviations (raw, centered in TP)
  real<lower=0> sigma_alpha_Se;               // RW innovation SD for alpha_Se_yr
  real<lower=0> sigma_delta_Se;               // SD for delta_Se_s

  // False positives (background reporting) on log scale: alpha_pr + delta_s + delta_y + noise
  array[R, P] real alpha_FP_pr;               // region-period baselines
  array[S] real raw_delta_FP_s;               // province deviations (raw, centered in TP)
  array[Y] real delta_FP_y;                   // global year deviations (centered via prior)
  real<lower=0> sigma_delta_FP;               // SD for delta_FP_s
  real<lower=0> sigma_FP_noise;               // SD for province-year noise
  array[Y, S] real eta_FP_ys;                 // province-year noise

  // Reported case-fatality ratio (CFR) on logit (Lee–Carter)
  array[A] real beta_CFR_a;                   // baseline age-specific fatality
  array[R] real gamma_CFR_pr_base;            // period-1 anchor per region
  array[P-1, R] real gamma_CFR_pr_innov;      // RW innovations for periods 2..P
  ordered[A] raw_h_CFR_a;                     // monotone-increasing raw age-loading

  // Hyperparameters for kappa evolution
  real<lower=0> sigma_gamma_kappa;          // RW innovation SD (gamma_kappa_pr)
  real<lower=0> sigma_delta_kappa;          // SD for delta_kappa_ys
}

transformed parameters {
  // Center deltas by region (Se, FP), region-year (Fourier), region-year (kappa)
  array[S] real delta_Se_s;
  for (r_ind in 1:R) {
    real mean_r = 0; int cnt = 0;
    for (s_ind in 1:S) if (region_id[s_ind] == r_ind) { mean_r += raw_delta_Se_s[s_ind]; cnt += 1; }
    mean_r = mean_r / cnt;
    for (s_ind in 1:S) if (region_id[s_ind] == r_ind) delta_Se_s[s_ind] = raw_delta_Se_s[s_ind] - mean_r;
  }

  array[S] real delta_FP_s;
  for (r_ind in 1:R) {
    real mean_r = 0; int cnt = 0;
    for (s_ind in 1:S) if (region_id[s_ind] == r_ind) { mean_r += raw_delta_FP_s[s_ind]; cnt += 1; }
    mean_r = mean_r / cnt;
    for (s_ind in 1:S) if (region_id[s_ind] == r_ind) delta_FP_s[s_ind] = raw_delta_FP_s[s_ind] - mean_r;
  }

  array[H, Y, S] real delta_cos;
  array[H, Y, S] real delta_sin;
  for (h in 1:H) for (y_ind in 1:Y) {
    for (r_ind in 1:R) {
      real mean_cos = 0; real mean_sin = 0; int cnt = 0;
      for (s_ind in 1:S) if (region_id[s_ind] == r_ind) { mean_cos += raw_delta_cos[h, y_ind, s_ind]; mean_sin += raw_delta_sin[h, y_ind, s_ind]; cnt += 1; }
      mean_cos = mean_cos / cnt; mean_sin = mean_sin / cnt;
      for (s_ind in 1:S) if (region_id[s_ind] == r_ind) {
        delta_cos[h, y_ind, s_ind] = raw_delta_cos[h, y_ind, s_ind] - mean_cos;
        delta_sin[h, y_ind, s_ind] = raw_delta_sin[h, y_ind, s_ind] - mean_sin;
      }
    }
  }

  array[Y, S] real delta_kappa_ys;
  for (y_ind in 1:Y) for (r_ind in 1:R) {
    real mean_r = 0; int cnt = 0;
    for (s_ind in 1:S) if (region_id[s_ind] == r_ind) { mean_r += raw_delta_kappa_ys[y_ind, s_ind]; cnt += 1; }
    mean_r = mean_r / cnt;
    for (s_ind in 1:S) if (region_id[s_ind] == r_ind) delta_kappa_ys[y_ind, s_ind] = raw_delta_kappa_ys[y_ind, s_ind] - mean_r;
  }

  // Latent monthly totals (means for NB2 likelihoods)
  array[T, S] real<lower=0> mu_case_ts;
  array[T, S] real<lower=0> mu_death_ts;

  // Reused arrays
  array[Y, M, S] real<lower=0> kappa_month_yms;   // seasonal modifier
  array[Y, A, S] real<lower=0> kappa_age_yas;     // age modifier
  array[Y, S] real log_lambda_base;               // log baseline infection rate

  // Age sensitivity for Lee–Carter structures
  vector[A] h_kappa_a;
  array[P, R] real gamma_kappa_pr;

  // Se region–year baselines
  array[Y, R] real alpha_Se_yr;

  // CFR
  array[Y, A, S] real<lower=0, upper=1> CFR_yas;
  vector[A] h_CFR_a;
  array[P, R] real gamma_CFR_pr;

  // Age composition
  array[Y, S] vector[A] q_case_ys;
  array[Y, S] vector[A] q_death_ys;
  array[Y, A, S] real C_yas;
  array[Y, A, S] real D_yas;

  // FOI baseline: additive structure
  for (s_ind in 1:S) for (y_ind in 1:Y)
    log_lambda_base[y_ind, s_ind] = alpha_lambda + delta_lambda_s[s_ind] + delta_lambda_y[y_ind] + epsilon_lambda_ys[y_ind, s_ind];

  // Seasonality Fourier series
  for (y_ind in 1:Y) for (s_ind in 1:S) {
    vector[M] log_kappa_yms;
    for (m_ind in 1:M) {
      real theta = 2 * pi() * m_ind / M; real f = 0;
      for (h in 1:H)
        f += (a_cos[h, period_id[y_ind], region_id[s_ind]] + delta_cos[h, y_ind, s_ind]) * cos(h * theta)
           + (a_sin[h, period_id[y_ind], region_id[s_ind]] + delta_sin[h, y_ind, s_ind]) * sin(h * theta);
      log_kappa_yms[m_ind] = f;
    }
    real sum_exp_m = 0; for (m_ind in 1:M) sum_exp_m += exp(log_kappa_yms[m_ind]);
    for (m_ind in 1:M) kappa_month_yms[y_ind, m_ind, s_ind] = M * exp(log_kappa_yms[m_ind]) / sum_exp_m;
  }

  // Build alpha_Se_yr RW
  for (r_ind in 1:R) {
    alpha_Se_yr[1, r_ind] = alpha_Se_yr_base[r_ind];
    for (y_ind in 2:Y) alpha_Se_yr[y_ind, r_ind] = alpha_Se_yr[y_ind - 1, r_ind] + alpha_Se_yr_innov[y_ind - 1, r_ind];
  }

  // h_kappa_a (ordered -> center -> normalize)
  {
    vector[A] h_kappa_a_raw = raw_h_kappa_a;
    vector[A] h_kappa_a_centered = h_kappa_a_raw - rep_vector(mean(h_kappa_a_raw), A);
    real h_kappa_a_norm2 = dot_product(h_kappa_a_centered, h_kappa_a_centered);
    h_kappa_a = h_kappa_a_centered / sqrt(h_kappa_a_norm2 + 1e-12);
  }

  // gamma_kappa_pr RW
  for (r_ind in 1:R) {
    gamma_kappa_pr[1, r_ind] = gamma_kappa_pr_base[r_ind];
    for (p_ind in 2:P) gamma_kappa_pr[p_ind, r_ind] = gamma_kappa_pr[p_ind - 1, r_ind] + gamma_kappa_pr_innov[p_ind - 1, r_ind];
  }

  // kappa_age_yas Lee–Carter
  for (y_ind in 1:Y) for (s_ind in 1:S) {
    vector[A] log_kappa_age;
    for (a_ind in 1:A)
      log_kappa_age[a_ind] = beta_kappa_a[a_ind]
        + (gamma_kappa_pr[period_id[y_ind], region_id[s_ind]] + delta_kappa_ys[y_ind, s_ind]) * h_kappa_a[a_ind];
    real sum_exp_a = 0; for (a_ind in 1:A) sum_exp_a += exp(log_kappa_age[a_ind]);
    for (a_ind in 1:A) kappa_age_yas[y_ind, a_ind, s_ind] = A * exp(log_kappa_age[a_ind]) / sum_exp_a;
  }

  // h_CFR_a (ordered -> center -> normalize)
  {
    vector[A] h_CFR_a_raw = raw_h_CFR_a;
    vector[A] h_CFR_a_centered = h_CFR_a_raw - rep_vector(mean(h_CFR_a_raw), A);
    real h_CFR_a_norm2 = dot_product(h_CFR_a_centered, h_CFR_a_centered);
    h_CFR_a = h_CFR_a_centered / sqrt(h_CFR_a_norm2 + 1e-12);
  }

  // gamma_CFR_pr RW
  for (r_ind in 1:R) {
    gamma_CFR_pr[1, r_ind] = gamma_CFR_pr_base[r_ind];
    for (p_ind in 2:P) gamma_CFR_pr[p_ind, r_ind] = gamma_CFR_pr[p_ind - 1, r_ind] + gamma_CFR_pr_innov[p_ind - 1, r_ind];
  }

  // CFR_yas
  for (y_ind in 1:Y) for (s_ind in 1:S) for (a_ind in 1:A)
    CFR_yas[y_ind, a_ind, s_ind] = inv_logit(beta_CFR_a[a_ind] + gamma_CFR_pr[period_id[y_ind], region_id[s_ind]] * h_CFR_a[a_ind]);

  // Initialize aggregates
  for (y_ind in 1:Y) for (s_ind in 1:S) for (a_ind in 1:A) { C_yas[y_ind, a_ind, s_ind] = 0.0; D_yas[y_ind, a_ind, s_ind] = 0.0; }

  // Monthly totals and accumulate yearly aggregates
  for (s_ind in 1:S) for (t_ind in 1:T) {
    int y_ind = year_id[t_ind];
    int m_ind = month_id[t_ind];
    int r_ind = region_id[s_ind];

    real Se_local = inv_logit(alpha_Se_yr[y_ind, r_ind] + delta_Se_s[s_ind]);
    real mu_case_sum = 0; real mu_death_sum = 0;
    for (a_ind in 1:A) {
      real log_FP = alpha_FP_pr[r_ind, period_id[y_ind]] + delta_FP_s[s_ind] + delta_FP_y[y_ind] + eta_FP_ys[y_ind, s_ind];
      real FP_local = exp(log(N_tas[t_ind, a_ind, s_ind] + 1e-9) + log_FP);
      real lambda_local = exp(log_lambda_base[y_ind, s_ind]) * kappa_month_yms[y_ind, m_ind, s_ind] * kappa_age_yas[y_ind, a_ind, s_ind];
      real I_local = nu * S_tas[t_ind, a_ind, s_ind] * N_tas[t_ind, a_ind, s_ind] * (1 - exp(-lambda_local));
      real C_local = I_local * Se_local + FP_local;
      real D_local = C_local * CFR_yas[y_ind, a_ind, s_ind];
      mu_case_sum += C_local; mu_death_sum += D_local;
      C_yas[y_ind, a_ind, s_ind] += C_local; D_yas[y_ind, a_ind, s_ind] += D_local;
    }
    mu_case_ts[t_ind, s_ind] = mu_case_sum; mu_death_ts[t_ind, s_ind] = mu_death_sum;
  }

  // q_case_ys and q_death_ys
  for (s_ind in 1:S) for (y_ind in 1:Y) {
    vector[A] C_age; vector[A] D_age;
    for (a_ind in 1:A) { C_age[a_ind] = C_yas[y_ind, a_ind, s_ind]; D_age[a_ind] = D_yas[y_ind, a_ind, s_ind]; }
    if (sum(C_age) > 0) q_case_ys[y_ind, s_ind] = C_age / sum(C_age); else q_case_ys[y_ind, s_ind] = rep_vector(1.0 / A, A);
    if (sum(D_age) > 0) q_death_ys[y_ind, s_ind] = D_age / sum(D_age); else q_death_ys[y_ind, s_ind] = rep_vector(1.0 / A, A);
  }
}

model {
  // Dispersion priors
  phi_c ~ exponential(0.1);
  phi_d ~ exponential(200);

  // Clinical fraction prior
  nu ~ beta(1, 9);

  // FOI priors (additive baseline)
  alpha_lambda ~ normal(0, 1);
  delta_lambda_s ~ normal(0, 1);
  delta_lambda_y ~ normal(0, 1);
  sigma_lambda_noise ~ exponential(1);
  for (s_ind in 1:S) for (y_ind in 1:Y) epsilon_lambda_ys[y_ind, s_ind] ~ normal(0, sigma_lambda_noise);

  // Seasonality priors
  for (h in 1:H) for (p_ind in 1:P) for (r_ind in 1:R) { a_cos[h, p_ind, r_ind] ~ normal(0, 0.5); a_sin[h, p_ind, r_ind] ~ normal(0, 0.5); }
  for (h in 1:H) for (y_ind in 1:Y) for (s_ind in 1:S) { raw_delta_cos[h, y_ind, s_ind] ~ normal(0, 0.2); raw_delta_sin[h, y_ind, s_ind] ~ normal(0, 0.2); }

  // Kappa (Lee–Carter)
  beta_kappa_a ~ normal(0, 1);
  raw_h_kappa_a ~ normal(0, 1);
  sigma_gamma_kappa ~ exponential(1);
  sigma_delta_kappa ~ exponential(2);
  for (r_ind in 1:R) {
    gamma_kappa_pr_base[r_ind] ~ normal(0, 0.5);
    for (p_ind in 1:(P-1)) gamma_kappa_pr_innov[p_ind, r_ind] ~ normal(0, sigma_gamma_kappa);
  }
  for (s_ind in 1:S) for (y_ind in 1:Y) raw_delta_kappa_ys[y_ind, s_ind] ~ normal(0, sigma_delta_kappa);

  // Se RW
  sigma_alpha_Se ~ exponential(2);
  sigma_delta_Se ~ exponential(0.5);
  for (r_ind in 1:R) {
    alpha_Se_yr_base[r_ind] ~ normal(0, 0.5);
    for (y_ind in 1:(Y-1)) alpha_Se_yr_innov[y_ind, r_ind] ~ normal(0, sigma_alpha_Se);
  }
  raw_delta_Se_s ~ normal(0, sigma_delta_Se);

  // FP (alpha_pr + delta_s + delta_y + noise)
  for (r_ind in 1:R) for (p_ind in 1:P) alpha_FP_pr[r_ind, p_ind] ~ normal(0, 1);
  raw_delta_FP_s ~ normal(0, sigma_delta_FP);
  delta_FP_y ~ normal(0, 1);
  sigma_FP_noise ~ exponential(1);
  for (s_ind in 1:S) for (y_ind in 1:Y) eta_FP_ys[y_ind, s_ind] ~ normal(0, sigma_FP_noise);

  // CFR (Lee–Carter)
  beta_CFR_a ~ normal(0, 1);
  raw_h_CFR_a ~ normal(0, 1);
  for (r_ind in 1:R) {
    gamma_CFR_pr_base[r_ind] ~ normal(0, 0.5);
    for (p_ind in 1:(P-1)) gamma_CFR_pr_innov[p_ind, r_ind] ~ normal(0, sigma_gamma_CFR);
  }

  // Monthly totals likelihood (NB2)
  for (s_ind in 1:S) for (t_ind in 1:T) {
    R_case_ts[t_ind, s_ind] ~ neg_binomial_2(mu_case_ts[t_ind, s_ind], phi_c);
    R_death_ts[t_ind, s_ind] ~ neg_binomial_2(mu_death_ts[t_ind, s_ind], phi_d);
  }

  // Yearly age composition likelihood (multinomial)
  for (s_ind in 1:S) for (y_ind in 1:Y) {
    array[A] int R_case_age; array[A] int R_death_age; int N_case_ys = 0; int N_death_ys = 0;
    for (a_ind in 1:A) { R_case_age[a_ind] = R_case_yas[y_ind, a_ind, s_ind]; R_death_age[a_ind] = R_death_yas[y_ind, a_ind, s_ind]; N_case_ys += R_case_age[a_ind]; N_death_ys += R_death_age[a_ind]; }
    if (N_case_ys > 0) R_case_age ~ multinomial(q_case_ys[y_ind, s_ind]);
    if (N_death_ys > 0) R_death_age ~ multinomial(q_death_ys[y_ind, s_ind]);
  }
}

generated quantities {
  array[T, S] real log_lik_case;
  array[T, S] real log_lik_death;
  array[Y, S] real log_lik_case_age;
  array[Y, S] real log_lik_death_age;

  for (s_ind in 1:S) for (t_ind in 1:T) {
    log_lik_case[t_ind, s_ind] = neg_binomial_2_lpmf(R_case_ts[t_ind, s_ind] | mu_case_ts[t_ind, s_ind], phi_c);
    log_lik_death[t_ind, s_ind] = neg_binomial_2_lpmf(R_death_ts[t_ind, s_ind] | mu_death_ts[t_ind, s_ind], phi_d);
  }
  for (s_ind in 1:S) for (y_ind in 1:Y) {
    array[A] int R_case_age; array[A] int R_death_age; int N_case_ys = 0; int N_death_ys = 0;
    for (a_ind in 1:A) { R_case_age[a_ind] = R_case_yas[y_ind, a_ind, s_ind]; R_death_age[a_ind] = R_death_yas[y_ind, a_ind, s_ind]; N_case_ys += R_case_age[a_ind]; N_death_ys += R_death_age[a_ind]; }
    if (N_case_ys > 0) log_lik_case_age[y_ind, s_ind] = multinomial_lpmf(R_case_age | q_case_ys[y_ind, s_ind]); else log_lik_case_age[y_ind, s_ind] = 0;
    if (N_death_ys > 0) log_lik_death_age[y_ind, s_ind] = multinomial_lpmf(R_death_age | q_death_ys[y_ind, s_ind]); else log_lik_death_age[y_ind, s_ind] = 0;
  }
}
