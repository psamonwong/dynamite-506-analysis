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
  array[Y] int<lower=1, upper=P> period_id;       // p(y(t)): period block corresponding to month t
  array[S] int<lower=1, upper=R> region_id;       // r(s): region corresponding to province s

  // Population size
  array[T,A,S] int<lower=0> N_tas;          // N_t,a,s population size at time t, age a, province s
  // Susceptible proportions
  array[T, A, S] real<lower=0, upper=1> S_tas;   // S_t,a,s proportion susceptible at time t, age a, province s

  // Reported cases and deaths
  array[T, S] int<lower=0> R_case_ts;          // R_case_t,s reported cases at time t, province s
  array[T, S] int<lower=0> R_death_ts;         // R_death_t,s reported deaths at time t, province s
  array[Y, A, S] int<lower=0> R_case_yas;      // R_case_y,a,s reported cases at year y, age a, province s
  array[Y, A, S] int<lower=0> R_death_yas;     // R_death_y,a,s reported deaths at year y, age a, province s

  // Seasonality configurations
  int<lower=1> H;                 // number of Fourier harmonics
}

transformed data {
    real nu = 0.1;                // fixed clinical fraction (10%)
    real beta_FP = 0.0;           // fixed: no transmission effect on FP (was unidentified)
}

parameters{
    // Overdispersion parameters for monthly totals (NB2)
    real<lower=1e-09> phi_c;           // dispersion for reported cases totals
    real<lower=1e-09> phi_d;           // dispersion for reported deaths totals

    // Baseline FOI: regional intercept + province deviation
    array[R] real alpha_lambda_r;          // region-specific baseline log-FOI
    array[S] real alpha_lambda_s;          // province deviations from region
    // Seasonality (Fourier): region–period baseline only
    array[H, P, R] real a_cos;        // cosine coefficients by harmonic/period/region
    array[H, P, R] real a_sin;        // sine coefficients by harmonic/period/region
    // Age transmission modifier (sum-to-zero age effects)
    vector[A-1] alpha_kappa_a;                   // age effects (ages 1..A-1; age A constrained to -sum)
    // Period trend (random walk)
    array[R] real<lower=0> gamma_kappa_pr_base;   // period-1 anchor (strictly positive) per region
    array[P-1, R] real gamma_kappa_pr_innov;      // RW innovations for periods 2..P

    // Case detection sensitivity (Se) on logit scale: region + province + year (additive)
    array[R] real alpha_Se_r;                   // regional baseline (logit scale)
    array[S] real delta_Se_s;                   // province deviations (i.i.d.)
    array[Y] real delta_Se_y;                   // year deviations (i.i.d.)
    real<lower=0> sigma_delta_Se;               // SD for delta_Se_s
    real<lower=0> sigma_delta_Se_y;             // SD for delta_Se_y

    // False positives (background reporting): region + province + year (additive on log scale)
    array[R] real alpha_FP_r;                   // regional baseline (log scale)
    array[S] real delta_FP_s;                   // province deviations (i.i.d.)
    array[Y] real delta_FP_y;                   // year deviations (i.i.d.)
    real<lower=0> sigma_delta_FP;               // SD for delta_FP_s
    real<lower=0> sigma_delta_FP_y;             // SD for delta_FP_y

    // Reported case-fatality ratio (CFR) (sum-to-zero age effects)
    vector[A-1] alpha_CFR_a;                     // age effects (ages 1..A-1; age A constrained to -sum)
    // Period trend (random walk)
    array[R] real<lower=0> gamma_CFR_pr_base;     // period-1 anchor (strictly positive) per region
    array[P-1, R] real gamma_CFR_pr_innov;        // RW innovations for periods 2..P

    // Hyperparameters for kappa evolution
    real<lower=0> sigma_gamma_kappa;          // RW innovation SD (gamma_kappa_pr)

    // Hyperparameters for CFR evolution
    real<lower=0> sigma_gamma_CFR;            // RW innovation SD (gamma_CFR_pr)

    // Hyperparameters for seasonality (Fourier)
    real<lower=0> sigma_seasonality_pr;       // SD for region-period baseline Fourier coefficients (a_cos, a_sin)

    // Hyperparameters for baseline FOI
    real<lower=0> sigma_lambda_r;              // SD for regional baseline FOI
    real<lower=0> sigma_lambda_s;              // SD for province deviations from region
}

transformed parameters {
    // Latent monthly totals (means for NB2 likelihoods) — KEPT for likelihood
    array[T, S] real<lower=0> mu_case_ts;        // mu_case_t,s mean reported cases at time t, province s
    array[T, S] real<lower=0> mu_death_ts;       // mu_death_t,s mean reported deaths at time t, province s
    
    // Smaller arrays that are reused — KEPT
    array[Y, M, S] real<lower=0> kappa_month_yms;   // kappa_month_y,m,s seasonal modifier
    array[Y, A, S] real<lower=0> kappa_age_yas;     // kappa_age_y,a,s age modifier
    array[Y, S] real log_lambda_base;               // log baseline infection rate for year y, province s

    // Age effects and period trends
    vector[A] alpha_kappa_age;                 // sum-to-zero age effects (computed from alpha_kappa_a)
    array[P, R] real gamma_kappa_pr;           // period trend (random walk)
    
    // Case fatality rate — KEPT
    array[Y, A, S] real<lower=0, upper=1> CFR_yas;    // CFR_y,a,s case fatality rate

    // Age effects and period trends for CFR
    vector[A] alpha_CFR_age;                   // sum-to-zero age effects (computed from alpha_CFR_a)
    array[P, R] real gamma_CFR_pr;             // period trend (random walk)

    // Proportion of cases and deaths by age — KEPT for multinomial likelihood
    array[Y, S] vector[A] q_case_ys;  // q_case_y,s: proportion of cases by age
    array[Y, S] vector[A] q_death_ys;  // q_death_y,s: proportion of deaths by age
    
    // Yearly aggregates for age composition (computed inline)
    array[Y, A, S] real C_yas;  // yearly cases by age (aggregated from monthly)
    array[Y, A, S] real D_yas;  // yearly deaths by age (aggregated from monthly)

    // Baseline FOI: regional intercept + province deviation
    for(s_ind in 1:S){
        for(y_ind in 1:Y){
                log_lambda_base[y_ind, s_ind] = alpha_lambda_r[region_id[s_ind]] + 
                                                 alpha_lambda_s[s_ind];
        }
    }

    // Seasonality: Monthly Fourier series for kappa_month_yms
    for(y_ind in 1:Y){
        for(s_ind in 1:S){
            vector[M] log_kappa_yms; // log seasonal modifier for month m, year y, province s
            for(m_ind in 1:M){
                real theta = 2 * pi() * m_ind / M; // angle for month m
                real f = 0;
                for(h in 1:H){
                    // Fourier series expansion for month m with harmonics h
                    f += a_cos[h, period_id[y_ind], region_id[s_ind]] * cos(h * theta) +
                         a_sin[h, period_id[y_ind], region_id[s_ind]] * sin(h * theta);
                }
                log_kappa_yms[m_ind] = f;
            }
            // Normalize across months to mean 1 (identifiability)
            real sum_exp_m = 0; 
            for(m_ind in 1:M) sum_exp_m += exp(log_kappa_yms[m_ind]); // sum over months 
            for(m_ind in 1:M) kappa_month_yms[y_ind, m_ind, s_ind] = M * exp(log_kappa_yms[m_ind])/ sum_exp_m; // normalize to mean 1
        }
    }

    // Sum-to-zero age effects for kappa
    {
    real sum_alpha = sum(alpha_kappa_a);
    for (a_ind in 1:(A-1)) {
        alpha_kappa_age[a_ind] = alpha_kappa_a[a_ind];
    }
    alpha_kappa_age[A] = -sum_alpha;  // last age constrained so sum = 0
    }

    // Build gamma_kappa_pr via anchored random walk (period 1 positive)
    for (r_ind in 1:R) {
        gamma_kappa_pr[1, r_ind] = gamma_kappa_pr_base[r_ind];
        for (p_ind in 2:P) {
            gamma_kappa_pr[p_ind, r_ind] = gamma_kappa_pr[p_ind - 1, r_ind] + gamma_kappa_pr_innov[p_ind - 1, r_ind];
        }
    }


    // Age modifier: sum-to-zero age effects + period trend, normalized to mean 1
    for(y_ind in 1:Y){
        for(s_ind in 1:S){
            vector[A] log_kappa_age; // log age modifier for age a, year y, province s
            for(a_ind in 1:A){
                log_kappa_age[a_ind] = alpha_kappa_age[a_ind] + 
                                       gamma_kappa_pr[period_id[y_ind], region_id[s_ind]] * alpha_kappa_age[a_ind];
            }
            // Normalize across ages to mean 1 (identifiability)
            real sum_exp_a = 0;
            for(a_ind in 1:A) sum_exp_a += exp(log_kappa_age[a_ind]); // sum over ages
            for(a_ind in 1:A) kappa_age_yas[y_ind, a_ind, s_ind] = A * exp(log_kappa_age[a_ind]) / sum_exp_a; // normalize to mean 1
        }
    }

    {
    // Sum-to-zero age effects for CFR
    real sum_alpha = sum(alpha_CFR_a);
    for (a_ind in 1:(A-1)) {
        alpha_CFR_age[a_ind] = alpha_CFR_a[a_ind];
    }
    alpha_CFR_age[A] = -sum_alpha;  // last age constrained so sum = 0
    }

    // Build gamma_CFR_pr via anchored random walk (period 1 positive)
    for (r_ind in 1:R) {
        gamma_CFR_pr[1, r_ind] = gamma_CFR_pr_base[r_ind];
        for (p_ind in 2:P) {
            gamma_CFR_pr[p_ind, r_ind] = gamma_CFR_pr[p_ind - 1, r_ind] + gamma_CFR_pr_innov[p_ind - 1, r_ind];
        }
    }


    // Case fatality rate CFR_yas (sum-to-zero age effects + period trend on logit)
    for(y_ind in 1:Y){
        for(s_ind in 1:S){
            for(a_ind in 1:A){
                CFR_yas[y_ind, a_ind, s_ind] = inv_logit(
                    alpha_CFR_age[a_ind] + 
                    gamma_CFR_pr[period_id[y_ind], region_id[s_ind]] * alpha_CFR_age[a_ind]
                );
            }
        }
    }
    
    // Initialize yearly aggregates for age composition
    for(y_ind in 1:Y){
        for(s_ind in 1:S){
            for(a_ind in 1:A){
                C_yas[y_ind, a_ind, s_ind] = 0.0;
                D_yas[y_ind, a_ind, s_ind] = 0.0;
            }
        }
    }
    
    // Compute monthly totals and accumulate yearly aggregates
    // (lambda_tas, I_tas, C_tas, D_tas, Se_ts, FP_tas computed inline as locals)
    for(s_ind in 1:S){
        for(t_ind in 1:T){
            int y_ind = year_id[t_ind];
            int m_ind = month_id[t_ind];
            int r_ind = region_id[s_ind];
            
            // Case detection sensitivity: region + province + year (additive on logit scale)
            real Se_local = inv_logit(alpha_Se_r[r_ind] + delta_Se_s[s_ind] + delta_Se_y[y_ind]);
            
            real mu_case_sum = 0;
            real mu_death_sum = 0;
            for(a_ind in 1:A){
                // False positives: population offset + baseline transmission effect + region + province + year
                real FP_local = exp(log(N_tas[t_ind, a_ind, s_ind] + 1e-9) +
                                    beta_FP * log_lambda_base[y_ind, s_ind] +
                                    alpha_FP_r[r_ind] + delta_FP_s[s_ind] + delta_FP_y[y_ind]);
                
                // Force of infection (local, was lambda_tas)
                real lambda_local = exp(log_lambda_base[y_ind, s_ind]) *
                                    kappa_month_yms[y_ind, m_ind, s_ind] *
                                    kappa_age_yas[y_ind, a_ind, s_ind];
                
                // Expected new clinical infections (local, was I_tas)
                real I_local = nu *
                               S_tas[t_ind, a_ind, s_ind] *
                               N_tas[t_ind, a_ind, s_ind] *
                               (1 - exp(-lambda_local));
                
                // Expected reported cases (local, was C_tas)
                real C_local = I_local * Se_local + FP_local;
                
                // Expected reported deaths (local, was D_tas)
                real D_local = C_local * CFR_yas[y_ind, a_ind, s_ind];
                
                // Aggregate to monthly totals
                mu_case_sum += C_local;
                mu_death_sum += D_local;
                
                // Accumulate to yearly totals for age composition
                C_yas[y_ind, a_ind, s_ind] += C_local;
                D_yas[y_ind, a_ind, s_ind] += D_local;
            }
            mu_case_ts[t_ind, s_ind] = mu_case_sum;
            mu_death_ts[t_ind, s_ind] = mu_death_sum;
        }
    }

    // Proportion of cases and deaths by age q_case_ys, q_death_ys 
    for(s_ind in 1:S){
        for(y_ind in 1:Y){
            vector[A] C_age;
            vector[A] D_age;
            for(a_ind in 1:A){
                C_age[a_ind] = C_yas[y_ind, a_ind, s_ind];
                D_age[a_ind] = D_yas[y_ind, a_ind, s_ind];
            }
            if(sum(C_age) > 0)
                    q_case_ys[y_ind, s_ind] = C_age / sum(C_age);
                else
                    q_case_ys[y_ind, s_ind] = rep_vector(1.0 / A, A); // uniform if no cases
                if(sum(D_age) > 0)
                    q_death_ys[y_ind, s_ind] = D_age / sum(D_age);
                else
                    q_death_ys[y_ind, s_ind] = rep_vector(1.0 / A, A); // uniform if no deaths
        }
    }

}

model {
    sigma_lambda_r ~ exponential(1); // prior for regional baseline FOI SD
    sigma_lambda_s ~ exponential(2); // prior for province deviation FOI SD

    // Priors for dispersion of monthly totals
    phi_c ~ normal(0, 5);          // cases NB2 dispersion
    phi_d ~ normal(0, 5);          // deaths NB2 dispersion

    // Baseline FOI priors (hierarchical: region -> province -> year)
    alpha_lambda_r ~ normal(0, sigma_lambda_r); // regional baseline intercepts
    for (s_ind in 1:S) {
        alpha_lambda_s[s_ind] ~ normal(0, sigma_lambda_s); // province deviations from region
    }

    // Se/FP additive decomposition priors
    for (r_ind in 1:R) {
        alpha_Se_r[r_ind] ~ normal(0, 1);   // regional baseline sensitivity
        alpha_FP_r[r_ind] ~ normal(0, 1);   // regional baseline false positives
    }
    for (s_ind in 1:S) {
        delta_Se_s[s_ind] ~ normal(0, sigma_delta_Se);   // province deviations
        delta_FP_s[s_ind] ~ normal(0, sigma_delta_FP);
    }
    for (y_ind in 1:Y) {
        delta_Se_y[y_ind] ~ normal(0, sigma_delta_Se_y); // year deviations
        delta_FP_y[y_ind] ~ normal(0, sigma_delta_FP_y);
    }
    
    // Hyperpriors for Se/FP
    sigma_delta_Se ~ exponential(2);     // mean = 0.5 (province heterogeneity)
    sigma_delta_Se_y ~ exponential(2);   // mean = 0.5 (year deviations)
    sigma_delta_FP ~ exponential(2);     // mean = 0.5
    sigma_delta_FP_y ~ exponential(2);   // mean = 0.5

    // Kappa (sum-to-zero age effects) priors
    alpha_kappa_a ~ normal(0, 2);       // age effects (sum to zero via constraint)
    sigma_gamma_kappa ~ exponential(1); // RW innovation SD

    // CFR (sum-to-zero age effects) priors
    alpha_CFR_a ~ normal(0, 2);         // age effects (sum to zero via constraint)
    sigma_gamma_CFR ~ exponential(1);   // RW innovation SD (analogous to kappa)

    // Seasonality hyperpriors
    sigma_seasonality_pr ~ exponential(2);  // mean = 0.5 for region-period baseline

    // Prior for seasonal Fourier coefficients (region–period baseline)
    for(h in 1:H){
        for(p_ind in 1:P){
            for(r_ind in 1:R){
                a_cos[h, p_ind, r_ind] ~ normal(0, sigma_seasonality_pr);
                a_sin[h, p_ind, r_ind] ~ normal(0, sigma_seasonality_pr);
            }
        }
    }

    // Region–period anchored random walk for gamma_kappa_pr (kappa)
    for (r_ind in 1:R){
        gamma_kappa_pr_base[r_ind] ~ normal(0, 0.5); // Prior for first period (half-normal via <lower=0>)
        for (p_ind in 1:(P-1)){
            gamma_kappa_pr_innov[p_ind, r_ind] ~ normal(0, sigma_gamma_kappa); // RW innovations for periods 2..P
        }
    }

    // Region–period anchored random walk for gamma_CFR_pr (CFR)
    for (r_ind in 1:R){
        gamma_CFR_pr_base[r_ind] ~ normal(0, 0.5); // Prior for first period (half-normal via <lower=0>)
        for (p_ind in 1:(P-1)){
            gamma_CFR_pr_innov[p_ind, r_ind] ~ normal(0, sigma_gamma_CFR); // RW innovations for periods 2..P
        }
    }

    // Monthly totals likelihood (NB2)
    for(s_ind in 1:S){
        for(t_ind in 1:T){
        R_case_ts[t_ind, s_ind] ~ neg_binomial_2(mu_case_ts[t_ind, s_ind], phi_c);
        R_death_ts[t_ind, s_ind] ~ neg_binomial_2(mu_death_ts[t_ind, s_ind], phi_d);
        }
    }

    // Yearly age composition likelihood (multinomial)
    for(s_ind in 1:S){
        for(y_ind in 1:Y){
            array[A] int R_case_age; // reported cases by age for year y, province s
            array[A] int R_death_age; // reported deaths by age for year y, province s
            int N_case_ys = 0; // total reported cases for year y, province s
            int N_death_ys = 0; // total reported deaths for year y, province s
            for(a_ind in 1:A){
                R_case_age[a_ind] = R_case_yas[y_ind, a_ind, s_ind]; // reported cases by age
                R_death_age[a_ind] = R_death_yas[y_ind, a_ind, s_ind]; // reported deaths by age 
                N_case_ys += R_case_age[a_ind]; // total reported cases
                N_death_ys += R_death_age[a_ind]; // total reported deaths
            }
            if(N_case_ys > 0)
                R_case_age ~ multinomial(q_case_ys[y_ind, s_ind]);
            if(N_death_ys > 0)
                R_death_age ~ multinomial(q_death_ys[y_ind, s_ind]);
        }
    }
}


generated quantities {
    // Log-likelihood for leave-one-out cross-validation (loo package)
    // Pointwise log-likelihood for each observation
  
    array[T, S] real log_lik_case;    // log-likelihood for monthly cases
    array[T, S] real log_lik_death;   // log-likelihood for monthly deaths
    array[Y, S] real log_lik_case_age;  // log-likelihood for yearly age composition (cases)
    array[Y, S] real log_lik_death_age; // log-likelihood for yearly age composition (deaths)

    // Compute log-likelihood for monthly totals
    for(s_ind in 1:S){
        for(t_ind in 1:T){
            log_lik_case[t_ind, s_ind] = neg_binomial_2_lpmf(R_case_ts[t_ind, s_ind] | mu_case_ts[t_ind, s_ind], phi_c);
            log_lik_death[t_ind, s_ind] = neg_binomial_2_lpmf(R_death_ts[t_ind, s_ind] | mu_death_ts[t_ind, s_ind], phi_d);
        }
    }
  
    // Compute log-likelihood for yearly age composition
    for(s_ind in 1:S){
        for(y_ind in 1:Y){
            array[A] int R_case_age; 
            array[A] int R_death_age;
            int N_case_ys = 0;
            int N_death_ys = 0;
      
            for(a_ind in 1:A){
                R_case_age[a_ind] = R_case_yas[y_ind, a_ind, s_ind];
                R_death_age[a_ind] = R_death_yas[y_ind, a_ind, s_ind];
                N_case_ys += R_case_age[a_ind];
                N_death_ys += R_death_age[a_ind];
            }
      
            if(N_case_ys > 0)
                log_lik_case_age[y_ind, s_ind] = multinomial_lpmf(R_case_age | q_case_ys[y_ind, s_ind]);
            else
                log_lik_case_age[y_ind, s_ind] = 0;
        
            if(N_death_ys > 0)
                log_lik_death_age[y_ind, s_ind] = multinomial_lpmf(R_death_age | q_death_ys[y_ind, s_ind]);
            else
                log_lik_death_age[y_ind, s_ind] = 0;
        }
    }
}


