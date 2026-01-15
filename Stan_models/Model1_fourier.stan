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
    
}

parameters{
    // Overdispersion parameters for monthly totals (NB2)
    real<lower=0> phi_c;           // dispersion for reported cases totals
    real<lower=0> phi_d;           // dispersion for reported deaths totals

    // Baseline FOI: province intercept + yearly random effect
    array[S] real alpha_lambda_s;          // province-specific baseline log-FOI
    real<lower=0> sigma_epsilon_lambda;    // SD for yearly residuals
    array[Y, S] real epsilon_lambda_ys;    // yearly residuals per (y,s)
    // Seasonality (Fourier): region–period baseline + year–province deviations
    array[H, P, R] real a_cos;        // cosine coefficients by harmonic/period/region
    array[H, P, R] real a_sin;        // sine coefficients by harmonic/period/region
    array[H, Y, S] real delta_cos;    // year–province cosine deviations
    array[H, Y, S] real delta_sin;    // year–province sine deviations
    // Age transmission modifier (Lee–Carter)
    array[A] real beta_kappa_a;                  // baseline age-specific RR
    // Sign-anchored random-walk parameterization for gamma_kappa_pr
    array[R] real gamma_kappa_pr_base;   // period-1 anchor (positive) per region
    array[P-1, R] real gamma_kappa_pr_innov;      // RW innovations for periods 2..P
    array[Y, S] real delta_kappa_ys;             // province–year deviation
    ordered[A] raw_h_kappa_a;                    // monotone-increasing raw age-loading

    // Clinical fraction by age (probability an infection is clinical)
    array[A] real<lower=0, upper=1> nu_a;       // age-specific clinical fraction

    // Case detection sensitivity (Se) on logit scale
    array[Y, R] real alpha_Se_yr;               // region–year baseline (random walk)
    array[S] real delta_Se_s;                   // province deviations (i.i.d.)
    real<lower=0> sigma_alpha_Se;               // RW innovation SD for alpha_Se_yr
    real<lower=0> sigma_delta_Se;               // SD for delta_Se_s

    // False positives (background reporting)
    array[Y, R] real alpha_FP_yr;               // region–year baseline (random walk)
    array[S] real delta_FP_s;                   // province deviations (i.i.d.)
    real<lower=0> sigma_alpha_FP;               // RW innovation SD for alpha_FP_yr 
    real<lower=0> sigma_delta_FP;               // SD for delta_FP_s

    // Reported case-fatality ratio (CFR) on logit (Lee–Carter)
    array[A] real beta_CFR_a;                   // baseline age-specific fatality
    // Sign-anchored random-walk parameterization for gamma_CFR_pr
    array[R] real gamma_CFR_pr_base;     // period-1 anchor (positive) per region
    array[P-1, R] real gamma_CFR_pr_innov;        // RW innovations for periods 2..P
    ordered[A] raw_h_CFR_a;                     // monotone-increasing raw age-loading
    real<lower=0> sigma_gamma_CFR;              // RW innovation SD for gamma_CFR_pr

    // Hyperparameters for kappa evolution
    real<lower=0> sigma_gamma_kappa;          // RW innovation SD (gamma_kappa_pr)
    real<lower=0> sigma_delta_kappa;          // SD for delta_kappa_ys

    // Hyperparameter for baseline FOI
    real<lower=0> sigma_lambda_s;              // SD for province baseline FOI
}

transformed parameters {
    // Latent monthly totals (means for NB2 likelihoods)
    array[T, S] real<lower=0> mu_case_ts;        // mu_case_t,s mean reported cases at time t, province s
    array[T, S] real<lower=0> mu_death_ts;       // mu_death_t,s mean reported deaths at time t, province s
    
    // Infection rate and related quantities
    array[T, A, S] real<lower=0> lambda_tas;    // lambda_t,a,s infection rate at time t, age a, province s
    array[Y, M, S] real<lower=0> kappa_month_yms;   // kappa_month_y,m,s seasonal modifier for month m, year y, province s
    array[Y, A, S] real<lower=0> kappa_age_yas;     // kappa_age_y,a,s age modifier for age a, year y, province s

    // Transmission rate parameters
    array[Y, S] real log_lambda_base;        // log baseline infection rate for year y, province s

    // Age sensitivity for lee–carter structures
    vector[A] h_kappa_a;                       // h_kappa[a]: age sensitivity for age a (oriented, final parameter)
    array[P, R] real gamma_kappa_pr;           // oriented gamma for product (final parameter)

    // Expected new clinical infections and reported quantities
    array[T, A, S] real<lower=0> I_tas;      // I_t,a,s expected new clinical infections at time t, age a, province s
    // Expected reported cases and deaths
    array[T, A, S] real<lower=0> C_tas;      // C_t,a,s expected reported cases at time t, age a, province s
    array[T, A, S] real<lower=0> D_tas;      // D_t,a,s expected reported deaths at time t, age a, province s

    // Case detection sensitivity and false positives
    array[T, S] real<lower=0, upper=1> Se_ts;        // Se_t,s case detection sensitivity at time t, province s
    array[T, A, S] real<lower=0> FP_tas;            // FP_t,a,s false positives at time t, age a, province s 
    
    // Case fatality rate
    array[Y, A, S] real<lower=0, upper=1> CFR_yas;    // CFR_y,a,s case fatality rate at year y, age a, province s

    // Age sensitivity for lee–carter structures
    vector[A] h_CFR_a;                       // h_CFR[a]: age sensitivity for age a (oriented, final parameter)
    array[P, R] real gamma_CFR_pr;             // oriented gamma for product (final parameter)

    // Proportion of cases and deaths by age (yearly age composition)
    array[Y, S] vector[A] q_case_ys;  // q_case_y,s: proportion of cases by age a, year y, province s
    array[Y, S] vector[A] q_death_ys;  // q_death_y,s: proportion of deaths by age a, year y, province s

    // Baseline FOI: province intercept + yearly residual
    for(s_ind in 1:S){
        for(y_ind in 1:Y){
                log_lambda_base[y_ind, s_ind] = alpha_lambda_s[s_ind] + epsilon_lambda_ys[y_ind, s_ind];
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
                    f += (a_cos[h, period_id[y_ind], region_id[s_ind]] + delta_cos[h,y_ind,s_ind]) * cos(h * theta) +
                        (a_sin[h, period_id[y_ind], region_id[s_ind]] + delta_sin[h,y_ind,s_ind]) * sin(h * theta);
                }
                log_kappa_yms[m_ind] = f;
            }
            // Normalize across months to mean 1 (identifiability)
            real sum_exp_m = 0; 
            for(m_ind in 1:M) sum_exp_m += exp(log_kappa_yms[m_ind]); // sum over months 
            for(m_ind in 1:M) kappa_month_yms[y_ind, m_ind, s_ind] = M * exp(log_kappa_yms[m_ind])/ sum_exp_m; // normalize to mean 1
        }
    }

    // Age sensitivity h_kappa_a (ordered -> center -> L2-normalize)
    {
    vector[A] h_kappa_a_raw = raw_h_kappa_a; // `ordered[A]` ensures monotone increase
    vector[A] h_kappa_a_centered = h_kappa_a_raw - rep_vector(mean(h_kappa_a_raw), A);
    real h_kappa_a_norm2 = dot_product(h_kappa_a_centered, h_kappa_a_centered);  // sum(centered^2)
    h_kappa_a = h_kappa_a_centered / sqrt(h_kappa_a_norm2 + 1e-12); // preserves monotonicity
    }

    // Build gamma_kappa_pr via anchored random walk (period 1 positive)
    for (r_ind in 1:R) {
        gamma_kappa_pr[1, r_ind] = gamma_kappa_pr_base[r_ind];
        for (p_ind in 2:P) {
            gamma_kappa_pr[p_ind, r_ind] = gamma_kappa_pr[p_ind - 1, r_ind] + gamma_kappa_pr_innov[p_ind - 1, r_ind];
        }
    }


    // Age modifier: kappa_age_yas Lee–Carter on log scale, normalized to mean 1 over age
    for(y_ind in 1:Y){
        for(s_ind in 1:S){
            vector[A] log_kappa_age; // log age modifier for age a, year y, province s
            for(a_ind in 1:A){
                log_kappa_age[a_ind] =
                    beta_kappa_a[a_ind] + // age effect 
                    (gamma_kappa_pr[period_id[y_ind], region_id[s_ind]] + // region–period interaction
                     delta_kappa_ys[y_ind, s_ind]) * // yearly deviation for year y, province s
                     h_kappa_a[a_ind]; // age sensitivity for age a
            }
            // Normalize across ages to mean 1 (identifiability)
            real sum_exp_a = 0;
            for(a_ind in 1:A) sum_exp_a += exp(log_kappa_age[a_ind]); // sum over ages
            for(a_ind in 1:A) kappa_age_yas[y_ind, a_ind, s_ind] = A * exp(log_kappa_age[a_ind]) / sum_exp_a; // normalize to mean 1
        }
    }

    // Case detection sensitivity Se_ts (logit)
    for(s_ind in 1:S){
        for(t_ind in 1:T){
            Se_ts[t_ind, s_ind] = 
                inv_logit(alpha_Se_yr[year_id[t_ind], region_id[s_ind]] + // baseline sensitivity by year and region
                          delta_Se_s[s_ind]); // province deviation
        }
    }
    // False positives FP_tas (age-specific scaling via log(N_tas[t,a,s]); use total pop if desired)
    for(s_ind in 1:S){
        for(t_ind in 1:T){
            for(a_ind in 1:A){
                FP_tas[t_ind, a_ind, s_ind] = exp(log(N_tas[t_ind, a_ind, s_ind]+ 1e-9) +
                                                alpha_FP_yr[year_id[t_ind], region_id[s_ind]] +
                                                delta_FP_s[s_ind]);
                                                }
        }
    }

    
    {
    vector[A] h_CFR_a_raw = raw_h_CFR_a; // ordered ascending
    vector[A] h_CFR_a_centered = h_CFR_a_raw - rep_vector(mean(h_CFR_a_raw), A);
    real h_CFR_a_norm2 = dot_product(h_CFR_a_centered, h_CFR_a_centered);
    h_CFR_a = h_CFR_a_centered / sqrt(h_CFR_a_norm2 + 1e-12);
    }

    // Build gamma_CFR_pr via anchored random walk (period 1 positive)
    for (r_ind in 1:R) {
        gamma_CFR_pr[1, r_ind] = gamma_CFR_pr_base[r_ind];
        for (p_ind in 2:P) {
            gamma_CFR_pr[p_ind, r_ind] = gamma_CFR_pr[p_ind - 1, r_ind] + gamma_CFR_pr_innov[p_ind - 1, r_ind];
        }
    }


    // Case fatality rate CFR_yas (Lee–Carter on logit)
    for(y_ind in 1:Y){
        for(s_ind in 1:S){
            for(a_ind in 1:A){
                CFR_yas[y_ind, a_ind, s_ind] = inv_logit(
                    beta_CFR_a[a_ind] + // age effect
                    gamma_CFR_pr[period_id[y_ind], region_id[s_ind]] * // period-region trend
                    h_CFR_a[a_ind] // age sensitivity
                );
            }
        }
    }
    for(s_ind in 1:S){
        for(t_ind in 1:T){
            real mu_case_sum = 0;
            real mu_death_sum = 0;
            for(a_ind in 1:A){
                // Force of infection at time t, age a, province s 
                lambda_tas[t_ind, a_ind, s_ind] = exp(log_lambda_base[year_id[t_ind], s_ind]) * // baseline infection rate
                                                    kappa_month_yms[year_id[t_ind], month_id[t_ind], s_ind] * // seasonal modifier
                                                    kappa_age_yas[year_id[t_ind], a_ind, s_ind]; // age modifier
                // Expected number of new clinical infections at time t, age a, province s
                I_tas[t_ind, a_ind, s_ind] = nu_a[a_ind] * // clinical fraction by age a 
                        S_tas[t_ind, a_ind, s_ind] *  // proportion susceptible at time t, age a, province s
                        N_tas[t_ind, a_ind, s_ind] *  // population size at time t, age a, province s
                        (1 - exp(-lambda_tas[t_ind, a_ind, s_ind])); // infection probability at time t, age a, province s
                // Expected reported cases and deaths at time t, age a, province s
                C_tas[t_ind, a_ind, s_ind] = (I_tas[t_ind, a_ind, s_ind] * Se_ts[t_ind, s_ind]) + // Case detection sensitivity
                                                 FP_tas[t_ind, a_ind, s_ind]; // False positives
                D_tas[t_ind, a_ind, s_ind] = C_tas[t_ind, a_ind, s_ind] * // Calculated from reported cases (not infections)
                                                CFR_yas[year_id[t_ind], a_ind, s_ind]; // case fatality rate by year, age, province 
                // Aggregate to monthly totals for likelihood (over age)
                mu_case_sum += C_tas[t_ind, a_ind, s_ind];
                mu_death_sum += D_tas[t_ind, a_ind, s_ind];
            }
            mu_case_ts[t_ind, s_ind] = mu_case_sum;
            mu_death_ts[t_ind, s_ind] = mu_death_sum;
        }
    }

    // Proportion of cases and deaths by age q_case_ys, q_death_ys 
    for(s_ind in 1:S){
        for(y_ind in 1:Y){
            vector[A] C_age = rep_vector(0.0, A);
            vector[A] D_age = rep_vector(0.0, A);
            for(t_ind in 1:T){
                if(year_id[t_ind] == y_ind){
                    for(a_ind in 1:A){
                        C_age[a_ind] += C_tas[t_ind, a_ind, s_ind]; 
                        D_age[a_ind] += D_tas[t_ind, a_ind, s_ind];
                    }
                
                }
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
    sigma_lambda_s ~ exponential(1); // prior for province baseline FOI SD
    sigma_epsilon_lambda ~ exponential(1); // prior for yearly FOI residual SD

    // Priors for dispersion of monthly totals
    phi_c ~ normal(0, 5);          // cases NB2 dispersion
    phi_d ~ normal(0, 5);          // deaths NB2 dispersion

    // Baseline FOI priors (keep transmission trend relatively stable)
    alpha_lambda_s ~ normal(0, sigma_lambda_s); // province intercepts

    // Case detection Se hyperpriors 
    sigma_alpha_Se ~ exponential(0.5);
    sigma_delta_Se ~ exponential(2);

    // False positives FP hyperpriors
    sigma_alpha_FP ~ exponential(0.5);
    sigma_delta_FP ~ exponential(2);

    // CFR hyperpriors
    sigma_gamma_CFR ~ exponential(1);

    // Kappa (Lee–Carter) priors
    beta_kappa_a ~ normal(0, 1);        // baseline age RR
    raw_h_kappa_a ~ normal(0, 1);       // raw loadings (center+normalize in TP)
    sigma_gamma_kappa ~ exponential(1); // RW innovation SD
    sigma_delta_kappa ~ exponential(2); // SD for delta_kappa_ys

    // CFR priors
    beta_CFR_a ~ normal(0, 1);
    raw_h_CFR_a ~ normal(0, 1);

    // Clinical fraction priors (weakly informative to favor low proportion ~10%)
    nu_a ~ beta(1, 9);  // mean = 1/(1+9) = 0.1 (10%)

    // Prior for seasonal Fourier coefficients (region–period baseline)
    for(h in 1:H){
        for(p_ind in 1:P){
            for(r_ind in 1:R){
                a_cos[h, p_ind, r_ind] ~ normal(0, 0.5);
                a_sin[h, p_ind, r_ind] ~ normal(0, 0.5);
            }
        }
    }

    // Prior for seasonal Fourier deviation coefficients (year–province deviations)
    for(h in 1:H){
        for(y_ind in 1:Y){
            for(s_ind in 1:S){
                delta_cos[h, y_ind, s_ind] ~ normal(0, 0.2);
                delta_sin[h, y_ind, s_ind] ~ normal(0, 0.2);
            }
        }
    }

    // Yearly i.i.d. random effects for baseline FOI residual epsilon_lambda_ys 
    for(s_ind in 1:S){
        for(y_ind in 1:Y){
            epsilon_lambda_ys[y_ind, s_ind] ~ normal(0, sigma_epsilon_lambda);
        }
    }

    // Region–year random walk for alpha_Se_yr (case detection sensitivity)
    for (r_ind in 1:R) {
        alpha_Se_yr[1, r_ind] ~ normal(0, 0.5); // Prior for first year
        for (y_ind in 2:Y) {
            alpha_Se_yr[y_ind, r_ind] ~ normal(alpha_Se_yr[y_ind - 1, r_ind], sigma_alpha_Se); // Random walk prior
        }
    }
    
    // Province i.i.d. random effects for delta_Se_s (case detection sensitivity)
    delta_Se_s ~ normal(0, sigma_delta_Se);

    // Region–year random walk for alpha_FP_yr (false positives)
    for (r_ind in 1:R) {
        alpha_FP_yr[1, r_ind] ~ normal(0, 1); // Prior for first year
        for (y_ind in 2:Y) {
            alpha_FP_yr[y_ind, r_ind] ~ normal(alpha_FP_yr[y_ind - 1, r_ind], sigma_alpha_FP);
        }
    }

    // Province i.i.d. random effects for delta_FP_s (false positives)
    delta_FP_s ~ normal(0, sigma_delta_FP);

    // Region–period anchored random walk for gamma_kappa_pr (kappa)
    for (r_ind in 1:R){
        gamma_kappa_pr_base[r_ind] ~ normal(0, 0.5); // Prior for first period (half-normal via <lower=0>)
        for (p_ind in 1:(P-1)){
            gamma_kappa_pr_innov[p_ind, r_ind] ~ normal(0, sigma_gamma_kappa); // RW innovations for periods 2..P
        }
    }

    // Province–year i.i.d. random effects for delta_kappa_ys (kappa)
    for (s_ind in 1:S){
        for (y_ind in 1:Y){
            delta_kappa_ys[y_ind, s_ind] ~ normal(0, sigma_delta_kappa);
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
  
}


