// Model 1: Fouriers for seasonality in infection rate
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

parameters{
    // Overdispersion parameters for cases and deaths negative binomial likelihoods
    real<lower=0> phi_c;           // overdispersion parameter for reported cases
    real<lower=0> phi_d;           // overdispersion parameter for reported deaths

     // Baseline infection rate parameters 
    array[S] real alpha_lambda_s;          // alpha_lambda_s[s]: baseline log infection rate for province s
    real<lower=-1, upper=1> rho_lambda;    // AR(1) correlation for yearly deviations in log infection rate
    real<lower=0> eta_lambda;              // AR(1) innovation std for yearly deviations in log infection rate
    array[Y, S] real epsilon_lambda_ys;    // epsilon_lambda_ys[y, s]: yearly deviation in log infection rate for year y, province s
    // Seasonal Fourier coefficients
    array[H, P, R] real a_cos;        // a_cos[h, p, r]: cosine coefficient for harmonic h, period p, region r
    array[H, P, R] real a_sin;        // a_sin[h, p, r]: sine coefficient for harmonic h, period p, region r
    array[H, Y, S] real delta_cos;    // delta_cos[h, y, s]: deviation cosine coefficient for harmonic h, year y, province s
    array[H, Y, S] real delta_sin;    // delta_sin[h, y, s]: deviation sine coefficient for harmonic h, year y, province s
    // Age modifier parameters (lee–carter structure)
    array[A] real beta_kappa_a;                  // beta_kappa_a[a]: age effect for age a
    array[P, R] real gamma_kappa_pr;             // gamma_kappa_pr[p, r]: period-region interaction effect for period p, region r with age sensitivity
    array[A] real h_kappa;                       // h_kappa[a]: age sensitivity for age a 
    array[Y, S] real delta_kappa_ys;            // delta_kappa_ys[y, s]: yearly deviation for year y, province s
    array[Y, A, S] real u_kappa_yas;            // u_kappa_yas[y, a, s]: random effect for year y, age a, province s
 

}

transformed parameters {
    // Latent monthly totals
    array[T, S] real<lower=0> mu_case_ts;        // mu_case_t,s mean reported cases at time t, province s
    array[T, S] real<lower=0> mu_death_ts;       // mu_death_t,s mean reported deaths at time t, province s
    
    // Infection rate and related quantities
    array[T, A, S] real<lower=0> lambda_tas;    // lambda_t,a,s infection rate at time t, age a, province s
    array[Y, M, S] real<lower=0> kappa_month_yms;   // kappa_month_y,m,s seasonal modifier for month m, year y, province s
    array[Y, A, S] real<lower=0> kappa_age_yas;     // kappa_age_y,a,s age modifier for age a, year y, province s

    // Transmission rate parameters
    array[Y, S] real log_lambda_base;        // log baseline infection rate for year y, province s


    // Evolution of baseline FOI, log_lambda_base, with AR(1) process could be added here
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
            // Centering to have mean 1 over months
            real sum_exp_m = 0; 
            for(m_ind in 1:M) sum_exp_m += exp(log_kappa_yms[m_ind]); // sum over months 
            for(m_ind in 1:M) kappa_month_yms[y_ind, m_ind, s_ind] = M * exp(log_kappa_yms[m_ind])/ sum_exp_m; // normalize to mean 1
        }
    }

    // Age modifier: kappa_age_yas Lee–Carter on log scale, normalized to mean 1 over age
    for(y_ind in 1:Y){
        for(s_ind in 1:S){
            vector[A] log_kappa_age; // log age modifier for age a, year y, province s
            for(a_ind in 1:A){
                log_kappa_age[a_ind] =
                    beta_kappa_a[a_ind] // age effect 
                  + gamma_kappa_pr[period_id[y_ind], region_id[s_ind]] * h_kappa[a_ind] // period-region interaction with age sensitivity
                  + delta_kappa_ys[y_ind, s_ind] // yearly deviation for year y, province s
                  + u_kappa_yas[y_ind, a_ind, s_ind]; // random effect for year y, age a, province s
            }
            // Centering to have mean 1 over age
            real sum_exp_a = 0;
            for(a_ind in 1:A) sum_exp_a += exp(log_kappa_age[a_ind]); // sum over ages
            for(a_ind in 1:A) kappa_age_yas[y_ind, a_ind, s_ind] = A * exp(log_kappa_age[a_ind]) / sum_exp_a; // normalize to mean 1
        }
    }




     // Calculating I_tas, lambda_tas, C_tas and D_tas
     for(s_ind in 1:S){
        for(t_ind in 1:T){


            for(a_ind in 1:A){
                // Force of infection at time t, age a, province s 
                lambda_tas[t_ind, a_ind, s_ind] = exp(log_lambda_base[year_id[t_ind], s_ind]) * // baseline infection rate
                                                    kappa_month_yms[year_id[t_ind], month_id[t_ind], s_ind] * // seasonal modifier
                                                    kappa_age_yas[year_id[t_ind], a_ind, s_ind]; // age modifier







            }





        }
     }





}

model {
    // Priors 
    phi_c ~ normal(0, 5);          // Prior for overdispersion in reported cases
    phi_d ~ normal(0, 5);          // Prior for overdispersion in reported deaths

    // Priors for baseline infection rate parameters
    alpha_lambda_s ~ normal(0, 1); // Prior for baseline log infection rate by province
    rho_lambda ~ uniform(-1, 1);   // Prior for AR(1) correlation
    eta_lambda ~ normal(0, 1);     // Prior for AR(1) innovation std

    // Prior for seasonal Fourier coefficients
    for(h in 1:H){
        for(p_ind in 1:P){
            for(r_ind in 1:R){
                a_cos[h, p_ind, r_ind] ~ normal(0, 0.5);
                a_sin[h, p_ind, r_ind] ~ normal(0, 0.5);
            }
        }
    }

    // Prior for seasonal Fourier deviation coefficients
    for(h in 1:H){
        for(y_ind in 1:Y){
            for(s_ind in 1:S){
                delta_cos[h, y_ind, s_ind] ~ normal(0, 0.2);
                delta_sin[h, y_ind, s_ind] ~ normal(0, 0.2);
            }
        }
    }


    // AR(1) process for yearly deviations in log infection rate by province epsilon_lambda_ys
    for(s_ind in 1:S){
        epsilon_lambda_ys[1, s_ind] ~ normal(0, 0.5); // Prior for first year deviation 
        for(y_ind in 2:Y){
            epsilon_lambda_ys[y_ind, s_ind] ~ normal(rho_lambda * epsilon_lambda_ys[y_ind - 1, s_ind], eta_lambda);
        }
    }






    // Likelihood: Reported cases and deaths by month and province R_case_ts, R_death_ts
    for(s_ind in 1:S){
        for(t_ind in 1:T){
        R_case_ts[t_ind, s_ind] ~ neg_binomial_2(mu_case_ts[t_ind, s_ind], phi_c);
        R_death_ts[t_ind, s_ind] ~ neg_binomial_2(mu_death_ts[t_ind, s_ind], phi_d);
        }
    }

  
}

generated quantities {
  
}


