library(here)
library(tidyverse)
library(rstan)
load(here("Rout", "stan_data_01.RData"))

stan_data$S_t <- 1
stan_data$N_tas <- stan_data$N_tas  %>% apply(c(1,3), sum)

system(paste('Rscript --vanilla Scripts/run_cmd_simple.R'))


mod <- stan_model(here("Stan_models", "Model2_simple.stan"), verbose = FALSE)

load(here("Rout", "fit_02.RData"))

rstan::traceplot(fit, pars = c("phi_case", "Se_y[1]", "Sp_y[1]", "Se_y[2]", "Sp_y[2]",
                              "log_lambda_0", "sigma_lambda_s", "sigma_lambda_y"))

xx













fits <- replicate(20,
 optimizing(mod, data = stan_data, iter = 10000, verbose = TRUE, 
                      hessian = FALSE, as_vector = FALSE),
  simplify = FALSE
)

# 1. Check objective function values - should be similar if well-identified
values <- sapply(fits, function(f) f$value)
hist(values, breaks = 20, main = "Optimization objective values", 
     xlab = "Negative log-posterior")
# Large spread or multiple clusters → multimodality/poor identification

# 2. Compare specific parameters across runs
# Extract parameter estimates from each run
extract_param <- function(fits, param_name) {
  sapply(fits, function(f) {
    p <- f$par[[param_name]]
    if(is.null(p)) return(NA)
    if(length(p) > 1) return(mean(p))  # average if vector/array
    return(p)
  })
}

# Key parameters to check
params_to_check <- c("phi_case", "log_lambda_0", "sigma_lambda_s", "sigma_lambda_y")
# Create diagnostic plots
par(mfrow = c(2, 1))
for(p in params_to_check) {
  vals <- extract_param(fits, p)
  hist(vals, main = p, xlab = "Value", breaks = 10)
  abline(v = median(vals, na.rm = TRUE), col = "red", lwd = 2)
}
par(mfrow = c(1, 1))


# 3. Coefficient of variation - high CV indicates poor identification
cv_summary <- data.frame(
  parameter = params_to_check,
  mean = sapply(params_to_check, function(p) mean(extract_param(fits, p), na.rm = TRUE)),
  median = sapply(params_to_check, function(p) median(extract_param(fits, p), na.rm = TRUE)),
  sd = sapply(params_to_check, function(p) sd(extract_param(fits, p), na.rm = TRUE))
) %>%
  mutate(cv = sd / abs(mean)) %>%
  arrange(desc(cv))

print(cv_summary)

