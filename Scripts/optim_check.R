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
params_to_check <- c("phi_c", "phi_d",
                      "sigma_lambda_r", "sigma_lambda_s", "alpha_CFR_a_cv",
                      "sigma_delta_Se", "sigma_delta_Se_y", 
                      "sigma_delta_FP", "sigma_delta_FP_y",
                      "sigma_gamma_kappa", "sigma_gamma_CFR",
                      "sigma_seasonality_pr")

# Create diagnostic plots
par(mfrow = c(1, 1))
for(p in params_to_check) {
  vals <- extract_param(fits, p)
  hist(vals, main = p, xlab = "Value", breaks = 10)
  abline(v = median(vals, na.rm = TRUE), col = "red", lwd = 2)
}

# 3. Coefficient of variation - high CV indicates poor identification
cv_summary <- data.frame(
  parameter = params_to_check,
  mean = sapply(params_to_check, function(p) mean(extract_param(fits, p), na.rm = TRUE)),
  sd = sapply(params_to_check, function(p) sd(extract_param(fits, p), na.rm = TRUE))
) %>%
  mutate(cv = sd / abs(mean)) %>%
  arrange(desc(cv))

print(cv_summary)
# CV > 0.5 suggests weak identification; CV > 1.0 severe problems

# 4. Pairwise correlations between hyperparameters
# Extract all scalar hyperparameters
hyperpars <- c("sigma_lambda_r", "sigma_lambda_s", "alpha_CFR_a_cv",
               "sigma_delta_Se", "sigma_delta_Se_y", 
               "sigma_delta_FP", "sigma_delta_FP_y",
               "sigma_gamma_kappa", "sigma_gamma_CFR",
               "sigma_seasonality_pr")

hyper_matrix <- sapply(hyperpars, function(p) extract_param(fits, p))
cor(hyper_matrix, use = "complete.obs")
# High correlations (|r| > 0.7) suggest confounding

# 5. Check best vs worst fits
best_idx <- which.min(values)
worst_idx <- which.max(values)

compare_runs <- data.frame(
  parameter = hyperpars,
  best = sapply(hyperpars, function(p) extract_param(list(fits[[best_idx]]), p)),
  worst = sapply(hyperpars, function(p) extract_param(list(fits[[worst_idx]]), p)),
  median_all = sapply(hyperpars, function(p) median(extract_param(fits, p), na.rm = TRUE))
) %>%
  mutate(ratio = best / worst)

print(compare_runs)
# Large ratios indicate parameter is driving fit quality

# 6. For arrays/vectors - check element-wise variation
check_array_param <- function(fits, param_name) {
  # Get dimensions from first fit
  first_val <- fits[[1]]$par[[param_name]]
  if(is.null(first_val) || length(first_val) == 1) return(NULL)
  
  # Extract from all fits
  all_vals <- lapply(fits, function(f) as.vector(f$par[[param_name]]))
  mat <- do.call(rbind, all_vals)
  
  # Compute CV for each element
  cv_by_element <- apply(mat, 2, function(x) sd(x) / abs(mean(x)))
  
  return(data.frame(
    element = 1:length(cv_by_element),
    cv = cv_by_element
  ))
}

# Check FOI regional baselines
alpha_lambda_r_cv <- check_array_param(fits, "alpha_lambda_r")
if(!is.null(alpha_lambda_r_cv)) {
  plot(alpha_lambda_r_cv$element, alpha_lambda_r_cv$cv, 
       type = "h", lwd = 3, main = "alpha_lambda_r: CV by region",
       xlab = "Region", ylab = "Coefficient of Variation")
  abline(h = 0.5, col = "red", lty = 2)
}

# Check alpha_kappa_a age effects
alpha_kappa_a_cv <- check_array_param(fits, "alpha_kappa_a")
if(!is.null(alpha_kappa_a_cv)) {
  plot(alpha_kappa_a_cv$element, alpha_kappa_a_cv$cv, 
       type = "h", lwd = 3, main = "alpha_kappa_a: CV by age group",
       xlab = "Age group", ylab = "Coefficient of Variation")
  abline(h = 0.5, col = "red", lty = 2)
}

# Check alpha_CFR_a age effects
alpha_CFR_a_cv <- check_array_param(fits, "alpha_CFR_a")
if(!is.null(alpha_CFR_a_cv)) {
  plot(alpha_CFR_a_cv$element, alpha_CFR_a_cv$cv, 
       type = "h", lwd = 3, main = "alpha_CFR_a: CV by age group",
       xlab = "Age group", ylab = "Coefficient of Variation")
  abline(h = 0.5, col = "red", lty = 2)
}

# Print first few fits for inspection
cat("\n=== alpha_kappa_a from first 6 fits ===\n")
for(i in 1:min(20, length(fits))) {
  cat(sprintf("Fit %d: %s\n", i, 
              paste(round(fits[[i]]$par$alpha_kappa_a, 3), collapse = ", ")))
}
