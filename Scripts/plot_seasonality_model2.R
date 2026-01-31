# Plot seasonal patterns (kappa_seasonal_yms) from Model 2
library(tidyverse)
library(here)

plot_seasonality_model2 <- function(fit, stan_data, probs = c(0.025, 0.5, 0.975), n_points = 100) {
  # Extract dimensions
  Y <- stan_data$Y
  M <- stan_data$M
  S <- stan_data$S
  H <- stan_data$H
  
  # Extract Fourier coefficients
  a_cos_draws <- rstan::extract(fit, pars = "a_cos", permuted = TRUE)$a_cos  # [n_draws, H, S]
  a_sin_draws <- rstan::extract(fit, pars = "a_sin", permuted = TRUE)$a_sin  # [n_draws, H, S]
  n_draws <- dim(a_cos_draws)[1]
  
  # Create smooth monthly grid (continuous over the year)
  month_grid <- seq(1, M + 1, length.out = n_points)  # Wrap around to show continuity
  
  # Reconstruct kappa_seasonal at fine resolution from Fourier coefficients
  kappa_seasonal <- array(0.0, dim = c(n_draws, n_points, S))
  
  for (i in 1:n_draws) {
    for (s in 1:S) {
      log_kappa <- numeric(n_points)
      for (m_idx in 1:n_points) {
        m <- month_grid[m_idx]
        theta <- 2 * pi * m / M
        f <- 0
        for (h in 1:H) {
          f <- f + a_cos_draws[i, h, s] * cos(h * theta) + 
                   a_sin_draws[i, h, s] * sin(h * theta)
        }
        log_kappa[m_idx] <- f
      }
      # Normalize to mean 1
      sum_exp <- sum(exp(log_kappa))
      kappa_seasonal[i, , s] <- n_points * exp(log_kappa) / sum_exp
    }
  }
  
  # Compute quantiles and mean
  qfun <- function(x) quantile(x, probs = probs, names = FALSE, type = 7)
  
  # Build data frame with smooth curves
  results <- expand.grid(
    m_grid = month_grid,
    s_idx = 1:S
  ) %>%
    mutate(
      m_idx = floor(m_grid),
      grid_idx = rep(1:n_points, times = S),
      q025 = NA_real_,
      q500 = NA_real_,
      q975 = NA_real_,
      mean = NA_real_
    )
  
  for (i in 1:nrow(results)) {
    g_idx <- results$grid_idx[i]
    s <- results$s_idx[i]
    
    vals <- kappa_seasonal[, g_idx, s]
    qs <- qfun(vals)
    
    results$q025[i] <- qs[1]
    results$q500[i] <- qs[2]
    results$q975[i] <- qs[3]
    results$mean[i] <- mean(vals)
  }
  
  # Add province labels if available
  if (exists("site_levels", inherits = TRUE)) {
    results <- results %>%
      mutate(province = factor(s_idx, levels = 1:S, labels = site_levels))
  } else {
    results <- results %>%
      mutate(province = factor(s_idx, levels = 1:S, labels = paste0("Province ", 1:S)))
  }
  
  results
}

plot_seasonality_by_province <- function(seasonality_df) {
  # In Model2, seasonality is constant across years (no temporal variation)
  # Use smooth geom_line since data is now at fine resolution
  ggplot(seasonality_df, aes(x = m_grid, y = mean, group = 1)) +
    geom_line(size = 1.2, color = "#1f78b4") +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.3, fill = "#1f78b4") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    facet_wrap(~ province, scales = "free_y") +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +
    labs(x = "Month", y = "Seasonal modifier (kappa)",
         title = "Seasonal patterns by province (Model 2)",
         subtitle = "Red dashed line shows baseline (kappa = 1). Smooth Fourier curves with 95% credible bands.") +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_seasonality_overall <- function(seasonality_df) {
  # Average across all provinces (no year variation in Model2)
  seasonality_avg <- seasonality_df %>%
    group_by(m_grid) %>%
    summarize(
      q025 = mean(q025),
      mean = mean(mean),
      q975 = mean(q975),
      .groups = "drop"
    )
  
  ggplot(seasonality_avg, aes(x = m_grid, y = mean, group = 1)) +
    geom_line(size = 1.5, color = "#1f78b4") +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.3, fill = "#1f78b4") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +
    labs(x = "Month", y = "Seasonal modifier (kappa)",
         title = "Overall seasonal pattern (averaged across provinces)",
         subtitle = "Red dashed line shows baseline (kappa = 1). Smooth Fourier curve.") +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Example usage:
seasonality_df <- plot_seasonality_model2(fit, stan_data)

# Facet by province
p1 <- plot_seasonality_by_province(seasonality_df)
print(p1)

# # Save
# ggsave(here("Figures", "seasonality_by_province_model2.png"), p1, 
#        width = 12, height = 8, dpi = 300)

