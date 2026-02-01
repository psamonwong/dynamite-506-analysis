# Plot seasonal patterns (kappa_seasonal_yms) from Model with period-region Fourier
library(tidyverse)
library(here)

plot_seasonality_model2 <- function(fit, stan_data, probs = c(0.025, 0.5, 0.975), n_points = 100) {
  # Extract dimensions
  Y <- stan_data$Y
  M <- stan_data$M
  P <- stan_data$P
  R <- stan_data$R
  S <- stan_data$S
  H <- stan_data$H
  
  # Extract mappings
  period_id <- stan_data$period_id   # [Y] year to period
  region_id <- stan_data$region_id   # [S] province to region
  
  # Extract Fourier coefficients (now period-region indexed)
  a_cos_draws <- rstan::extract(fit, pars = "a_cos", permuted = TRUE)$a_cos  # [n_draws, H, P, R]
  a_sin_draws <- rstan::extract(fit, pars = "a_sin", permuted = TRUE)$a_sin  # [n_draws, H, P, R]
  n_draws <- dim(a_cos_draws)[1]
  
  # Create smooth monthly grid (continuous over the year)
  month_grid <- seq(1, M + 1, length.out = n_points)  # Wrap around to show continuity
  
  # Reconstruct kappa_seasonal at fine resolution from Fourier coefficients
  # Now we need year-province specific, using period_id and region_id
  kappa_seasonal <- array(0.0, dim = c(n_draws, n_points, Y, S))
  
  for (i in 1:n_draws) {
    for (y in 1:Y) {
      p <- period_id[y]  # period for this year
      for (s in 1:S) {
        r <- region_id[s]  # region for this province
        log_kappa <- numeric(n_points)
        for (m_idx in 1:n_points) {
          m <- month_grid[m_idx]
          theta <- 2 * pi * m / M
          f <- 0
          for (h in 1:H) {
            f <- f + a_cos_draws[i, h, p, r] * cos(h * theta) + 
                     a_sin_draws[i, h, p, r] * sin(h * theta)
          }
          log_kappa[m_idx] <- f
        }
        # Normalize to mean 1
        sum_exp <- sum(exp(log_kappa))
        kappa_seasonal[i, , y, s] <- n_points * exp(log_kappa) / sum_exp
      }
    }
  }
  
  # Compute quantiles and mean
  qfun <- function(x) quantile(x, probs = probs, names = FALSE, type = 7)
  
  # Build data frame with smooth curves (year-province-month)
  results <- expand.grid(
    m_grid = month_grid,
    y_idx = 1:Y,
    s_idx = 1:S
  ) %>%
    mutate(
      grid_idx = rep(1:n_points, times = Y * S),
      q025 = NA_real_,
      q500 = NA_real_,
      q975 = NA_real_,
      mean = NA_real_
    )
  
  for (i in 1:nrow(results)) {
    g_idx <- results$grid_idx[i]
    y <- results$y_idx[i]
    s <- results$s_idx[i]
    
    vals <- kappa_seasonal[, g_idx, y, s]
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
  
  # Add period labels
  results <- results %>%
    mutate(period = period_id[y_idx])
  
  results
}

plot_seasonality_by_province <- function(seasonality_df, facet_by_period = FALSE) {
  # Seasonality now varies by period-region
  # Can facet by province or by province+period
  if (facet_by_period) {
    ggplot(seasonality_df, aes(x = m_grid, y = mean, color = factor(period), group = period)) +
      geom_line(size = 1.2) +
      geom_ribbon(aes(ymin = q025, ymax = q975, fill = factor(period)), alpha = 0.2) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      facet_wrap(~ province, scales = "free_y") +
      scale_x_continuous(breaks = 1:12, labels = month.abb) +
      labs(x = "Month", y = "Seasonal modifier (kappa)",
           title = "Seasonal patterns by province and period",
           subtitle = "Red dashed line shows baseline (kappa = 1). Smooth Fourier curves with 95% credible bands.",
           color = "Period", fill = "Period") +
      theme_bw(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  } else {
    # Average across years within each province
    seasonality_avg <- seasonality_df %>%
      group_by(m_grid, s_idx, province) %>%
      summarize(
        q025 = mean(q025),
        mean = mean(mean),
        q975 = mean(q975),
        .groups = "drop"
      )
    
    ggplot(seasonality_avg, aes(x = m_grid, y = mean, group = 1)) +
      geom_line(size = 1.2, color = "#1f78b4") +
      geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.3, fill = "#1f78b4") +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      facet_wrap(~ province, scales = "free_y") +
      scale_x_continuous(breaks = 1:12, labels = month.abb) +
      labs(x = "Month", y = "Seasonal modifier (kappa)",
           title = "Seasonal patterns by province (averaged across years)",
           subtitle = "Red dashed line shows baseline (kappa = 1). Smooth Fourier curves with 95% credible bands.") +
      theme_bw(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
}

plot_seasonality_overall <- function(seasonality_df) {
  # Average across all provinces and years (period-region variation)
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
         title = "Overall seasonal pattern (averaged across provinces and years)",
         subtitle = "Red dashed line shows baseline (kappa = 1). Smooth Fourier curve.") +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_seasonality_by_period <- function(seasonality_df) {
  # Average within each period across all provinces
  seasonality_period <- seasonality_df %>%
    group_by(m_grid, period) %>%
    summarize(
      q025 = mean(q025),
      mean = mean(mean),
      q975 = mean(q975),
      .groups = "drop"
    )
  
  ggplot(seasonality_period, aes(x = m_grid, y = mean, color = factor(period), fill = factor(period), group = period)) +
    geom_line(size = 1.2) +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +
    labs(x = "Month", y = "Seasonal modifier (kappa)",
         title = "Seasonal patterns by period (averaged across provinces)",
         subtitle = "Red dashed line shows baseline (kappa = 1). Smooth Fourier curves.",
         color = "Period", fill = "Period") +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Example usage:
# seasonality_df <- plot_seasonality_model2(fit, stan_data)

# Facet by province (averaged across years)
# p1 <- plot_seasonality_by_province(seasonality_df)
# print(p1)

# Facet by province, colored by period
# p2 <- plot_seasonality_by_province(seasonality_df, facet_by_period = TRUE)
# print(p2)

# Overall average
# p3 <- plot_seasonality_overall(seasonality_df)
# print(p3)

# By period (averaged across provinces)
# p4 <- plot_seasonality_by_period(seasonality_df)
# print(p4)

# # Save
# ggsave(here("Figures", "seasonality_by_province_model3.png"), p1, 
#        width = 12, height = 8, dpi = 300)

