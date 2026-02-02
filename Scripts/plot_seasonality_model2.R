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
  kappa_seasonal <- array(0.0, dim = c(n_draws, n_points, P, R))
  
  for (i in 1:n_draws) {
    for (p in 1:P) {
      #p <- period_id[y]  # period for this year
      for (r in 1:R) {
        #r <- region_id[s]  # region for this province
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
        kappa_seasonal[i, , p, r] <- n_points * exp(log_kappa) / sum_exp
      }
    }
  }
  
  # Compute quantiles and mean
  qfun <- function(x) quantile(x, probs = probs, names = FALSE, type = 7)
  
  # Build data frame with smooth curves (year-province-month)
  results <- expand.grid(
    m_grid = month_grid,
    p_idx = 1:P,
    r_idx = 1:R
  ) %>%
    mutate(
      grid_idx = rep(1:n_points, times = P * R),
      q025 = NA_real_,
      q500 = NA_real_,
      q975 = NA_real_,
      mean = NA_real_
    )
  
  for (i in 1:nrow(results)) {
    g_idx <- results$grid_idx[i]
    p <- results$p_idx[i]
    sr <- results$r_idx[i]
    
    vals <- kappa_seasonal[, g_idx, p, sr]  # Fixed: use sr instead of r
    qs <- qfun(vals)
    
    results$q025[i] <- qs[1]
    results$q500[i] <- qs[2]
    results$q975[i] <- qs[3]
    results$mean[i] <- mean(vals)
  }
  

  
  # Add period labels
  results <- results %>%
    mutate(period = factor(p_idx, levels = 1:P, labels = c("2003 to 2012", "2013 to 2024")),
           region = factor(r_idx, levels = 1:R, labels = paste0("Zone ", 1:R))) 
  
  results
}

plot_seasonality_by_province <- function(seasonality_df, facet_by_period = FALSE) {
  # Seasonality now varies by period-region
  # Facet by region (13 regions), color by period
  
  # Add region to the dataframe
  seasonality_with_region <- seasonality_df 
  
  # Average across provinces within each region and period
  seasonality_avg <- seasonality_with_region %>%
    group_by(m_grid, r_idx, region, period) %>%
    summarize(
      mean = mean(mean),
      .groups = "drop"
    )
  
  # Color palette: blue for period 1, red for period 2
  period_colors <- c("#2166ac", "#b2182b")  # Nice blue and red
  
  ggplot(seasonality_avg, aes(x = m_grid, y = mean, color = factor(period), group = period)) +
    geom_line(size = 1.2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
    facet_wrap(~ region, scales = "free_y", ncol = 4) +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +
    scale_color_manual(values = period_colors, name = "Period") +
    labs(x = "Month", y = "Seasonal modifier (kappa)",
         title = "Seasonal patterns by region and period",
         subtitle = "Gray dashed line shows baseline (kappa = 1). Averaged across provinces within each region.") +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_seasonality_heatmap <- function(seasonality_df) {
  # Create heatmap with months on x-axis, regions on y-axis
  # Separate panels for each period
  
  # Round month to integer for cleaner heatmap
  seasonality_monthly <- seasonality_df %>%
    mutate(month = round(m_grid)) %>%
    filter(month >= 1 & month <= 12) %>%
    group_by(month, r_idx, region, p_idx, period) %>%
    summarize(mean = mean(mean), .groups = "drop")
  
  # Create heatmap
  ggplot(seasonality_monthly, aes(x = month, y = region, fill = mean)) +
    geom_tile(color = "white", linewidth = 0.5) +
    facet_wrap(~ period, ncol = 2) +
    scale_fill_gradient2(
      low = "#2166ac", 
      mid = "white", 
      high = "#b2182b",
      midpoint = 1,
      name = "Seasonal\nmodifier",
      limits = c(NA, NA)
    ) +
    scale_x_continuous(breaks = 1:12, labels = month.abb, expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
      x = "Month",
      y = "Region",
      title = "Seasonal patterns by region and period",
      subtitle = "White indicates baseline (kappa = 1); blue = below average, red = above average"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "gray90", color = NA),
      strip.text = element_text(face = "bold", size = 12)
    )
}

# Example usage:
seasonality_df <- plot_seasonality_model2(fit, stan_data)

# Line plot: Facet by region, color by period
p2 <- plot_seasonality_by_province(seasonality_df, facet_by_period = TRUE)
print(p2)

# Heatmap: months x regions, facet by period
p3 <- plot_seasonality_heatmap(seasonality_df)
print(p3)

pname <- here("Plots", "seasonality_by_region_model3_line.png")
png(pname, width = 15, height = 10, units = "in", res = 300)
print(p2)
dev.off()

pname <- here("Plots", "seasonality_by_region_model3_heatmap.png")
png(pname, width = 15, height = 8, units = "in", res = 300)
print(p3)
dev.off()


# # Save
# ggsave(here("Figures", "seasonality_by_province_model3.png"), p1, 
#        width = 12, height = 8, dpi = 300)

