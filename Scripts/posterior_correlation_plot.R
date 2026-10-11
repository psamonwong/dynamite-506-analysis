library(rstan)
library(tidyverse)

# 1. Extract posterior samples (using your confirmed dimensions)
post <- rstan::extract(fit, pars = c("Se_y", "lambda_base_ys"))
Se_y_raw <- post$Se_y             # Dims: [draws, 13 (regions), 22 (years)]
lambda_raw <- post$lambda_base_ys # Dims: [draws, 22 (years), 77 (provinces)]

n_draws <- dim(Se_y_raw)[1]
n_prov <- 77
n_regions <- 13
n_years <- 22                     
region_id <- stan_data$region_id

# 2. Build the tidy data frame for all draws, provinces, and years
set.seed(42)
draw_sample <- sample(n_draws, min(200, n_draws))

df_list <- list()
for (y in 1:n_years) {
  for (p in 1:n_prov) {
    reg <- region_id[p]
    
    df_list[[length(df_list) + 1]] <- data.frame(
      draw = 1:length(draw_sample),
      province = factor(p),
      region = factor(reg),
      year = y + 2002,
      Se_y = Se_y_raw[draw_sample, reg, y],          
      lambda_base_ys = lambda_raw[draw_sample, y, p] 
    )
  }
}
df_all <- bind_rows(df_list)

# 3. Plot: 1 Figure per Region (Faceted by Year)
for (r in 1:n_regions) {
  p_region <- df_all %>%
    filter(region == r) %>%
    ggplot(aes(x = log10(lambda_base_ys), y = Se_y, color = province)) +
    geom_point(alpha = 0.4, size = 1) +
    facet_wrap(~ year, scales = 'free_y') +
    labs(
      title = paste("Posterior Correlation — Region", r),
      subtitle = "Faceted by Year (Points represent different provinces in this region)",
      x = "Province log10 baseline FOI",
      y = "Region reporting proportions"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none") # Hides legend since multiple provinces share the region
  
  print(p_region)
  
  # Optional: Save each region's plot to disk
   ggsave(filename = paste0("Plots/Posterior_correlation/correlation_region_", r, ".png"), plot = p_region, width = 12, height = 10)
}
