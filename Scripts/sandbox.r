CFR_mat <- as.matrix(fit, pars = "CFR_pa")
idx <- stringr::str_match(colnames(CFR_mat), "CFR_pa\\[(\\d+),(\\d+)\\]")

CFR_summ <- tibble(par = colnames(CFR_mat), p_idx = as.integer(idx[, 2]), a_idx = as.integer(idx[, 3])) |>
 left_join(
   as_tibble(CFR_mat) |> mutate(.draw = row_number()) |>
     pivot_longer(cols = -.draw, names_to = "par", values_to = "CFR"),
   by = "par"
 ) |>
 group_by(p_idx, a_idx) |>
 summarize(
   mean = mean(CFR), q025 = quantile(CFR, 0.025), q975 = quantile(CFR, 0.975),
   .groups = "drop"
 ) |>
 mutate(
   period = paste0("Period ", p_idx),
   age_group = factor(age_levels[a_idx], levels = age_levels)
 )

CFR_summ |>
  mutate(
      age_group = factor(age_group, levels = age_levels),
      period = factor(period, levels = paste0("Period ", 1:4), 
                     labels = c("2003 to 2007", "2008 to 2012", "2013 to 2018", "2019 to 2024"))
    ) |>
  ggplot() +
  geom_abline(intercept = 0, slope = 0, linetype = "21", color = "gray30", linewidth = 1, alpha = 0.7) +
  geom_point(aes(x = age_group, y = mean * 100, color = period), size = 3, alpha = 0.8,
             position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(x = age_group, ymin = q025 * 100, ymax = q975 * 100, color = period), alpha = 0.8,
                width = 0, position = position_dodge(width = 0.5), linewidth = 1.1) +
  scale_color_viridis_d(option = "magma", end = 0.8, name = "") +
  labs(x = "", y = "Case-Fatality Rate (%)", 
       title = "Case-Fatality Rate") +
  theme_bw(base_size = 14) +
  coord_flip()




cfr_plot <- function(fit, age_levels, period_labels = c("2003 to 2007", "2008 to 2012", "2013 to 2018", "2019 to 2024")) {
  CFR_mat <- as.matrix(fit, pars = "CFR_pa")
  idx <- stringr::str_match(colnames(CFR_mat), "CFR_pa\\[(\\d+),(\\d+)\\]")
  
  CFR_summ <- tibble(par = colnames(CFR_mat), p_idx = as.integer(idx[, 2]), a_idx = as.integer(idx[, 3])) |>
    left_join(
      as_tibble(CFR_mat) |> mutate(.draw = row_number()) |>
        pivot_longer(cols = -.draw, names_to = "par", values_to = "CFR"),
      by = "par"
    ) |>
    group_by(p_idx, a_idx) |>
    summarize(
      mean = mean(CFR), q025 = quantile(CFR, 0.025), q975 = quantile(CFR, 0.975),
      .groups = "drop"
    ) |>
    mutate(
      period = paste0("Period ", p_idx),
      age_group = factor(age_levels[a_idx], levels = age_levels)
    )
  
  CFR_summ |>
    mutate(
      age_group = factor(age_group, levels = age_levels),
      period = factor(period, levels = paste0("Period ", 1:4), labels = period_labels)
    ) |>
    ggplot() +
    geom_abline(intercept = 0, slope = 0, linetype = "21", color = "gray30", linewidth = 1, alpha = 0.7) +
    geom_point(aes(x = age_group, y = mean * 100, color = period), size = 3, alpha = 0.8,
               position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(x = age_group, ymin = q025 * 100, ymax = q975 * 100, color = period), alpha = 0.8,
                  width = 0, position = position_dodge(width = 0.5), linewidth = 1.1) +
    scale_color_viridis_d(option = "magma", end = 0.8, name = "") +
    labs(x = "", y = "Case-Fatality Rate (%)", title = "Case-Fatality Rate") +
    theme_bw(base_size = 14) +
    coord_flip()
}

cfr_plot(fit, age_levels)
