
cfr_plot <- function(fit, age_levels, period_labels = c("2003 to 2007", "2008 to 2012", "2013 to 2018", "2019 to 2024"), ref_age = "25 to 35") {
  CFR_mat <- as.matrix(fit, pars = "CFR_pa")
  idx <- stringr::str_match(colnames(CFR_mat), "CFR_pa\\[(\\d+),(\\d+)\\]")
  
  ref_age_idx <- which(age_levels == ref_age)

  CFR_summ <-    tibble(par = colnames(CFR_mat), p_idx = as.integer(idx[, 2]), a_idx = as.integer(idx[, 3])) |>
      left_join(
        as_tibble(CFR_mat) |> mutate(.draw = row_number()) |>
          pivot_longer(cols = -.draw, names_to = "par", values_to = "CFR"),
        by = "par"
      ) |>
      group_by(.draw, p_idx) |>
      mutate(ref_CFR = CFR[a_idx == ref_age_idx],
             rel_CFR = CFR/ref_CFR) |>
      ungroup() |>
      group_by(p_idx, a_idx) |>
      summarize(
        mean_rel_CFR = mean(rel_CFR), q025_ref_CFR = quantile(rel_CFR, 0.025), q975_rel_CFR = quantile(rel_CFR, 0.975),
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
    filter(age_group != ref_age) |>
    ggplot() +
    geom_abline(intercept = 1, slope = 0, linetype = "21", color = "gray30", linewidth = 1, alpha = 0.7) +
    geom_point(aes(x = age_group, y =mean_rel_CFR, color = period), size = 3, alpha = 0.8,
               position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(x = age_group, ymin =q025_ref_CFR, ymax = q975_rel_CFR, color = period), alpha = 0.8,
                  width = 0, position = position_dodge(width = 0.5), linewidth = 1.1) +
    scale_color_viridis_d(option = "magma", end = 0.8, name = "") +
    labs(x = "", y = "Relative case-fatality rate", title = "Relative case-Fatality Rate\n(relative to 25 to 35 years)") +
    theme_bw(base_size = 14) +
    scale_y_continuous(breaks = seq(0, 10, by = 1)) +
    coord_flip()
}
