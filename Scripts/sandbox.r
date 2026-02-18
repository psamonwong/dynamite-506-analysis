
# Extract and plot Se_y over time
plot_Se_over_time <- function(probs = c(0.025, 0.5, 0.975)) {
  stopifnot(exists("fit"), exists("stan_data"), exists("calendar"))
  Y <- stan_data$Y
  
  # Extract Se_y draws [Y]
  Se_mat <- as.matrix(fit, pars = "Se_y")
  if (is.null(colnames(Se_mat))) stop("Se_y columns have no names; re-run sampling with pars='Se_y'.")
  
  # Parse indices from column names [y]
  idx <- stringr::str_match(colnames(Se_mat), "Se_y\\[(\\d+),(\\d+)\\]")
  if (any(is.na(idx))) stop("Failed to parse Se_y indices from column names.")
  y_idx_col <- as.integer(idx[, 3])
  r_idx_col <- as.integer(idx[, 2])
 
  # Create parameter map
  param_map <- tibble::tibble(
    par = colnames(Se_mat), 
    y_idx = y_idx_col,
    r_idx = r_idx_col
  )
  
  # Convert to long format
  Se_long <- tibble::as_tibble(Se_mat) |>
    dplyr::mutate(.draw = dplyr::row_number()) |>
    tidyr::pivot_longer(cols = -.draw, names_to = "par", values_to = "Se") |>
    dplyr::left_join(param_map, by = "par")
  
  # Summarize by year
  Se_summ <- Se_long |>
    dplyr::group_by(y_idx, r_idx) |>
    dplyr::summarize(
      mean = mean(Se),
      q025 = stats::quantile(Se, probs[1], names = FALSE),
      q500 = stats::quantile(Se, probs[2], names = FALSE),
      q975 = stats::quantile(Se, probs[3], names = FALSE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      year = min(calendar$year) + y_idx - 1,
      Zone = sprintf("Zone%02d", r_idx),
      region = factor(
        paste0("Health Region ", r_idx),
        levels = paste0("Health Region ", sort(unique(r_idx)))
      )
    )

# Plot
se_trend <- ggplot(Se_summ, aes(x = year)) +
  annotate(
    "rect",
    xmin = 2020, xmax = 2022,
    ymin = 0, ymax = 1,
    fill = "#CB80AB", alpha = 0.6
  ) +
  geom_ribbon(aes(ymin = q025, ymax = q975), fill = "#295F98", alpha = 0.3) +
  geom_line(aes(y = q500), color = "#295F98", linewidth = 1) +
  geom_point(aes(y = q500), color = "#295F98", size = 1.2) +
  labs(
    x = "Year",
    y = "Case Detection Sensitivity",
    title = "Temporal Evolution of Case Detection Sensitivity",
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_blank()) +
  facet_wrap(.~region, ncol = 3) 
  
   map_se <- map_data %>%
    left_join(
      Se_summ %>% filter(year == 2024),
      by = "Zone"
    ) %>%
    ggplot() +
    geom_sf(aes(fill = q500), color = "white", linewidth = 0.1) +
    scale_fill_gradientn(
      name = "",
      colours = c("#f7fbff", "#6baed6", "#08306b"),
      limits = c(0, 1)
    ) +   
     coord_sf(xlim = c(96, 106), ylim = c(5, 21), expand = FALSE) +
    theme_void() +
    theme(axis.text = element_blank(),
          legend.position = "right") +
    labs(caption = "Estimated Case Detection Sensitivity\nby Health Region (2024)") +
    theme(plot.caption = element_text(hjust = 0.5, vjust = 20, size = 12, face = "bold"))  

  combined_se_plot <- se_trend + map_se + plot_layout(ncol = 2, widths  = c(3, 1)) 
  combined_se_plot + plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))
  
}

plot_Se_over_time()
map_data




