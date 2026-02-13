# Function to export traceplots for key parameters
export_traceplots <- function(dataset) {
  # Create output directory
  plot_dir <- here("Plots", "traceplots", dataset)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Key parameters to examine
  params_to_plot <- list(
    "dispersion" = c("phi_case", "phi_death"),
    "detection" = c("Se_y", "sigma_Se"),
    "foi_global" = c("mu_log_lambda_0"),
    "age_susceptibility" = c("rho_pa", "sigma_rho"),
    "cfr" = c("mu_logit_CFR", "sigma_CFR", "CFR_pa")
  )
  
  # Export traceplots
  for (group_name in names(params_to_plot)) {
    pars <- params_to_plot[[group_name]]
    
    # Check which parameters exist in the fit
    pars_exist <- pars[pars %in% names(fit)]
    
    if (length(pars_exist) > 0) {
      p <- traceplot(fit, pars = pars_exist, inc_warmup = FALSE) +
        ggtitle(paste("Traceplots:", group_name))
      
      png(
        filename = file.path(plot_dir, paste0("trace_", group_name, ".png")),
        width = 10,
        height = 2 + 2 * length(pars_exist),
        units = "in",
        res = 350
      )
      print(p)
      dev.off()
      
      message(sprintf("Saved: trace_%s.png", group_name))
    }
  }
  
  message(sprintf("Traceplots saved to: %s", plot_dir))
  
  invisible(plot_dir)
}
####################################################################################
plot_health_regions <- function(map_data, 
                                show_numbers = TRUE, 
                                base_size = 14,
                                title = "Health Regions of Thailand") {
    colors_pallete <- c(
      "Health Region 1" = "#F06292",
      "Health Region 2" = "#F49BAB",
      "Health Region 3" = "#DA6C6C",
      "Health Region 4" = "#5C6BC0",
      "Health Region 5" = "#42A5F5",
      "Health Region 6" = "#4DD0E1",
      "Health Region 7" = "#4DB6AC",
      "Health Region 8" = "#80CBC4",
      "Health Region 9" = "#81C784",
      "Health Region 10" = "#AED581",
      "Health Region 11" = "#FF8A65",
      "Health Region 12" = "#FFD54F",
      "Health Region 13" = "#BA68C8"
    )
  
  map_data <- map_data %>%
    mutate(Zone = paste("Health Region", str_remove(Zone, "Zone0?")))
  
  # Calculate centroids for labels
  zone_centroids <- map_data %>%
    group_by(Zone) %>%
    summarise(geometry = st_union(geometry), .groups = 'drop') %>%
    st_centroid() %>%
    mutate(
      coords = st_coordinates(geometry),
      lon = coords[, 1],
      lat = coords[, 2],
      zone_number = gsub("Health Region ", "", Zone)
    )
  
  # Base plot
  p <- ggplot() +
    geom_sf(data = map_data, aes(fill = Zone), colour = "white", linewidth = 0.2) +
    annotation_scale(location = "bl",
                     width_hint = 0.2,
                     bar_cols = c("black", "white")) +
    annotation_north_arrow(location = "bl",
                           which_north = "true",
                           pad_y = unit(0.3, "in"),
                           height = unit(0.8, "cm"),
                           width = unit(0.8, "cm"),
                           style = north_arrow_fancy_orienteering) +
    coord_sf(xlim = c(96, 106), ylim = c(5, 21), expand = FALSE) +
    scale_fill_manual(values = colors_pallete, guide = "none") +  
    theme_void(base_size = base_size) +
    labs(x = "", y = "", title = title) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.margin = margin(1, 1, 1, 1, unit = "pt")
    )
  
  # Add numbers with white border
  if (show_numbers) {
    p <- p +
      geom_text(data = zone_centroids, 
                aes(x = lon, y = lat, label = zone_number),
                size = 7, fontface = "bold", color = "white", vjust = 0.5) +
      geom_text(data = zone_centroids, 
                aes(x = lon, y = lat, label = zone_number),
                size = 6.5, fontface = "bold", color = "black", vjust = 0.5)
  }
  
  return(p)
}
##################################################################################
plot_area_cases_map <- function(case_data, map_data) {
  area_summary <- case_data |>
    group_by(Area) |>
    summarize(
      cases = sum(case, na.rm = TRUE)
    )

  map_joined <- map_data |>
    left_join(area_summary, by = c("Provinc" = "Area"))

  ggplot() +
    geom_sf(data = map_joined, aes(fill = cases), color = "white", size = 0.15) +
    scale_fill_gradientn(
      colours = c('#DEDED1','#94B4C1','#547792','#213448'),
      na.value = "grey90",
      labels = scales::comma,
      name = "Cumulative reported cases\n(2003-2024)",
      limits = c(0, 18000)
    ) +
    coord_sf(xlim = c(96, 106), ylim = c(5, 21), expand = FALSE) +
    theme_void(base_size = 14) +
    theme(
      legend.position = c(0.80, 0.25),
      legend.title = element_text(size = 11, lineheight = 0.95),
      legend.text = element_text(size = 10),
      legend.key.height = unit(0.55, "cm"),
      legend.key.width = unit(0.35, "cm"),
      legend.spacing.y = unit(0.1, "cm"),
      plot.margin = margin(1, 1, 1, 1, unit = "pt")
    )
}
##################################################################################
plot_area_deaths_map <- function(case_data, map_data) {
  area_summary <- case_data |>
    group_by(Area) |>
    summarize(
      deaths = sum(death, na.rm = TRUE)
    )

  map_joined <- map_data |>
    left_join(area_summary, by = c("Provinc" = "Area"))

  ggplot() +
    geom_sf(data = map_joined, aes(fill = deaths), color = "white", size = 0.15) +
    scale_fill_gradientn(
      colours = c('#DEDED1','#B6AE9F','#D34E4E','#BF1A1A'),
      na.value = "grey90",
      labels = scales::comma,
      name = "Cumulative reported deaths\n(2003-2024)",
      limits = c(0, 30)
    ) +
    coord_sf(xlim = c(96, 106), ylim = c(5, 21), expand = FALSE) +
    theme_void(base_size = 14) +
    theme(
      legend.position = c(0.80, 0.25),
      legend.title = element_text(size = 11, lineheight = 0.95),
      legend.text = element_text(size = 10),
      legend.key.height = unit(0.55, "cm"),
      legend.key.width = unit(0.35, "cm"),
      legend.spacing.y = unit(0.1, "cm"),
      plot.margin = margin(1, 1, 1, 1, unit = "pt")
    )
}
##################################################################################
plot_case_bar <- function(case_data, base_size = 14) {
  plot_data <- case_data %>%
    group_by(year) %>%
    summarise(case = sum(case), .groups = "drop")

  plot_data %>%
    ggplot(aes(x = year)) +
    geom_bar(aes(y = case), stat = "identity", fill = "#213448", alpha = 0.8, width = 0.95) +
    theme_bw(base_size = base_size) +
    labs(x = "", y = "") +
    scale_x_continuous(breaks = seq(2003, 2024, by = 2)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(labels = comma) +
    ggtitle("Annual reported cases")
}
##################################################################################
plot_death_bar <- function(case_data, base_size = 14) {
  plot_data <- case_data %>%
    group_by(year) %>%
    summarise(death = sum(death), .groups = "drop")

  plot_data %>%
    ggplot(aes(x = year)) +
    geom_bar(aes(y = death), stat = "identity", fill = "#BF1A1A", alpha = 0.8, width = 0.95) +
    theme_bw(base_size = base_size) +
    labs(x = "", y = "") +
    scale_x_continuous(breaks = seq(2003, 2024, by = 2)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(labels = comma) +
    ggtitle("Annual reported deaths")
}
##################################################################################
plot_agegroup_proportion <- function(case_age_data) {
  age_levels <- c("0 to 5",  "5 to 15", "15 to 25", "25 to 35", "35 to 45", "45 to 55", "55 to 65", "above 65", "unknown")
  age_summary <- case_age_data %>%
    group_by(age_group_new) %>%
    summarise(
      cases = sum(case, na.rm = TRUE),
      deaths = sum(death, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(cases, deaths), names_to = "type", values_to = "count") %>%
    group_by(type) %>%
    mutate(proportion = count / sum(count)) %>%
    mutate(age_group_new = factor(age_group_new, levels = age_levels))
  
  ggplot(age_summary, aes(x = age_group_new, y = proportion, fill = type)) +
    geom_col(position = "dodge", alpha = 0.8) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "", y = "", fill = "", title = "Age distribution of cases and deaths") +
    scale_fill_manual(values = c("cases" = "#213448", "deaths" = "#BF1A1A")) +
    theme_bw(base_size = 14) +
    coord_flip() +
    theme(legend.position = c(0.85, 0.3))
}


