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
                                base_size = 12,
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
    annotation_scale(location = "br",      
                     width_hint = 0.3,     
                     bar_cols = c("black", "white")) +
    annotation_north_arrow(location = "br",
                           which_north = "true",
                           pad_x = unit(0.2, "in"),
                           pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering) +
    scale_fill_manual(values = colors_pallete, guide = "none") +  
    theme_minimal(base_size = base_size) +
    labs(x = "", y = "", title = title) +
    theme(plot.title = element_text(face = "bold"))
  
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
################################################################################
plot_case_death <- function(case_data, base_size = 14) {
  plot_data <- case_data %>%
    group_by(year) %>%
    summarise(
      case = sum(case),
      death = sum(death))

  plot_case <- plot_data %>%
    ggplot(aes(x = year)) +
    geom_bar(aes(y = case), stat = "identity", fill = "gray30", alpha = 0.8, width = 0.95) +
    theme_bw(base_size = base_size) +
    labs(x = "", y = "") +
    scale_x_continuous(breaks = seq(2003, 2024, by = 2)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(labels = comma) +
    ggtitle("Reported cases")

  plot_death <- plot_data %>%
    ggplot(aes(x = year)) +
    geom_bar(aes(y = death), stat = "identity", fill = "#C3110C", alpha = 0.8, width = 0.95) +
    theme_bw(base_size = base_size) +
    labs(x = "", y = "") +
    scale_x_continuous(breaks = seq(2003, 2024, by = 2)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(labels = comma) +
    ggtitle("Reported deaths")

  ggarrange(plot_case, plot_death, ncol = 1, nrow = 2, align = "v")
}






