export_traceplots <- function(dataset) {
  library(here)
  library(rstan)
  library(ggplot2)
  
  # Load fitted model
  load(here("Rout", paste0("fit_", dataset, ".RData")))
  
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
      
      ggsave(
        filename = file.path(plot_dir, paste0("trace_", group_name, ".png")),
        plot = p,
        width = 10,
        height = 2 + 2 * length(pars_exist),
        dpi = 150
      )
      
      message(sprintf("Saved: trace_%s.png", group_name))
    }
  }
  
  message(sprintf("Traceplots saved to: %s", plot_dir))
  
  invisible(plot_dir)
}
