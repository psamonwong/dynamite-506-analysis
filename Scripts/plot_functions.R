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


#############################################################################
extract_lambda_base <- function(probs = c(0.025, 0.5, 0.975)) {
  stopifnot(exists("fit"), exists("stan_data"), exists("site_levels"), exists("calendar"))
  S <- stan_data$S
  Y <- stan_data$Y
  R <- stan_data$R

  # Extract lambda_base_ys draws [Y, S]
  lambda_mat <- as.matrix(fit, pars = "lambda_base_ys")
  if (is.null(colnames(lambda_mat))) stop("lambda_base_ys columns have no names; re-run sampling with pars='lambda_base_ys'.")

  # Parse indices from column names [y,s]
  idx <- stringr::str_match(colnames(lambda_mat), "lambda_base_ys\\[(\\d+),(\\d+)\\]")
  if (any(is.na(idx))) stop("Failed to parse lambda_base_ys indices from column names.")
  y_idx_col <- as.integer(idx[, 2])
  s_idx_col <- as.integer(idx[, 3])

  # Create parameter map
  param_map <- tibble::tibble(
    par = colnames(lambda_mat),
    y_idx = y_idx_col,
    s_idx = s_idx_col
  )

  # Convert to long format
  lambda_long <- tibble::as_tibble(lambda_mat) |>
    dplyr::mutate(.draw = dplyr::row_number()) |>
    tidyr::pivot_longer(cols = -.draw, names_to = "par", values_to = "lambda_base") |>
    dplyr::left_join(param_map, by = "par")

  # Summarize by year and province
  lambda_summ <- lambda_long |>
    dplyr::group_by(y_idx, s_idx) |>
    dplyr::summarize(
      mean = mean(lambda_base),
      q025 = 12*stats::quantile(lambda_base, probs[1], names = FALSE),
      q500 = 12*stats::quantile(lambda_base, probs[2], names = FALSE),
      q975 = 12*stats::quantile(lambda_base, probs[3], names = FALSE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      province = site_levels[s_idx],
      year = min(calendar$year) + y_idx - 1  # Map y_idx to actual year
    )

  # Extract region info from stan_data
  region_id <- stan_data$region_id  # [S] province to region mapping
  region_map <- tibble::tibble(
    s_idx = seq_len(S),
    r_idx = region_id[s_idx],
    region = paste0("Zone", sprintf("%02d", region_id[s_idx]))
  )
  
  lambda_summ <- lambda_summ |>
    dplyr::left_join(region_map, by = "s_idx")

  lambda_summ
}
################################################################################
  # Helper: convert Zone to Health Region
  zone_to_hr <- function(z) {
    z_num <- as.integer(sub("Zone0?", "", z))
    if (z_num %in% c(4, 13)) return("Health Region 4 and 13")
    paste0("Health Region ", z_num)
  }


plot_foi_trends_inc_groups <- function(lambda_data,
                                      year0 = 2003,
                                      end_year = 2024,
                                       inc_breaks = c(-Inf, 100, 1000, Inf),
                                       inc_labels = c("< 100", "100-1000", "> 1000"),
                                       highlight_combos = c("100-1000-> 1000", "> 1000-> 1000", "> 1000-> 100-1000"),
                                       grey_combo = "Other") {
  stopifnot(length(inc_breaks) == length(inc_labels) + 1)

  # Province abbreviation lookup
  abbr_lookup <- c(
    "Amnat Charoen" = "ACR", "Ang Thong" = "ATG", "Bangkok" = "BKK",
    "Bueng Kan" = "BKN", "Bungkan" = "BKN",
    "Buri Ram" = "BRM", "Chachoengsao" = "CCO",
    "Chai Nat" = "CNT", "Chaiyaphum" = "CPM", "Chanthaburi" = "CTI",
    "Chiang Mai" = "CMI", "Chiang Rai" = "CRI", "Chon Buri" = "CBI",
    "Chumphon" = "CPN", "Kalasin" = "KSN", "Kamphaeng Phet" = "KPT",
    "Kanchanaburi" = "KRI", "Khon Kaen" = "KKN", "Krabi" = "KBI",
    "Lampang" = "LPG", "Lamphun" = "LPN", "Loei" = "LEI",
    "Lopburi" = "LRI", "Lop Buri" = "LRI",
    "Mae Hong Son" = "MSN", "Maha Sarakham" = "MKM",
    "Mukdahan" = "MDH", "Nakhon Nayok" = "NYK", "Nakhon Pathom" = "NPT",
    "Nakhon Phanom" = "NPM", "Nakhon Ratchasima" = "NMA",
    "Nakhon Sawan" = "NSN", "Nakhon Si Thammarat" = "NRT",
    "Nan" = "NAN", "Narathiwat" = "NWT", "Nong Bua Lam Phu" = "NBP",
    "Nong Khai" = "NKI", "Nonthaburi" = "NBI", "Pathum Thani" = "PTE",
    "Pattani" = "PTN",
    "Phang Nga" = "PNA", "Phangnga" = "PNA",
    "Phatthalung" = "PLG",
    "Phayao" = "PYO", "Phetchabun" = "PNB", "Phetchaburi" = "PBI",
    "Phichit" = "PCT", "Phitsanulok" = "PLK",
    "P.Nakhon S.Ayutthaya" = "AYA", "Phra Nakhon Si Ayutthaya" = "AYA",
    "Phrae" = "PRE", "Phuket" = "PKT",
    "Prachinburi" = "PRI", "Prachin Buri" = "PRI",
    "Prachuap Khiri Khan" = "PKN", "Ranong" = "RNG", "Ratchaburi" = "RBR",
    "Rayong" = "RYG", "Roi Et" = "RET", "Sa Kaeo" = "SKW",
    "Sakon Nakhon" = "SNK", "Samut Prakan" = "SPK", "Samut Sakhon" = "SKN",
    "Samut Songkhram" = "SKM", "Saraburi" = "SRI", "Satun" = "STN",
    "Sing Buri" = "SBR", "Si Sa Ket" = "SSK", "Songkhla" = "SKA",
    "Sukhothai" = "STI", "Suphan Buri" = "SPB", "Surat Thani" = "URT",
    "Surin" = "SRN", "Tak" = "TAK", "Trang" = "TRG", "Trat" = "TRA",
    "Ubon Ratchathani" = "UBN", "Udon Thani" = "UDN",
    "Uthai Thani" = "UTI", "Uttaradit" = "UTD", "Yala" = "YLA",
    "Yasothon" = "YST"
  )

  region_levels <- c(
    "Health Region 1", "Health Region 2", "Health Region 3", "Health Region 4 and 13",
    "Health Region 5", "Health Region 6", "Health Region 7", "Health Region 8",
    "Health Region 9", "Health Region 10", "Health Region 11", "Health Region 12"
  )

  lambda_inc <- lambda_data |>
    mutate(
      region = sapply(region, zone_to_hr),
      foi_ann = q500,
      inc_100k = (1 - exp(-foi_ann)) * 100000
    )

  # end_year <- max(lambda_inc$year, na.rm = TRUE)
  # year0 <- min(lambda_inc$year, na.rm = TRUE)

 
class_tbl <- lambda_inc |>
    filter(year %in% c(year0, end_year)) |>
    transmute(
      region, province, year,
      inc_val = inc_100k,
      inc_group = cut(inc_val, breaks = inc_breaks, labels = inc_labels, include.lowest = TRUE, right = TRUE)
    ) |>
  mutate(year = ifelse(year == year0, "y0", "y1")) |>
    tidyr::pivot_wider(
      names_from = year,
      values_from = c(inc_val, inc_group),
      names_sep = "_"
    ) |>
    mutate(
      inc_group_combo = paste0(inc_group_y0, "-", inc_group_y1),
      is_highlight = inc_group_combo %in% highlight_combos,
      inc_group_combo_plot = inc_group_combo#ifelse(is_highlight, inc_group_combo, grey_combo)
    )


  plot_df <- lambda_inc |>
    left_join(class_tbl, by = c("region", "province")) |>
    mutate(
      region = factor(region, levels = region_levels),
      inc_group_combo_plot = factor(inc_group_combo_plot)
    )

  label_df <- class_tbl |>
    filter(is_highlight) |>
    mutate(
      province_abbr = abbr_lookup[province],
      region = factor(region, levels = region_levels)
    )

  pals <- expand.grid(inc_labels, inc_labels) |>
    dplyr::mutate(inc_group_combo = paste0(Var1, "-", Var2))

  group_labels <- pals$inc_group_combo

  core_cols <- c(
    "#e8e8e8", "#ace4e4", "#5ac8c8",
    "#dfb0d6", "#a5add3", "#5698b9",
    "#be64ac", "#8c62aa", "#3b4994"
  )

  pal_cols <- setNames(core_cols[seq_along(group_labels)], group_labels)

  line_size <- ifelse(plot_df$inc_group_combo_plot %in% highlight_combos, 1, 0.3)
  line_alpha <- ifelse(plot_df$inc_group_combo_plot %in% highlight_combos, 0.85, 0.85)

  ggplot(plot_df, aes(year, foi_ann, group = province)) +
    geom_line(
      aes(linewidth = I(line_size), alpha = I(line_alpha)),
      colour = "grey60"
    ) +
    geom_line(
      data = dplyr::filter(plot_df, inc_group_combo_plot != grey_combo),
      aes(colour = inc_group_combo_plot),
      linewidth = 1.2,
      alpha = 1
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(
        x = end_year,
        y = -log1p(-pmin(inc_val_y1, 99999) / 100000),
        label = province_abbr,
        colour = inc_group_combo_plot
      ),
      direction = "y",
      hjust = 0,
      nudge_x = 1,
      size =4,
      fontface = "bold",
      segment.size = 0.3,
      segment.color = "grey50",
      segment.linetype = "dotted",
      min.segment.length = 0,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    facet_wrap(~ region) +
    scale_y_log10(
      name = "Estimated annual force of infection\n(per at-risk person-year)",
      labels = scales::label_scientific(digits = 1),
      sec.axis = sec_axis(
        transform = ~ (1 - exp(-.)) * 100000,
        name = "Estimated annual incidence\n(per 100,000 at-risk population)",
        labels = scales::label_comma(accuracy = 1)
      )
    ) +
    scale_colour_manual(values = pal_cols, name = NULL) +
    scale_fill_manual(values = pal_cols, name = NULL) +
    coord_cartesian(clip = "off") +
    scale_x_continuous(
      expand = expansion(mult = c(0.02, 0.2)),
      breaks = seq(min(plot_df$year, na.rm = TRUE), max(plot_df$year, na.rm = TRUE), 5)
    ) +
    theme_bw(base_size = 16) +
    theme(
      panel.grid.major = element_line(linewidth = 0.4),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      panel.spacing = unit(0.5, "lines"),
      plot.margin = margin(5.5, 50, 5.5, 5.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank()
    ) +
    labs(x = NULL) +
    ggtitle("Temporal trends in estimated force of infection and annual incidence") 
}

######################################################################################
make_bivar_map <- function(lambda_data,
                                              map_data,
                                              year0 = 2003,
                                              year1 = 2024,
                                              foi_breaks = c(-Inf, 100, 1000, Inf),
                                              foi_labels = c("< 100", "100-1000", "> 1000"),
                                              pal = "DkBlue",
                                              subtitle = "Bivariate choropleth: 2003 vs 2024 (median FOI)",
                                              legend_bbox = c(xmin = 101.0, xmax = 104.3, ymin = 8.0, ymax = 11.2)) {
  stopifnot(length(foi_breaks) == length(foi_labels) + 1)

  foi_wide <- lambda_data |>
    mutate(
      foi_ann = q500,
      inc_100k = (1 - exp(-foi_ann)) * 100000
    ) |>
    filter(year %in% c(year0, year1)) |>
    select(province, year, inc_100k) |>
    tidyr::pivot_wider(names_from = year, values_from = inc_100k, names_prefix = "y") |>
    dplyr::rename(inc_0 = !!paste0("y", year0), inc_1 = !!paste0("y", year1))

  foi_wide <- foi_wide |>
    mutate(
      inc_0_cat = cut(inc_0, breaks = foi_breaks, labels = foi_labels, include.lowest = TRUE, right = TRUE),
      inc_1_cat = cut(inc_1, breaks = foi_breaks, labels = foi_labels, include.lowest = TRUE, right = TRUE),
      bi_ix_x = as.integer(inc_0_cat),
      bi_ix_y = as.integer(inc_1_cat),
      bi_class = paste0(bi_ix_x, "-", bi_ix_y)
    )

  map_bi <- map_data |>
    left_join(foi_wide, by = c("Provinc" = "province")) |>
    mutate(
      bi_class = factor(
        bi_class,
        levels = as.vector(outer(seq_along(foi_labels), seq_along(foi_labels), \(x, y) paste0(x, "-", y)))
      )
    )

  
  p_map <- ggplot(map_bi) +
    geom_sf(aes(fill = bi_class), color = "white", linewidth = 0.3, show.legend = FALSE) +
    biscale::bi_scale_fill(pal = pal, dim = length(foi_labels), na.value = "grey90") +
    coord_sf(xlim = c(96, 106), ylim = c(5, 21), expand = FALSE) +
    theme_void(base_size = 14) +
    labs(title = subtitle) +
    theme(
      plot.title = element_blank(),
      axis.title = element_blank(),
      plot.margin = margin(0.1, 0.1, 0.1, 0.1)
    )

  options(error = NULL)
legend_plot <- biscale::bi_legend(
    pal = pal,
    dim = length(foi_labels),
    xlab = paste0("Incidence in ", year0),
    ylab = paste0("Incidence in ", year1),
  ) +
    scale_x_continuous(breaks = seq_along(foi_labels),
                       labels = foi_labels,
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = seq_along(foi_labels),
                       labels = foi_labels,
                       expand = c(0, 0)) +
    theme(
      axis.text.x  = element_text(size = 11, angle = 45, hjust = 1),
      axis.text.y  = element_text(size = 11),
      axis.title = element_text(size = 13),
      plot.margin = margin(0, 0, 0, 0)
    )

  # stack vertically
  p_map / legend_plot + patchwork::plot_layout(heights = c(5, 1))
}
######################################################################################
eval_seasonality_curves <- function(
  fit, H, P, S,
  entity = c("province", "region"),
  region_id = NULL,
  region_labels = NULL,
  month_grid = seq(1, 12, length.out = 241),
  probs = c(0.025, 0.5, 0.975)
) {
  entity <- match.arg(entity)

  dr <- rstan::extract(fit, pars = c("a_cos", "a_sin"))
  a_cos <- dr$a_cos  # iter x H x P x S
  a_sin <- dr$a_sin  # iter x H x P x S

  stopifnot(length(dim(a_cos)) == 4, dim(a_sin)[1] == dim(a_cos)[1])

  n_iter <- dim(a_cos)[1]
  theta <- 2 * pi * month_grid / 12

  # Precompute basis: (grid x H)
  cos_basis <- sapply(seq_len(H), \(h) cos(h * theta))
  sin_basis <- sapply(seq_len(H), \(h) sin(h * theta))

  # Result: iter x grid x P x S
  f <- array(0, dim = c(n_iter, length(month_grid), P, S))
for (p in seq_len(P)) {
  for (s in seq_len(S)) {
    A_cos <- a_cos[, , p, s, drop = FALSE][, , 1, 1]  # n_iter x H
    A_sin <- a_sin[, , p, s, drop = FALSE][, , 1, 1]  # n_iter x H

    f[, , p, s] <- A_cos %*% t(cos_basis) + A_sin %*% t(sin_basis)
  }
}

  mult <- exp(f) # iter x grid x P x S

  # Optionally aggregate provinces -> regions for each draw
  if (entity == "region") {
    stopifnot(!is.null(region_id))
    R <- max(region_id)
    mult_reg <- array(NA_real_, dim = c(n_iter, length(month_grid), P, R))
    for (r in seq_len(R)) {
      s_idx <- which(region_id == r)
      mult_reg[, , , r] <- apply(mult[, , , s_idx, drop = FALSE], c(1, 2, 3), mean)
    }
    mult <- mult_reg
    S <- dim(mult)[4]

    if (is.null(region_labels)) region_labels <- paste0("Region ", seq_len(S))
    entity_names <- region_labels
  } else {
    entity_names <- paste0("Province ", seq_len(S))
  }

  # Summarize across iterations
  out <- tidyr::crossing(
    x = month_grid,
    period = factor(seq_len(P), labels = paste0("Period ", seq_len(P))),
    entity = factor(entity_names, levels = entity_names)
  ) |>
    mutate(q025 = NA_real_, q500 = NA_real_, q975 = NA_real_)

  # Create index matching
  out <- out |>
    mutate(
      idx = match(
        paste(x, period, entity, sep = "|"),
        paste(month_grid[match(x, month_grid)], period, entity, sep = "|")
      )
    )
  
  for (e in seq_len(S)) {
    for (p in seq_len(P)) {
      idx_mask <- out$period == levels(out$period)[p] & out$entity == entity_names[e]
      qs_all <- apply(mult[, , p, e], 2, \(col) stats::quantile(col, probs = probs, names = FALSE))
      out$q025[idx_mask] <- qs_all[1, ]
      out$q500[idx_mask] <- qs_all[2, ]
      out$q975[idx_mask] <- qs_all[3, ]
    }
  }

  out
}


compute_seasonality_indices <- function(season_curves) {
  season_curves |>
    rename(region = entity, month = x, mult = q500) |>
    group_by(region, period) |>
    summarize(
      peak_month = month[which.max(mult)],
      trough_month = month[which.min(mult)],
      peak_mult = max(mult),
      trough_mult = min(mult),
      amp_ratio = peak_mult / trough_mult,
      amp_log = log(amp_ratio),
      .groups = "drop"
    )
}

cluster_seasonality <- function(indices, k, 
  exclude_cols = c("region", "period", "trough_month", "trough_mult", "amp_ratio")) {
  selected <- indices |> select(-all_of(exclude_cols))
  km <- kmeans(selected, centers = k, nstart = 25)
  # Order clusters by frequency
  freq <- sort(table(km$cluster), decreasing = TRUE)
  new_labels <- match(km$cluster, as.integer(names(freq)))
  pattern_labels <- paste0("Pattern ", seq_along(freq))
  indices$clusters <- factor(new_labels, levels = seq_along(freq), labels = pattern_labels)
  indices
}

plot_seasonality_by_cluster <- function(season_curves, indices) {
  cluster_summaries <- indices |>
    group_by(clusters) |>
    summarize(
      median_peak_month = median(peak_month),
      median_amp_ratio = median(amp_ratio),
      .groups = "drop"
    ) |>
    mutate(
      label = paste0(
        "Peak month: ", formatC(median_peak_month, format = "f", digits = 1), "\n",
        "Amplitude ratio: ", formatC(median_amp_ratio, format = "f", digits = 1), "\n"
      ),
      x = 1, y = 2.0
    )

  season_curves %>%
    left_join(indices %>% select(region, period, clusters), by = c("entity" = "region", "period")) %>%
    group_by(x, clusters) %>%
    summarize(q_med = median(q500),
              q_low = quantile(q500, 0.25),
              q_up = quantile(q500, 0.77),
              .groups = "drop_last") %>%
    ggplot(aes(x = x)) +
    geom_ribbon(aes(ymin = q_low, ymax = q_up, fill = as.factor(clusters)), alpha = 0.3, color = NA) +
    geom_line(linewidth = 1.5, aes(y = q_med, color = as.factor(clusters))) +
    geom_hline(yintercept = 1, linetype = "31", alpha = 0.75) +
    facet_wrap(~ clusters, ncol = 3) +
    geom_text(data = cluster_summaries, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, hjust = 0, vjust = 1, size = 4) +
    theme_bw(base_size = 14) +
    scale_x_continuous(breaks = seq(1, 12, by = 3), minor_breaks = NULL,
                       labels = month.abb[seq(1, 12, by = 3)]) +
    scale_y_continuous(breaks = seq(0.5, 3, by = 0.5)) +
    theme(legend.position = "none") +
    xlab("") + ylab("Seasonal multiplier") +
    scale_color_manual(values = c("#F88A40", "#3E8359", "#5B96C5")) +
    scale_fill_manual(values = c("#F88A40", "#3E8359", "#5B96C5"))
}

elbow_plot <- function(indices, max_k = 10, exclude_cols = c("region", "period", "clusters", "trough_month", "trough_mult", "amp_ratio")) {
  selected <- indices |> select(-all_of(exclude_cols))
  wss <- sapply(1:max_k, function(k) {
    kmeans(selected, centers = k, nstart = 25)$tot.withinss
  })
  tibble::tibble(k = 1:max_k, wss = wss) %>%
    ggplot(aes(x = k, y = wss)) +
    geom_point(size = 2) +
    geom_line() +
    labs(
      x = "Number of clusters (k)",
      y = "Total within-cluster sum of squares",
      title = "Elbow Plot for K-means Clustering"
    ) +
    theme_bw(base_size = 14)
}
###############################################################################
plot_seasonality_map <- function(seasonality_indices, map_data, 
  site_levels, period_i) {
  seasonality_indices <- seasonality_indices %>%
    mutate(region = factor(region, levels = paste0("Province ", seq_along(site_levels)), labels = site_levels))
  
  seasonality_indices_i <- seasonality_indices %>%
    filter(period == period_i)
  
  map_seasonality_i <- map_data %>%
    left_join(
      seasonality_indices_i %>%
        rename(geometry_id = region),
      by = c("Provinc" = "geometry_id")
    )
  
  ggplot(map_seasonality_i) +
    geom_sf(aes(fill = as.factor(clusters)), color = "white", linewidth = 0.1, alpha = 0.85) +
    scale_fill_manual(name = NULL, values = c("#F88A40", "#3E8359", "#5B96C5")) +
    coord_sf(xlim = c(96, 106), ylim = c(5, 21), expand = FALSE) +
    theme_void() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none") 
}

#################################################################################
plot_seasonality_alluvial <- function(seasonality_indices,
                                      period_labels = c("2003-2007", "2008-2012", "2013-2018", "2019-2024"),
                                      cluster_colors =c("#F88A40", "#3E8359", "#5B96C5")) {
  library(dplyr)
  library(tidyr)
  library(ggalluvial)
  library(ggplot2)

  # Prepare data
  seasonality_trend <- seasonality_indices %>%
    select(region, period, clusters)

  # Wide format: one column per period's cluster
  seasonality_trend_wide <- seasonality_trend %>%
    mutate(period = as.character(period)) %>%
    pivot_wider(
      id_cols = region,
      names_from = period,
      values_from = clusters,
      names_prefix = "cluster_"
    ) %>%
    mutate(freq = 1)

  # Plot
  ggplot(seasonality_trend_wide,
         aes(axis1 = `cluster_Period 1`,
             axis2 = `cluster_Period 2`,
             axis3 = `cluster_Period 3`,
             axis4 = `cluster_Period 4`,
             y = freq)) +
    geom_alluvium(aes(fill = `cluster_Period 4`), width = 1/3, alpha = 0.7) +
    geom_stratum(width = 1/3, fill = "grey90", color = "black") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5) +
    scale_x_continuous(
      breaks = 1:4,
      labels = period_labels,
      expand = c(0, 0)
    ) +
    labs(x = "Period", y = "", fill = "Cluster (Period 1)") +
    scale_fill_manual(values = cluster_colors) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(color = "black"),
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0)
    )
}
################################################################################
# Extract and plot Se_y over time
plot_Se_over_time <- function(probs = c(0.025, 0.5, 0.975), map_data) {
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
        axis.title.x = element_blank()) +
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
    theme(plot.caption = element_text(hjust = 0.5, vjust = 10, size = 12, face = "bold")) + 
     plot_annotation(theme = theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5)))


  combined_se_plot <- se_trend + map_se + plot_layout(ncol = 2, widths  = c(2, 1)) 
  combined_se_plot + plot_annotation(theme = theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2)))
  
}
################################################################################

plot_relative_susceptibility <- function(fit, age_levels, ref_age = "25 to 35") {
  rho_mat <- as.matrix(fit, pars = "rho_pa")
  idx <- stringr::str_match(colnames(rho_mat), "rho_pa\\[(\\d+),(\\d+)\\]")
  
  rho_summ <- tibble(par = colnames(rho_mat), p_idx = as.integer(idx[, 2]), a_idx = as.integer(idx[, 3])) |>
    left_join(
      as_tibble(rho_mat) |> mutate(.draw = row_number()) |>
        pivot_longer(cols = -.draw, names_to = "par", values_to = "rho"),
      by = "par"
    ) |>
    group_by(p_idx, a_idx) |>
    summarize(
      mean = mean(rho), q025 = quantile(rho, 0.025), q975 = quantile(rho, 0.975),
      .groups = "drop"
    ) |>
    mutate(
      period = paste0("Period ", p_idx),
      age_group = age_levels[a_idx]
    )
  
  rho_summ |>
    mutate(
      age_group = factor(age_group, levels = age_levels),
      period = factor(period, levels = paste0("Period ", 1:4), 
                     labels = c("2003 to 2007", "2008 to 2012", "2013 to 2018", "2019 to 2024"))
    ) |>
    filter(age_group != ref_age) |>
    ggplot() +
    geom_hline(yintercept = 1, linetype = "21", color = "gray30", linewidth = 1, alpha = 0.7) +
    geom_point(aes(x = age_group, y = mean, col = period), alpha = 0.8,
                  position = position_dodge(width = 0.4), size = 3) +
    geom_errorbar(aes(x = age_group, ymin = q025, ymax = q975, col = period), 
                  position = position_dodge(width = 0.4), alpha = 0.8,
                  width = 0, linewidth = 1.1) +
     scale_color_viridis_d(option = "magma", end = 0.8, name="") +
    labs(
      x = "",
      y = "Relative Susceptibility",
      title = paste0("Relative Susceptibility\n(relative to ", ref_age, " years)")) +
    theme_bw(base_size = 14) +
    scale_y_continuous(limits = c(0.2, 1.8)) +
    coord_flip()
}
################################################################################

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
    labs(x = "", y = "Case-Fatality Rate (%)", title = "\nCase-Fatality Rate") +
    theme_bw(base_size = 14) +
    coord_flip()
}
