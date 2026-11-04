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
# Average annual reported incidence (per 100,000) and mortality (per 1,000,000)
# by province, using the model's population data as the denominator.
# Province-years without population data (Bungkan, 2011) are excluded from
# both numerator and denominator.
summarise_reported_rates <- function(case_data, population_data, by = c("Area", "Zone")) {
  pop_year <- population_data |>
    group_by(Province, year) |>
    summarise(population = sum(Total), .groups = "drop")

  case_data |>
    group_by(Area, Zone, year) |>
    summarise(cases = sum(case), deaths = sum(death), .groups = "drop") |>
    left_join(pop_year, by = c("Area" = "Province", "year")) |>
    filter(!is.na(population)) |>
    group_by(across(all_of(by))) |>
    summarise(
      cases = sum(cases),
      deaths = sum(deaths),
      person_years = sum(population),
      .groups = "drop"
    ) |>
    mutate(
      incidence_100k = cases / person_years * 1e5,
      mortality_1m = deaths / person_years * 1e6
    )
}

# Directly age-standardised reported incidence and mortality by province,
# using the pooled national population over the same years as the standard.
summarise_age_standardised_rates <- function(case_age_data, population_data) {
  age_map <- c("0 to 4" = "0 to 5", "5 to 14" = "5 to 15", "15 to 24" = "15 to 25",
               "25 to 34" = "25 to 35", "35 to 44" = "35 to 45", "45 to 54" = "45 to 55",
               "55 to 64" = "55 to 65", "above 65" = "above 65")
  years <- range(case_age_data$year)

  pop_age <- population_data |>
    filter(between(year, years[1], years[2])) |>
    mutate(age_group_new = age_map[age_group]) |>
    group_by(Province, age_group_new) |>
    summarise(person_years = sum(Total), .groups = "drop")

  standard <- pop_age |>
    group_by(age_group_new) |>
    summarise(weight = sum(person_years), .groups = "drop") |>
    mutate(weight = weight / sum(weight))

  case_age_data |>
    filter(age_group_new != "unknown") |>
    group_by(Area, age_group_new) |>
    summarise(cases = sum(case), deaths = sum(death), .groups = "drop") |>
    left_join(pop_age, by = c("Area" = "Province", "age_group_new")) |>
    left_join(standard, by = "age_group_new") |>
    group_by(Area) |>
    summarise(
      asr_incidence_100k = sum(weight * cases / person_years) * 1e5,
      asr_mortality_1m = sum(weight * deaths / person_years) * 1e6,
      .groups = "drop"
    )
}

plot_rate_map <- function(rates, map_data, value, name, colours, breaks) {
  map_joined <- map_data |>
    left_join(rates, by = c("Provinc" = "Area"))

  ggplot() +
    geom_sf(data = map_joined, aes(fill = .data[[value]]), color = "white", size = 0.15) +
    scale_fill_gradientn(
      colours = colours,
      na.value = "grey90",
      trans = "sqrt",
      breaks = breaks,
      labels = scales::label_number(accuracy = 0.1, drop0trailing = TRUE),
      name = name
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

plot_area_cases_map <- function(case_data, population_data, map_data) {
  plot_rate_map(
    summarise_reported_rates(case_data, population_data), map_data,
    value = "incidence_100k",
    name = "Average annual reported\nincidence per 100,000\n(2003-2024)",
    colours = c('#DEDED1','#94B4C1','#547792','#213448'),
    breaks = c(1, 10, 50, 100, 200)
  )
}

plot_area_deaths_map <- function(case_data, population_data, map_data) {
  plot_rate_map(
    summarise_reported_rates(case_data, population_data), map_data,
    value = "mortality_1m",
    name = "Average annual reported\ndeaths per 1,000,000\n(2003-2024)",
    colours = c('#DEDED1','#B6AE9F','#D34E4E','#BF1A1A'),
    breaks = c(0, 0.1, 0.5, 1, 1.5)
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
    labs(x = "", y = "", fill = "", title = "Age distribution of\nreported cases and deaths") +
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
    "Sukhothai" = "STI", "Suphan Buri" = "SPB", "Surat Thani" = "SRT",
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

  # core_cols <- c(
  #   "#e8e8e8", "#ace4e4", "#5ac8c8",
  #   "#dfb0d6", "#a5add3", "#5698b9",
  #   "#be64ac", "#8c62aa", "#3b4994"
  # )

      core_cols <- c(
    "#d3d3d3", "#8aa6c2", "#4279b0",
    "#ba8890", "#7a6b84", "#3a4e78",
    "#9e3547", "#682a41", "#311e3b"
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
        name = "Estimated annual infection incidence\n(per 100,000 at-risk population)",
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
                                              pal = "DkViolet",
                                              flip_axes = FALSE,
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
      bi_class = if (flip_axes) paste0(bi_ix_y, "-", bi_ix_x) else paste0(bi_ix_x, "-", bi_ix_y)
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
    xlab = if (flip_axes) paste0("High incidence in ", year1) else paste0("High incidence in ", year0),
    ylab = if (flip_axes) paste0("High incidence in ", year0) else paste0("High incidence in ", year1),
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

###############################################################################
# Choice of K for seasonal clustering (cluster validity and robustness)

# Same features as cluster_seasonality(): peak_month, peak_mult, amp_log
seasonality_feature_matrix <- function(indices,
  exclude_cols = c("region", "period", "trough_month", "trough_mult", "amp_ratio")) {
  indices |> select(-any_of(c(exclude_cols, "clusters"))) |> as.matrix()
}

adjusted_rand_index <- function(a, b) {
  tab <- table(a, b)
  sum_ij <- sum(choose(tab, 2))
  sum_a <- sum(choose(rowSums(tab), 2))
  sum_b <- sum(choose(colSums(tab), 2))
  expected <- sum_a * sum_b / choose(sum(tab), 2)
  (sum_ij - expected) / ((sum_a + sum_b) / 2 - expected)
}

# Index of the nearest centre (Euclidean) for each row of X
assign_nearest_centre <- function(X, centres) {
  apply(X, 1, function(x) which.min(colSums((t(centres) - x)^2)))
}

# Elbow (WSS), average silhouette width, gap statistic
# and bootstrap stability (ARI between the full-data partition and partitions
# refitted on bootstrap resamples) for each K
cluster_validity_metrics <- function(indices, k_range = 1:8, nstart = 25,
                                     n_boot = 200, gap_B = 100) {
  X <- seasonality_feature_matrix(indices)
  d <- dist(X)
  n <- nrow(X)

  gap_tab <- cluster::clusGap(X, FUNcluster = kmeans, K.max = max(k_range),
                              B = gap_B, nstart = nstart, iter.max = 50)$Tab

  purrr::map_dfr(k_range, function(k) {
    km <- kmeans(X, centers = k, nstart = nstart)
    if (k == 1) {
      sil <- NA_real_; boot_ari <- NA_real_
    } else {
      sil <- mean(cluster::silhouette(km$cluster, d)[, "sil_width"])
      boot_ari <- replicate(n_boot, {
        b <- sample(n, replace = TRUE)
        km_b <- kmeans(X[b, ], centers = k, nstart = 10)
        adjusted_rand_index(km$cluster, assign_nearest_centre(X, km_b$centers))
      })
    }
    tibble::tibble(
      k = k,
      wss = km$tot.withinss,
      silhouette = sil,
      gap = gap_tab[k, "gap"],
      gap_se = gap_tab[k, "SE.sim"],
      boot_ari_median = median(boot_ari),
      boot_ari_q05 = unname(quantile(boot_ari, 0.05, na.rm = TRUE)),
      cluster_sizes = paste(sort(table(km$cluster), decreasing = TRUE), collapse = "/")
    )
  })
}

plot_cluster_validity <- function(metrics, k_selected = 3) {
  metric_labels <- c(
    wss = "A. Elbow: within-cluster SS",
    silhouette = "B. Average silhouette width",
    boot_ari_median = "C. Bootstrap stability (ARI)"
  )
  plot_data <- metrics |>
    select(k, all_of(names(metric_labels)), boot_ari_q05) |>
    pivot_longer(all_of(names(metric_labels)), names_to = "metric", values_to = "value") |>
    mutate(
      lower = ifelse(metric == "boot_ari_median", boot_ari_q05, NA_real_),
      metric = factor(metric, levels = names(metric_labels), labels = metric_labels)
    ) |>
    filter(!is.na(value))

  ggplot(plot_data, aes(x = k, y = value)) +
    geom_vline(xintercept = k_selected, linetype = "31", colour = "grey50") +
    geom_linerange(aes(ymin = lower, ymax = value), colour = "grey60", na.rm = TRUE) +
    geom_line() +
    geom_point(size = 2) +
    facet_wrap(~ metric, scales = "free_y", ncol = 3) +
    scale_x_continuous(breaks = unique(metrics$k), minor_breaks = NULL) +
    labs(x = "Number of clusters (K)", y = NULL) +
    theme_bw(base_size = 14)
}

# Agreement (ARI) of the K-cluster partition with alternative feature encodings
cluster_encoding_sensitivity <- function(indices, k = 3, nstart = 25) {
  X <- seasonality_feature_matrix(indices)
  theta <- 2 * pi * X[, "peak_month"] / 12
  encodings <- list(
    # Peak month on the unit circle (Dec and Jan adjacent), radius in months
    "Circular peak month" = cbind(cbind(cos(theta), sin(theta)) * 12 / (2 * pi),
                                  X[, c("peak_mult", "amp_log")]),
    "Peak month only" = X[, "peak_month", drop = FALSE],
    "Standardised features" = scale(X)
  )
  baseline <- as.integer(indices$clusters)
  tibble::tibble(
    encoding = names(encodings),
    ari_vs_main = purrr::map_dbl(encodings, \(Z) {
      adjusted_rand_index(baseline, kmeans(Z, centers = k, nstart = nstart)$cluster)
    })
  )
}

# Assign every posterior draw of each province-period seasonal curve to the
# nearest pattern centroid of the main (posterior-median) clustering
seasonality_draw_membership <- function(fit, indices, H, P,
                                        month_grid = seq(1, 12, length.out = 241)) {
  dr <- rstan::extract(fit, pars = c("a_cos", "a_sin"))
  S <- dim(dr$a_cos)[4]
  n_iter <- dim(dr$a_cos)[1]
  theta <- 2 * pi * month_grid / 12
  cos_basis <- sapply(seq_len(H), \(h) cos(h * theta))
  sin_basis <- sapply(seq_len(H), \(h) sin(h * theta))

  X <- seasonality_feature_matrix(indices)
  centres <- apply(X, 2, \(col) tapply(col, indices$clusters, mean))

  purrr::map_dfr(seq_len(P), function(p) {
    purrr::map_dfr(seq_len(S), function(s) {
      mult <- exp(dr$a_cos[, , p, s] %*% t(cos_basis) + dr$a_sin[, , p, s] %*% t(sin_basis))
      peak_mult <- apply(mult, 1, max)
      draw_features <- cbind(
        peak_month = month_grid[max.col(mult, ties.method = "first")],
        peak_mult = peak_mult,
        amp_log = log(peak_mult / apply(mult, 1, min))
      )[, colnames(centres)]
      tibble::tibble(
        region = paste0("Province ", s),
        period = paste0("Period ", p),
        .draw = seq_len(n_iter),
        draw_cluster = levels(indices$clusters)[assign_nearest_centre(draw_features, centres)]
      )
    })
  })
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

summarise_Se_over_time <- function(probs = c(0.025, 0.5, 0.975)) {
  stopifnot(exists("fit"), exists("stan_data"), exists("calendar"))
  Y <- stan_data$Y

  # Extract Se_y draws [Y]
  Se_mat <- as.matrix(fit, pars = "Se_y")
  if (is.null(colnames(Se_mat))) stop("Se_y columns have no names; re-run sampling with pars='Se_y'.")

  # Parse indices from column names [y, r]
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

  # Summarize by year and region
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
  Se_summ
}

plot_Se_over_time <- function(Se_summ, map_data) {
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
      y = "Reporting proportion",
      title = "Temporal evolution of reporting proportion"
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_bw(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          axis.title.x = element_blank()) +
    facet_wrap(.~region, ncol = 3) 

  map_se <- map_data %>%
    left_join(
      Se_summ %>% filter(year == max(Se_summ$year)),
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
    labs(caption = "Estimated reporting proportion\nby health region (latest year)") +
    theme(plot.caption = element_text(hjust = 0.5, vjust = 10, size = 12, face = "bold")) + 
    plot_annotation(theme = theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5)))

  combined_se_plot <- se_trend + map_se + patchwork::plot_layout(ncol = 2, widths  = c(2, 1)) 
  combined_se_plot + plot_annotation(theme = theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2))) 
}

################################################################################
summarise_relative_susceptibility <- function(fit, age_levels) {
  rho_mat <- as.matrix(fit, pars = "rho_pa")
  idx <- stringr::str_match(colnames(rho_mat), "rho_pa\\[(\\d+),(\\d+)\\]")
  tibble(par = colnames(rho_mat), p_idx = as.integer(idx[, 2]), a_idx = as.integer(idx[, 3])) |>
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
}

plot_relative_susceptibility <- function(rho_summ, age_levels, ref_age = "25 to 35") {
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
      y = "Relative susceptibility",
      title = paste0("Relative susceptibility\n(relative to ", ref_age, " years)")) +
    theme_bw(base_size = 14) +
    scale_y_continuous(limits = c(0.2, 1.8)) +
    coord_flip()
}
################################################################################

summarise_cfr <- function(fit, age_levels, ref_age = "25 to 35") {
  CFR_mat <- as.matrix(fit, pars = "CFR_pa")
  idx <- stringr::str_match(colnames(CFR_mat), "CFR_pa\\[(\\d+),(\\d+)\\]")
  ref_age_idx <- which(age_levels == ref_age)
  tibble(par = colnames(CFR_mat), p_idx = as.integer(idx[, 2]), a_idx = as.integer(idx[, 3])) |>
    left_join(
      as_tibble(CFR_mat) |> mutate(.draw = row_number()) |>
        pivot_longer(cols = -.draw, names_to = "par", values_to = "CFR"),
      by = "par"
    ) |>
    group_by(.draw, p_idx) |>
    mutate(ref_CFR = CFR[a_idx == ref_age_idx],
           rel_CFR = CFR / ref_CFR) |>
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
}

plot_cfr <- function(CFR_summ, age_levels, 
  period_labels = c("2003 to 2007", "2008 to 2012", "2013 to 2018", "2019 to 2024"), 
  ref_age = "25 to 35") {
  CFR_summ |>
    mutate(
      age_group = factor(age_group, levels = age_levels),
      period = factor(period, levels = paste0("Period ", 1:4), labels = period_labels)
    ) |>
    filter(age_group != ref_age) |>
    ggplot() +
    geom_abline(intercept = 1, slope = 0, linetype = "21", color = "gray30", linewidth = 1, alpha = 0.7) +
    geom_point(aes(x = age_group, y = mean_rel_CFR, color = period), size = 3, alpha = 0.8,
               position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(x = age_group, ymin = q025_ref_CFR, ymax = q975_rel_CFR, color = period), alpha = 0.8,
                  width = 0, position = position_dodge(width = 0.5), linewidth = 1.1) +
    scale_color_viridis_d(option = "magma", end = 0.8, name = "") +
    labs(x = "", y = "Relative case-fatality rate", title = paste0("Relative case-fatality rate\n(relative to ", ref_age, " years)")) +
    theme_bw(base_size = 14) +
    scale_y_continuous(breaks = seq(0, 10, by = 1)) +
    coord_flip() +
    theme(panel.grid.minor = element_blank())
}

################################################################################
ppc_cases_by_province <- function(s_idx) {
  stopifnot(exists("fit"), exists("calendar"), exists("S"), exists("T"), exists("site_levels"))
  if (!(s_idx %in% seq_len(S))) stop("s_idx out of range 1..S")

  # Extract draws from Model2_simple: C_ts (expected observed cases)
  C_mat <- as.matrix(fit, pars = "C_ts")
  phi_vec <- as.numeric(as.matrix(fit, pars = "phi_case"))
  if (is.null(colnames(C_mat))) stop("C_ts columns have no names; re-run sampling with pars='C_ts'.")

  # Parse indices from column names [t,s]
  idx <- stringr::str_match(colnames(C_mat), "C_ts\\[(\\d+),(\\d+)\\]")
  if (any(is.na(idx))) stop("Failed to parse C_ts indices from column names.")
  t_idx_col <- as.integer(idx[, 2])
  s_idx_col <- as.integer(idx[, 3])

  param_map <- tibble::tibble(par = colnames(C_mat), t_idx = t_idx_col, s = s_idx_col) |>
    dplyr::filter(s == s_idx)

  C_long <- tibble::as_tibble(C_mat[, param_map$par, drop = FALSE]) |>
    dplyr::mutate(.draw = dplyr::row_number()) |>
    tidyr::pivot_longer(cols = - .draw, names_to = "par", values_to = "C_expected") |>
    dplyr::left_join(param_map, by = "par")

  phi_df <- tibble::tibble(.draw = seq_len(nrow(C_mat)), phi = phi_vec)

  # Summarize PPC for this province over time
  pred_summ <- C_long |>
    dplyr::left_join(phi_df, by = ".draw") |>
    dplyr::group_by(t_idx) |>
    dplyr::summarize(
      mu_mean   = mean(C_expected),
      mu_median = median(C_expected),
      mu_q025   = stats::quantile(C_expected, 0.025, names = FALSE),
      mu_q975   = stats::quantile(C_expected, 0.975, names = FALSE),
      yrep_q025 = stats::quantile(stats::rnbinom(n = dplyr::n(), size = phi, mu = C_expected), 0.025, names = FALSE),
      yrep_q500 = stats::quantile(stats::rnbinom(n = dplyr::n(), size = phi, mu = C_expected), 0.50, names = FALSE),
      yrep_q975 = stats::quantile(stats::rnbinom(n = dplyr::n(), size = phi, mu = C_expected), 0.975, names = FALSE),
      .groups = "drop"
    )

  # Observed cases for this province
  obs <- calendar |>
    dplyr::filter(s == s_idx) |>
    dplyr::select(t_idx, Area, year, Month, observed_case = case)

  # Final PPC dataset for the province
  obs |>
    dplyr::left_join(pred_summ, by = "t_idx") |>
    dplyr::arrange(year, Month) |>
    dplyr::mutate(province_index = s_idx, province = site_levels[s_idx]) |>
    dplyr::mutate(date = as.Date(sprintf("%d-%02d-01", year, Month)))
}

################################################################################
# Supplementary: comprehensive per-province diagnostic figure
# Panels: A) Cases PPC  B) Inset map  C) Seasonality  D) Age dist  E) Deaths PPC

build_province_summary <- function(lambda_data, Se_summ, seasonality_indices,
                                    map_data, calendar, stan_data,
                                    period_labels = c("2003-2007", "2008-2012",
                                                      "2013-2018", "2019-2024")) {
  S <- stan_data$S
  region_id <- stan_data$region_id

  purrr::map_dfr(seq_len(S), function(s) {
    prov <- site_levels[s]
    zone <- map_data |> dplyr::filter(Provinc == prov) |> dplyr::pull(Zone)
    hr <- zone_to_hr(zone)
    r <- region_id[s]

    # PI coverage from PPC
    ppc <- ppc_cases_by_province(s)
    pi_cov <- mean(ppc$observed_case >= ppc$yrep_q025 &
                   ppc$observed_case <= ppc$yrep_q975, na.rm = TRUE)

    # Reporting rate (Se) for the province's region, most recent year
    se_prov <- Se_summ |>
      dplyr::filter(r_idx == r) |>
      dplyr::filter(year == max(year))
    se_med <- 100 * se_prov$q500
    se_lo  <- 100 * se_prov$q025
    se_hi  <- 100 * se_prov$q975

    # FOI for latest period — median across years in that period
    foi_prov <- lambda_data |>
      dplyr::filter(province == prov, year >= 2019) |>
      dplyr::summarize(q500 = median(q500), q025 = median(q025), q975 = median(q975))

    # Seasonal pattern
    prov_label <- paste0("Province ", s)
    idx_row <- seasonality_indices |>
      dplyr::filter(
        trimws(as.character(region)) == prov_label,
        trimws(as.character(period)) == paste0("Period ", stan_data$P)
      )
    pattern <- if (nrow(idx_row) > 0) as.character(idx_row$clusters) else NA_character_

    # Total deaths
    total_deaths <- calendar |> dplyr::filter(s == !!s) |> dplyr::summarize(d = sum(death)) |> dplyr::pull(d)

    # Has age data
    has_age <- any(stan_data$has_age_data_ys[, s] == 1)

    tibble::tibble(
      s_idx = s, province = prov, health_region = hr, r_idx = r, zone = zone,
      pi_coverage = round(pi_cov, 3),
      se_median = round(se_med, 1), se_low = round(se_lo, 1), se_high = round(se_hi, 1),
      foi_median = signif(foi_prov$q500, 3),
      foi_low = signif(foi_prov$q025, 3),
      foi_high = signif(foi_prov$q975, 3),
      seasonal_pattern = pattern,
      total_deaths = total_deaths,
      has_age_data = has_age
    )
  })
}

################################################################################
ppc_deaths_by_province <- function(s_idx) {
  stopifnot(exists("fit"), exists("calendar"), exists("S"), exists("T"), exists("site_levels"))
  if (!(s_idx %in% seq_len(S))) stop("s_idx out of range 1..S")

  D_mat <- as.matrix(fit, pars = "D_ts")
  phi_vec <- as.numeric(as.matrix(fit, pars = "phi_death"))

  idx <- stringr::str_match(colnames(D_mat), "D_ts\\[(\\d+),(\\d+)\\]")
  t_idx_col <- as.integer(idx[, 2])
  s_idx_col <- as.integer(idx[, 3])

  param_map <- tibble::tibble(par = colnames(D_mat), t_idx = t_idx_col, s = s_idx_col) |>
    dplyr::filter(s == s_idx)

  D_long <- tibble::as_tibble(D_mat[, param_map$par, drop = FALSE]) |>
    dplyr::mutate(.draw = dplyr::row_number()) |>
    tidyr::pivot_longer(cols = -.draw, names_to = "par", values_to = "D_expected") |>
    dplyr::left_join(param_map, by = "par")

  phi_df <- tibble::tibble(.draw = seq_len(nrow(D_mat)), phi = phi_vec)

  pred_summ <- D_long |>
    dplyr::left_join(phi_df, by = ".draw") |>
    dplyr::group_by(t_idx) |>
    dplyr::summarize(
      yrep_q025 = stats::quantile(stats::rnbinom(n = dplyr::n(), size = phi, mu = D_expected), 0.025, names = FALSE),
      yrep_q500 = stats::quantile(stats::rnbinom(n = dplyr::n(), size = phi, mu = D_expected), 0.50, names = FALSE),
      yrep_q975 = stats::quantile(stats::rnbinom(n = dplyr::n(), size = phi, mu = D_expected), 0.975, names = FALSE),
      .groups = "drop"
    )

  obs <- calendar |>
    dplyr::filter(s == s_idx) |>
    dplyr::select(t_idx, year, Month, observed_death = death)

  obs |>
    dplyr::left_join(pred_summ, by = "t_idx") |>
    dplyr::arrange(year, Month) |>
    dplyr::mutate(province = site_levels[s_idx],
                  date = as.Date(sprintf("%d-%02d-01", year, Month)))
}

################################################################################
extract_age_proportions <- function(s_idx, age_levels,
                                     period_labels = c("2003-2007", "2008-2012",
                                                       "2013-2018", "2019-2024")) {
  stopifnot(exists("fit"), exists("stan_data"))
  Y <- stan_data$Y; A <- stan_data$A; P <- stan_data$P

  # Observed age-stratified cases [Y, A, S]
  obs_yas <- stan_data$R_case_yas[, , s_idx]   # Y x A matrix
  has_age <- stan_data$has_age_data_ys[, s_idx]  # Y-length vector

  # Predicted proportions from rho_pa × N_tas (averaged over months within year)
  rho_mat <- as.matrix(fit, pars = "rho_pa")
  rho_median <- matrix(apply(rho_mat, 2, median), nrow = P, ncol = A, byrow = FALSE)

  N_tas <- stan_data$N_tas  # [T, A, S]
  period_id <- stan_data$period_id  # Y-length
  Tmonths <- stan_data$T

  # Map month index → year index (12 months per year)
  month_to_year <- rep(seq_len(Y), each = 12)[seq_len(Tmonths)]

  # Average population per year-age-site
  pred_list <- list()
  for (y in seq_len(Y)) {
    p <- period_id[y]
    t_months <- which(month_to_year == y)
    # Mean pop across months in this year
    N_as <- colMeans(matrix(N_tas[t_months, , s_idx], ncol = A))
    expected_a <- rho_median[p, ] * N_as
    total <- sum(expected_a)
    pred_list[[y]] <- tibble::tibble(
      y_idx = y, a_idx = seq_len(A),
      pred_prop = if (total > 0) expected_a / total else rep(1 / A, A)
    )
  }
  pred_df <- dplyr::bind_rows(pred_list)

  # Observed proportions (aggregate by period)
  obs_long <- tidyr::crossing(y_idx = seq_len(Y), a_idx = seq_len(A)) |>
    dplyr::mutate(
      obs_cases = purrr::map2_dbl(y_idx, a_idx, ~ obs_yas[.x, .y]),
      has_data = has_age[y_idx],
      p_idx = period_id[y_idx]
    ) |>
    dplyr::filter(has_data == 1)

  combined <- obs_long |>
    dplyr::left_join(pred_df, by = c("y_idx", "a_idx")) |>
    dplyr::group_by(p_idx, a_idx) |>
    dplyr::summarize(obs = sum(obs_cases), pred_prop = mean(pred_prop), .groups = "drop") |>
    dplyr::group_by(p_idx) |>
    dplyr::mutate(obs_prop = obs / sum(obs)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      age_group = factor(age_levels[a_idx], levels = age_levels),
      period = factor(period_labels[p_idx], levels = period_labels)
    )

  combined
}

################################################################################
plot_province_diagnostic <- function(s_idx,
                                      province_summary,
                                      season_curves,
                                      map_data,
                                      config = list(
                                        use_log = FALSE,
                                        show_age = TRUE
                                      ),
                                      period_labels = c("2003-2007", "2008-2012",
                                                        "2013-2018", "2019-2024"),
                                      period_colors = c("#1b9e77", "#d95f02",
                                                        "#7570b3", "#e7298a"),
                                      base_size = 11,
                                      fig_number = NULL,
                                      output_dir = NULL) {

  prov_info <- province_summary |> dplyr::filter(s_idx == !!s_idx)
  province_name <- prov_info$province
  province_entity <- paste0("Province ", s_idx, " ")

  # ---------- common theme ----------
  base_theme <- theme_minimal(base_size = base_size) +
    theme(
      axis.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = base_size)
    )

  # ===== Panel A: Cases PPC =====
  ppc_data <- ppc_cases_by_province(s_idx)
  period_boundaries <- as.Date(c("2007-12-15", "2012-12-15", "2018-12-15"))

  p_ppc <- ggplot(ppc_data, aes(x = date)) +
    geom_vline(xintercept = period_boundaries, linetype = "dashed", color = "grey60", linewidth = 0.3) +
    geom_ribbon(aes(ymin = yrep_q025, ymax = yrep_q975), fill = "steelblue", alpha = 0.18) +
    geom_line(aes(y = yrep_q500), color = "steelblue", linewidth = 0.7) +
    geom_point(aes(y = observed_case), color = "black", size = 0.9, alpha = 0.6) +
    scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
    labs(x = NULL, y = "Reported cases", title = "A. Posterior predictive check (cases)") +
    base_theme

  if (isTRUE(config$use_log)) {
    p_ppc <- p_ppc + scale_y_log10(labels = scales::label_comma())
  }

  # ===== Panel B: Inset map =====
  zone_val <- prov_info$zone
  prov_fill <- "#e31a1c"

  map_highlight <- map_data |>
    dplyr::mutate(
      fill_cat = dplyr::case_when(
        Provinc == province_name ~ "Province",
        Zone == zone_val ~ "Health Region",
        TRUE ~ "Other"
      ),
      fill_cat = factor(fill_cat, levels = c("Province", "Health Region", "Other"))
    )

  p_map <- ggplot(map_highlight) +
    geom_sf(aes(fill = fill_cat), color = "white", linewidth = 0.15) +
    scale_fill_manual(values = c("Province" = prov_fill, "Health Region" = "#a6cee3", "Other" = "grey92"),
                      name = NULL) +
    coord_sf(xlim = c(96, 106), ylim = c(5, 21), expand = FALSE) +
    theme_void(base_size = base_size) +
    labs(title = paste0("B. ", prov_info$health_region)) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size)
    )

  # ===== Panel D: Seasonality =====
  entity_label <- trimws(province_entity)
  season_prov <- season_curves |>
    dplyr::filter(trimws(as.character(entity)) == entity_label)

  n_per <- length(unique(season_prov$period))
  legend_labels <- setNames(period_labels[seq_len(n_per)], levels(season_prov$period)[seq_len(n_per)])
  pal <- setNames(period_colors[seq_len(n_per)], levels(season_prov$period)[seq_len(n_per)])

  # Look up seasonal pattern for each period
  prov_label <- paste0("Province ", s_idx)
  pattern_by_period <- seasonality_indices |>
    dplyr::filter(trimws(as.character(region)) == prov_label) |>
    dplyr::mutate(p_idx = as.integer(gsub("Period ", "", trimws(as.character(period))))) |>
    dplyr::arrange(p_idx)

  p_season <- ggplot(season_prov, aes(x = x)) +
    geom_ribbon(aes(ymin = q025, ymax = q975, fill = period), alpha = 0.15) +
    geom_line(aes(y = q500, color = period), linewidth = 0.7) +
    geom_hline(yintercept = 1, linetype = "31", alpha = 0.6) +
    scale_x_continuous(breaks = seq(1, 12, by = 2), labels = month.abb[seq(1, 12, by = 2)],
                       minor_breaks = NULL) +
    scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, 0.5)) +
    scale_color_manual(values = pal, labels = legend_labels, name = NULL) +
    scale_fill_manual(values = pal, labels = legend_labels, name = NULL) +
    labs(x = NULL, y = "Seasonal multiplier", title = "D. Seasonality by period") +
    base_theme +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.4, "cm"),
          legend.text = element_text(size = base_size - 2)) +
    guides(color = guide_legend(nrow = 2), fill = guide_legend(nrow = 2))

  # Annotate seasonal pattern as a single text block (top-left)
  if (nrow(pattern_by_period) > 0) {
    annot_lines <- paste0(
      period_labels[pattern_by_period$p_idx], ": ",
      as.character(pattern_by_period$clusters)
    )
    annot_text <- paste(annot_lines, collapse = "\n")
    p_season <- p_season +
      annotate("text", x = 1, y = 4, label = annot_text,
               hjust = 0, vjust = 1, size = 2.8, fontface = "bold", color = "black")
  }

  # ===== Panel C: Age distribution (conditional) =====
  show_age <- isTRUE(config$show_age) && prov_info$has_age_data
  if (show_age) {
    age_props <- extract_age_proportions(s_idx, age_levels, period_labels)
    age_long <- age_props |>
      tidyr::pivot_longer(cols = c(obs_prop, pred_prop), names_to = "source", values_to = "proportion") |>
      dplyr::mutate(source = ifelse(source == "obs_prop", "Observed", "Predicted"))

    p_age <- ggplot(age_long, aes(x = age_group, y = proportion, fill = source)) +
      geom_col(position = "dodge", alpha = 0.8, width = 0.7) +
      facet_wrap(~ period, nrow = 1) +
      scale_fill_manual(values = c("Observed" = "#213448", "Predicted" = "#7fcdbb"), name = NULL) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(x = NULL, y = "Proportion", title = "C. Age distribution: observed vs predicted") +
      base_theme +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = base_size - 3),
            legend.position = "bottom",
            legend.key.size = unit(0.35, "cm"))
  }

  # ===== Layout composition =====
  row1 <- p_ppc | p_map
  row1 <- row1 + plot_layout(widths = c(3, 1))

  if (show_age) {
    row2 <- p_age | p_season
    row2 <- row2 + plot_layout(widths = c(3, 1))
    final <- row1 / row2 + plot_layout(heights = c(1, 1))
    fig_height <- 6
  } else {
    final <- row1 / p_season + plot_layout(heights = c(1, 1))
    fig_height <- 6
  }

  final <- final +
    plot_annotation(
      title = province_name,
      caption = "Data: Thailand MoPH Bureau of Epidemiology (2003\u20132024). Model: Bayesian latent process model (see Methods).",
      theme = theme(
        plot.title = element_text(face = "bold", size = base_size + 5, hjust = 0.5),
        plot.caption = element_text(size = base_size - 3, hjust = 0.5, color = "grey50",
                                    margin = margin(t = 8)),
        plot.margin = margin(8, 8, 8, 8)
      )
    )

  # ===== Save =====
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    safe_name <- gsub("[^A-Za-z0-9]", "_", province_name)
    fig_num <- if (is.null(fig_number)) s_idx else fig_number
    outfile <- file.path(output_dir, sprintf("SuppFig_%02d_%s.pdf", fig_num, safe_name))
    ggplot2::ggsave(outfile, final, width = 10, height = fig_height, units = "in", dpi = 300)
  }

  final
}
