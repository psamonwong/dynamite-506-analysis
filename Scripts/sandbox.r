

make_bivar_map <- function(lambda_data,
                                              map_data,
                                              year0 = 2003,
                                              year1 = 2024,
                                              foi_breaks = c(-Inf, 100, 1000, Inf),
                                              foi_labels = c("< 100", "100-1000", "> 1000"),
                                              pal = "DkViolet2",
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
    geom_sf(aes(fill = bi_class), color = "white", linewidth = 0.15, show.legend = FALSE) +
    biscale::bi_scale_fill(pal = pal, dim = length(foi_labels), na.value = "grey90") +
    coord_sf(xlim = c(96, 106), ylim = c(5, 21), expand = FALSE) +
    theme_void(base_size = 14) +
    labs(title = subtitle) +
    theme(
      plot.title = element_text(face = "bold", size = 12, margin = margin(b = 6)),
      plot.margin = margin(6, 6, 6, 6)
    )

  legend_plot <- biscale::bi_legend(
    pal = pal,
    dim = length(foi_labels),
    xlab = year0 %>% as.character(),
    ylab = year1 %>% as.character(),
    size = 15
  ) +
    scale_x_continuous(breaks = seq_along(foi_labels), labels = foi_labels, expand = c(0, 0)) +
    scale_y_continuous(breaks = seq_along(foi_labels), labels = foi_labels, expand = c(0, 0)) +
    theme(
      axis.text = element_text(size = 6),
      axis.title = element_text(size = 7, face = "bold"),
      plot.margin = margin(0, 0, 0, 0)
    )

  p_map +
    annotation_custom(
      grob = ggplotGrob(legend_plot),
      xmin = legend_bbox["xmin"], xmax = legend_bbox["xmax"],
      ymin = legend_bbox["ymin"], ymax = legend_bbox["ymax"]
    )
}

p_bivar <- make_bivar_map(
            lambda_data, 
            map_data, 
            year0 = 2003, year1 = 2024, 
            pal = "DkBlue", 
            subtitle = "Bivariate choropleth: 2003 vs 2024 (median FOI)")
p_bivar
