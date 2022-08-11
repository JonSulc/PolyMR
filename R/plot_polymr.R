#' @export
plot_polymr <- function(
    polymr_results,
    show_polymr = TRUE,
    show_observational = TRUE,
    show_binned_observations = TRUE,
    show_confidence_ribbons = TRUE,
    scaled = TRUE,
    output_filename = NULL,
    n_points = 1000,
    main_title = "Effects of exposure on outcome",
    subtitle = NULL,
    auto_subtitle = TRUE,
    xlab = "Exposure (scaled)",
    ylab = "Outcome (scaled)",
    ylims = NULL,
    show_legend = TRUE,
    scale_color_brewer_arguments = NULL,
    ggsave_options = list(width = 18,
                          height = 18,
                          units = "cm",
                          filename = output_filename)) {

  if (!scaled)
    stop("Not yet implemented")

  exposure_values <-
    seq(from = min(polymr_results$binned_observations$exposure_median -
                     polymr_results$phenotypes_summary$mean[1]) /
                 polymr_results$phenotypes_summary$sd[1],
        to   = max(polymr_results$binned_observations$exposure_median -
                     polymr_results$phenotypes_summary$mean[1]) /
                 polymr_results$phenotypes_summary$sd[1],
            length.out = n_points)

  p <- base_plot()

  if (show_binned_observations)
    p <- p + plot_binned_observations(polymr_results)

  if (show_observational) {
    p <- p + plot_curve(polymr_results$observational,
                        exposure_values,
                        method_name = "Observational",
                        show_confidence_ribbons = show_confidence_ribbons)
    if (auto_subtitle)
      subtitle <- c(subtitle,
                    get_subtitle(polymr_results$observational,
                                 method_name = "Observational")) |>
        paste(collapse = "\n")
  }

  if (show_polymr) {
    p <- p + plot_curve(polymr_results$polymr,
                        exposure_values,
                        method_name = "PolyMR",
                        show_confidence_ribbons = show_confidence_ribbons)
    if (auto_subtitle)
      subtitle <- c(subtitle,
                    get_subtitle(polymr_results$polymr,
                                 method_name = "PolyMR")) |>
        paste(collapse = "\n")
  }

  p <- p +
    ggplot2::labs(title = main_title,
                  subtitle = subtitle,
                  x = xlab,
                  y = ylab,
                  color = "Method",
                  linetype = "Method") +
    ggplot2::theme_classic()

  if (!is.null(ylims))
    p <- p + ggplot2::coord_cartesian(ylim = ylims)

  if (!is.null(scale_color_brewer_arguments))
    p <- p + do.call(ggplot2::scale_color_brewer, scale_color_brewer_arguments)

  if (!show_legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  if (!is.null(ggsave_options$filename)) {
    do.call(ggplot2::ggsave, ggsave_options)
  }

  p
}

get_subtitle <- function(eo_model,
                         method_name) {
  if (is.null(eo_model$pval_linear_model))
    return()

  sprintf("%s non-linearity p-value: %.2e, %sariance: %.2e",
            method_name,
            eo_model$pval_linear_model,
            ifelse(method_name == "observational", "V", "Causal v"),
            eo_model$r_squared)
}

plot_binned_observations <- function(polymr_results) {
  observations <- data.frame(
    exposure = (polymr_results$binned_observations$exposure_median -
      polymr_results$phenotypes_summary$mean[1]) /
      polymr_results$phenotypes_summary$sd[1],
    outcome = (polymr_results$binned_observations$outcome_median -
      polymr_results$phenotypes_summary$mean[2]) /
      polymr_results$phenotypes_summary$sd[2],
    method = "Observational")
  ggplot2::geom_point(data = observations)
}