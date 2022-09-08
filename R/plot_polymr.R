#' Plot results of PolyMR
#'
#' This function provides an easy interface to quickly visualize the results
#' from PolyMR using \code{ggplot2}.
#'
#' @param polymr_results A list containing the results from the \code{polymr}
#'   function.
#'
#' @param show_polymr Boolean indicating whether or not to plot the inferred
#'   causal function (default is TRUE).
#' @param show_observational Boolean indicating whether or not to plot the
#'   observed association. Requires the observational association to have been
#'   estimated by the \code{polymr} function using
#'   \code{return_observational_function = TRUE}. Default is TRUE.
#' @param show_binned_observations Boolean indicating whether or not to plot the
#'   binned observations, specifically the median exposure and outcome for each
#'   bin. Requires the binned values to have been returned by the \code{polymr}
#'   function using \code{return_binned_observations = TRUE}. Default is TRUE.
#' @param show_confidence_ribbons Boolean indicating whether or not to display
#'   the 95\% confidence hulls (default is TRUE).
#' @param scale_values Boolean indicating whether or not to leave phenotypes
#'   (exposure and outcome) standardized with mean 0 and variance 1 or convert
#'   them back to their original range using the values in
#'   \code{phenotypes_summary}. Default is TRUE, indicating that the values
#'   will be plotted on a standard scale.
#' @param output_filename The file name of the plot to be saved, if any. This
#'   will override any \code{filename} provided in \code{ggsave_options}.
#'   Default is \code{NULL}, in which case the plot will not be saved.
#' @param main_title Main title of the plot (default is "Effect of exposure on
#'   outcome").
#' @param subtitle Subtitle of the plot. If \code{auto_subtitle} is enabled,
#'   this will be prepended to the automatic subtitle. Default is \code{NULL}.
#' @param auto_subtitle Whether to create and include a summary of results as
#'   subtitle. This summary includes the p-values and variance explained for the
#'   plotted functions (total variance for the observational association,
#'   causal variance for PolyMR). Default is \code{TRUE}.
#' @param xlab Label for the X axis. Default is "Exposure (scaled)" if
#'   \code{scale_values} is \code{TRUE} and "Exposure" otherwise.
#' @param ylab Label for the Y axis. Default is "Outcome (scaled)" if
#'   \code{scale_values} is \code{TRUE} and "Outcome" otherwise.
#' @param ylims An optional vector of length 2 defining the boundaries of the
#'   plotted region on the Y axis. Default is NULL, relying on \code{ggplot2}
#'   standard behavior.
#' @param show_legend Whether or not to display the legend next to the plot
#'   (default is \code{TRUE}).
#' @param ggsave_options A list of options to provide to \code{ggplot2::ggsave}.
#'   Ignored if no filename is specified here or using \code{output_filename}.
#'   Default is a list specifying a \code{width} and \code{height} both equal to
#'   18 \code{cm}.
#' @param n_points Number of points to plot in rendering the curve. Default is
#'   1000, resulting in a smooth enough curve for most purposes.
#'
#' @return An object of the \code{gg} and \code{ggplot} classes.
#'
#' @details The plot is created using the \code{gpplot2} library. The resulting
#'   plot can be customized using the standard \code{ggplot2} functions.
#'
#' @export
plot_polymr <- function(
    polymr_results,
    show_polymr = TRUE,
    show_observational = TRUE,
    show_binned_observations = TRUE,
    show_confidence_ribbons = TRUE,
    scale_values = TRUE,
    output_filename = NULL,
    main_title = "Effects of exposure on outcome",
    subtitle = NULL,
    auto_subtitle = TRUE,
    xlab = ifelse(scale_values,
                  "Exposure (scaled)",
                  "Exposure"),
    ylab = ifelse(scale_values,
                  "Outcome (scaled)",
                  "Outcome"),
    ylims = NULL,
    show_legend = TRUE,
    ggsave_options = list(width = 18,
                          height = 18,
                          units = "cm"),
    n_points = 1000) {

  exposure_range <- range(polymr_results$binned_observations$exposure_median) |>
    scale(center = polymr_results$phenotypes_summary$mean[1],
          scale  = polymr_results$phenotypes_summary$sd[1])

  exposure_values <- seq(exposure_range[1],
                         exposure_range[2],
                         length.out = n_points)

  p <- base_plot()

  if (show_binned_observations)
    p <- p + plot_binned_observations(polymr_results,
                                      scale_values = scale_values)

  if (show_observational) {
    p <- p + plot_curve(polymr_results$observational,
                        exposure_values,
                        method_name = "Observational",
                        show_confidence_ribbons = show_confidence_ribbons,
                        unscale_values = !scale_values,
                        phenotypes_summary = polymr_results$phenotypes_summary)
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
                        show_confidence_ribbons = show_confidence_ribbons,
                        unscale_values = !scale_values,
                        phenotypes_summary = polymr_results$phenotypes_summary)
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

  if (!show_legend)
    p <- p + ggplot2::theme(legend.position = "none")

  if (!is.null(output_filename))
    ggsave_options$filename <- output_filename

  if (!is.null(ggsave_options$filename))
    do.call(ggplot2::ggsave, ggsave_options)

  p
}

get_subtitle <- function(eo_model,
                         method_name) {
  if (is.null(eo_model$pval_linear_model))
    return()

  sprintf("%s non-linearity p-value: %.2e, %sariance explained: %.2e",
            method_name,
            eo_model$pval_linear_model,
            ifelse(method_name == "Observational", "V", "Causal v"),
            eo_model$r_squared)
}

plot_binned_observations <- function(polymr_results,
                                     scale_values = TRUE,
                                     method_name = "Observational") {
  observations <- data.frame(
    exposure = polymr_results$binned_observations$exposure_median,
    outcome  = polymr_results$binned_observations$outcome_median,
    method = method_name
  )

  if (scale_values) {
    observations$exposure <- observations$exposure |>
      scale(center = polymr_results$phenotypes_summary$mean[1],
            scale  = polymr_results$phenotypes_summary$sd[1])
    observations$outcome <- observations$outcome |>
      scale(center = polymr_results$phenotypes_summary$mean[2],
            scale  = polymr_results$phenotypes_summary$sd[2])
  }
  ggplot2::geom_point(data = observations)
}
