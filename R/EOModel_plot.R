#' @import ggplot2
#' @importFrom MASS mvrnorm

get_exposure_powers <- function(eo_model) {
  eo_coefficients <- coef(eo_model)
  is_exposure <- grepl("^exposure[0-9]+$",
                       names(eo_coefficients))
  exposure_terms <- eo_coefficients[is_exposure] |>
    names()
  regmatches(exposure_terms,
             regexpr("[0-9]+", exposure_terms)) |>
    as.numeric()
}

#' @exportS3Method
plot.EOModel <- function(eo_model,
                         xlim = NULL,
                         n_points = 1000,
                         add = FALSE,
                         ...,
                         show_confidence_ribbons = TRUE) {
  if (is.null(xlim)) {
    warning("It is highly recommended to specify the exposure range (xlim).")
    xlim <- qnorm(c(1, n_points) / (n_points + 1))
  }

  exposure_values <- seq(xlim[1],
                         xlim[2],
                         length.out = n_points)

  eo_model_plot <- plot_curve(eo_model,
                              exposure_values,
                              ...,
                              show_confidence_ribbons = show_confidence_ribbons)
  if (add)
    return(eo_model_plot)

  base_plot() +
    eo_model_plot
}

#' @exportS3Method
predict.EOModel <- function(eo_model, exposure_values) {
  exposure_powers <- get_exposure_powers(eo_model)
  outcome_predictors <- create_power_table(exposure_values,
                                           exposure_powers,
                                           "exposure")
  predict(eo_model$outcome_model, outcome_predictors)
}

base_plot <- function() {
  data.frame(exposure = double(),
             outcome = double(),
             method = character()) |>
  ggplot2::ggplot(ggplot2::aes(x = exposure,
                               y = outcome,
                               group = method,
                               color = method))
}

plot_curve <- function(eo_model,
                       exposure_values,
                       method_name = "",
                       show_confidence_ribbons = TRUE,
                       unscale_values = FALSE,
                       phenotypes_summary = NULL,
                       ...) {
  if (is.null(coef(eo_model))) {
    to_plot <- data.frame(exposure = range(exposure_values),
                          outcome = c(0,0),
                          method = method_name)
  } else {
    to_plot <- data.frame(exposure = exposure_values,
                          outcome = predict(eo_model, exposure_values),
                          method = method_name)
  }
  if (unscale_values) {
    to_plot$exposure <- to_plot$exposure |>
      unscale(phenotypes_summary$mean[1],
              phenotypes_summary$sd[1])
    to_plot$outcome <- to_plot$outcome |>
      unscale(phenotypes_summary$mean[2],
              phenotypes_summary$sd[2])
  }
  list(ggplot2::geom_line(data = to_plot, ggplot2::aes(linetype = method)),
       if (show_confidence_ribbons)
         plot_confidence_ribbon(eo_model,
                                exposure_values,
                                method_name,
                                unscale_values = unscale_values,
                                phenotypes_summary = phenotypes_summary))
}

get_confidence_boundaries <- function(
    eo_model,
    exposure_values,
    n_sim = 1000,
    interval = c(0.025, 0.975),
    to_keep = "(Intercept)|(exposure)" |>
      grep(colnames(eo_model$vcov), value = TRUE),
    unscale_values = FALSE,
    phenotypes_summary = NULL) {

  exposure_powers <- get_exposure_powers(eo_model)

  if ("(Intercept)" %in% to_keep)
    exposure_powers <- c(0, exposure_powers)

  exposure_table <- outer(exposure_values, exposure_powers, "^") |>
    t()

  random_curves <- MASS::mvrnorm(n = n_sim,
                                 mu = coef(eo_model)[to_keep],
                                 Sigma = eo_model$vcov[to_keep, to_keep],
                                 empirical = TRUE) %*%
    exposure_table

  boundaries <- apply(random_curves, 2, quantile, interval) |>
    t() |>
    as.data.frame() |>
    setNames(c("lower", "upper"))
  boundaries$exposure <- exposure_values

  if (!unscale_values)
    return(boundaries)

  boundaries$exposure <- boundaries$exposure |>
    unscale(phenotypes_summary$mean[1],
            phenotypes_summary$sd[1])
  boundaries[c("lower", "upper")] <- boundaries[c("lower", "upper")] |>
    unscale(phenotypes_summary$mean[2],
            phenotypes_summary$sd[2])
  boundaries
}

plot_confidence_ribbon <- function(eo_model,
                                   exposure_values,
                                   method_name = "",
                                   geom_ribbon_arguments = list(alpha = .1),
                                   ...) {
  ribbon_boundaries <- get_confidence_boundaries(eo_model,
                                                 exposure_values,
                                                 ...)
  ribbon_boundaries$outcome <- 0
  ribbon_boundaries$method <- method_name

  do.call(ggplot2::geom_ribbon,
          modifyList(list(data = ribbon_boundaries,
                          ggplot2::aes(ymin = lower,
                                       ymax = upper)),
                     geom_ribbon_arguments))
}

unscale <- function(values, center, scale) {
  values * scale + center
}
