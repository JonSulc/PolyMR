get_exposure_powers <- function(eo_model) {
  eo_coefficients <- coef(eo_model)
  grep("[0-9]+",
       grep("^exposure([0-9]+)$", names(eo_coefficients))) |>
    as.numeric()
}

#' @exportS3Method
plot.EOModel <- function(eo_model,
                         ...,
                         show_confidence_ribbons = TRUE) {
  p <- base_plot(eo_model, ...) +
    plot_curve(eo_model, ...)

  if (show_confidence_ribbons)
    p <- p + plot_confidence_ribbon(eo_model, ...)
  p
}

#' @exportS3Method
predict.EOModel <- function(eo_model, exposure) {
  exposure_powers <- get_exposure_powers(eo_model)
  outcome_predictors <- create_power_table(exposure,
                                           exposure_powers,
                                           "exposure")
  predict(eo_model$outcome_model, outcome_predictors)
}

base_plot <- function(eo_model,
                      exposure,
                      method_name) {
  to_plot <- data.frame(x = double(),
                        estimate = double(),
                        method = character())
  ggplot2::ggplot(to_plot,
                  ggplot2::aes(x = x,
                               y = estimate,
                               group = method,
                               color = method))
}

plot_curve <- function(eo_model,
                       exposure,
                       method_name) {
  if (is.null(coef(eo_model))) {
    to_plot <- data.frame(x = range(exposure),
                          method = method_name,
                          estimate = c(0, 0))
  } else {
    to_plot <- data.frame(x = exposure,
                          method = method_name,
                          estimate = predict(eo_model, exposure))
  }
  ggplot2::geom_line(data = to_plot,
                     ggplot2::aes(linetype = method))
}

get_confidence_boundaries <- function(
    eo_model,
    exposure_values,
    n_sim = 1000,
    interval = c(0.025, 0.975),
    to_keep = "(Intercept)|(exposure)" |>
      grep(colnames(eo_model$vcov), value = TRUE)) {

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

  apply(random_curves, 2, quantile, interval) |>
    t() |>
    as.data.frame() |>
    setNames(c("lower", "upper"))
}

plot_confidence_ribbon <- function(eo_model,
                                   exposure_values,
                                   method_name = "Observational",
                                   geom_ribbon_arguments = list(alpha = .1),
                                   ...) {
  ribbon_boundaries <- get_confidence_boundaries(eo_model,
                                                 exposure_values,
                                                 ...)
  ribbon_boundaries$method <- method_name
  ribbon_boundaries$x <- exposure_values
  ribbon_boundaries$estimate <- 0

  do.call(ggplot2::geom_ribbon,
          modifyList(list(data = ribbon_boundaries,
                          ggplot2::aes(ymin = lower,
                                       ymax = upper)),
                     geom_ribbon_arguments))
}