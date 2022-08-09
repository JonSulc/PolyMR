get_exposure_powers <- function(eo_model) {
  eo_coefficients <- coef(eo_model)
  grep("[0-9]+",
       grep("^exposure([0-9]+)$", names(eo_coefficients))) |>
    as.numeric()
}



#' @exportS3Method
predict.EOModel <- function(eo_model, exposure) {
  exposure_powers <- get_exposure_powers(eo_model)
  outcome_predictors <- create_power_table(exposure,
                                           exposure_powers,
                                           "exposure")
  predict(eo_model$outcome_model, outcome_predictors)
}

get_confidence_hull <- function(eo_model,
                                exposure,
                                to_keep = "(Intercept)|(exposure)" |>
                                  grepl(colnames(eo_model$vcov))) {
  vcov <- eo_model$vcov[to_keep, to_keep]
}