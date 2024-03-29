#' @exportS3Method
predict.PolyMRModel <- function(polymr_model, exposure) {
  cf_columns <- grep("^cf[0-9]+$", names(coef(polymr_model)), value = TRUE)
  outcome_predictors <- array(0, dim = c(length(exposure),
                                         length(cf_columns)))|>
    as.data.frame()
  colnames(outcome_predictors) <- cf_columns

  exposure_powers <- get_exposure_powers(polymr_model)
  outcome_predictors <- create_power_table(exposure,
                                           exposure_powers,
                                           "exposure") |>
    cbind(outcome_predictors)
  predict(polymr_model$outcome_model, outcome_predictors)
}