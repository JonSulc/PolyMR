# EOModel represents the base for exposure-outcome models and is used to model
# the observed association. The PolyMRModel class extends this, redefining some
# of its methods where necessary.

new_EOModel <- function(eo_model = NULL,
                        exposure = eo_model$exposure,
                        outcome = eo_model$outcome,
                        exposure_powers = 1,
                        max_exposure_power = 10){
  if ((is.null(eo_model$exposure)  && is.null(exposure)) ||
      (is.null(eo_model$outcome)   && is.null(outcome)))
    stop("Exposure and outcome must both be provided.")

  if (length(exposure) != length(outcome))
    stop("Exposure and outcome sample sizes must match.")

  eo_model <- structure(
    list(
      exposure = exposure,
      outcome = outcome,
      exposure_powers = exposure_powers,
      max_exposure_power = 10
    ),
    class = "EOModel"
  )
  eo_model <- preprocess_data(eo_model)
  model(eo_model)
}


model <- function(eo_model, ...){
  UseMethod("model")
}
model.EOModel <- function(eo_model){
  model_outcome(eo_model)
}


return_model <- function(eo_model, ...){
  UseMethod("return_model")
}
return_model.EOModel <- function(eo_model, ...){
  eo_model <- calculate_vcov(eo_model)
    cleanup(eo_model, ...)
}

cleanup <- function(eo_model,
                    to_return = c("outcome_model",
                                  "vcov",
                                  "pval_null_model",
                                  "pval_linear_model",
                                  "r_squared"),
                            ...){
  values_to_remove <- setdiff(names(eo_model), to_return)
  eo_model[values_to_remove] <- NULL
  eo_model
}





