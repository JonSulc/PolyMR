get_outcome_model_coefficients.PolyMRModel <- function(
    polymr_model,
    include_exposure = TRUE,
    include_control_function = TRUE,
    include_intercept = FALSE,
    values_to_return = c("Estimate", "Pr(>|t|)")
){
  get_outcome_model_coefficients.EOModel(
    polymr_model,
    include_exposure = include_exposure,
    include_intercept = include_intercept,
    terms_to_include = if (include_control_function)
      paste0("ex", polymr_model$control_function_powers),
    values_to_return = values_to_return
  )
}

get_outcome_model_coefficients_exposure.PolyMRModel <- function(polymr_model,
                                                                ...){
  get_outcome_model_coefficients(polymr_model,
                                 include_exposure = TRUE,
                                 include_control_function = FALSE,
                                 ...)
}

get_outcome_model_coefficients_exposure_pvalues.PolyMRModel <- function(polymr_model){
  get_outcome_model_coefficients_exposure(polymr_model,
                                          values = "Pr(>|t|)")
}


get_null_model <- function(polymr_model){
  control_function_terms <- paste0("ex", polymr_model$control_function_powers) |>
    paste(collapse = " + ")
  null_model_formula <- as.formula(paste(". ~", control_function_terms))
  update(polymr_model$outcome_model,
         null_model_formula,
         data = polymr_model$outcome_predictors)
}

get_pval_vs_null_model.PolyMRModel <- function(polymr_model){
  null_model <- get_null_model(polymr_model)
  get_lrt_pval_alternate_model(polymr_model, null_model)
}

get_r_squared.PolyMRModel <- function(polymr_model){
  null_model <- get_null_model(polymr_model)
  summary(polymr_model$outcome_model)$r.squared - summary(null_model)$r.squared
}
