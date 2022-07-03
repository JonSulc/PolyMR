drop_least_significant_term.PolyMRModel <- function(
    polymr_model,
    drop_higher_control_function_powers = TRUE
){
  if (length(polymr_model$exposure_powers) == 1){
    warning("No significant terms left, returning null model.")
    return(update_outcome_model(polymr_model, NULL))
  }

  weakest_term <- get_outcome_model_coefficients_exposure_pvalues(polymr_model) |>
    which.max()

  polymr_model$exposure_powers <- polymr_model$exposure_powers[-weakest_term]

  if (drop_higher_control_function_powers)
    polymr_model <- update_control_function_powers(polymr_model)

  model_outcome(polymr_model)
}
