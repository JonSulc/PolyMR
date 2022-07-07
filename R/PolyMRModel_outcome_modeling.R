update_outcome_model.PolyMRModel <- function(polymr_model,
                                             exposure_powers,
                                             control_function_powers = NULL){
  if (any(exposure_powers > polymr_model$max_exposure_power))
    warning(sprintf("Modeling exposure powers above the defined max_exposure_power
max(exposure_powers) = %i, max_exposure_power = %i",
                    max(exposure_powers), polymr_model$max_exposure_power))

  polymr_model$exposure_powers <- exposure_powers

  polymr_model <- update_control_function_powers(
      polymr_model,
      control_function_powers = control_function_powers
    )
  model_outcome(polymr_model)
}


update_control_function_powers <- function(polymr_model,
                                           control_function_powers = NULL,
                                           update_outcome_model = FALSE){
  old_control_function_powers <- polymr_model$control_function_powers
  if (!is.null(control_function_powers))
    polymr_model$control_function_powers <- control_function_powers
  else
    polymr_model$control_function_powers <-
      1:min(max(polymr_model$exposure_powers, 1),
            polymr_model$max_control_function_power)

  if (update_outcome_model &&
      !setequal(old_control_function_powers,
                polymr_model$control_function_powers))
    polymr_model <- model_outcome(polymr_model)

  polymr_model
}



