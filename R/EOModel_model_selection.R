increase_order_of_outcome_model.EOModel <- function(eo_model,
                                                    power_step = 2){
  new_exposure_powers <- c(eo_model$exposure_powers,
                           max(eo_model$exposure_powers) + 1:power_step)
  new_exposure_powers <-
    new_exposure_powers[new_exposure_powers <= eo_model$max_exposure_power]

  update_outcome_model(
    eo_model,
    new_exposure_powers
  )
}


drop_least_significant_term.EOModel <- function(eo_model, ...){
  if (length(eo_model$exposure_powers) == 1){
    warning("No significant terms left, returning null model.")
    return(update_outcome_model(eo_model, NULL))
  }

  weakest_term <- which.max(
    get_outcome_model_coefficients_exposure_pvalues(eo_model)
    )

  update_outcome_model(
    eo_model,
    eo_model$exposure_powers[-weakest_term]
  )
}


are_terms_significant.EOModel <- function(
    eo_model,
    p_thr = NULL,
    only_new = length(eo_model$exposure_powers)
){

  if (is.null(p_thr))
    p_thr <- .05 / only_new

  p_values <- tail(
    get_outcome_model_coefficients_exposure_pvalues(eo_model),
    only_new
    )
  p_values <= p_thr
}
