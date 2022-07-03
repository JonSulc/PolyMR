model_exposure <- function(polymr_model){
  calculate_control_function(polymr_model)
}


calculate_control_function <- function(polymr_model,
                                       recompute = FALSE){
  if (!is.null(polymr_model$control_function) && !recompute)
    return(polymr_model)
  polymr_model <- calculate_beta_exposure(polymr_model)
  polymr_model$control_function <-
    c(polymr_model$exposure -
        polymr_model$genotypes %*% polymr_model$beta_exposure)
  polymr_model
}


calculate_beta_exposure <- function(polymr_model,
                                    recompute = FALSE){
  if (!is.null(polymr_model$beta_exposure) &&
      !recompute &&
      (is.null(polymr_model$reverse_t_thr) || !is.null(polymr_model$se_exposure)))
    return(polymr_model)

  full_exposure_model <- lm(polymr_model$exposure ~ polymr_model$genotypes)
  polymr_model$beta_exposure <- coef(full_exposure_model)[-1]

  if (!is.null(polymr_model$reverse_t_thr))
    polymr_model$se_exposure <- summary(full_exposure_model)$coefficients[-1, 2]

  polymr_model
}
