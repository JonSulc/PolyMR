new_PolyMRModel <- function(
    exposure_outcome_model_data = NULL,
    exposure = exposure_outcome_model_data$exposure,
    outcome = exposure_outcome_model_data$outcome,
    genotypes = exposure_outcome_model_data$genotypes,
    beta_exposure = exposure_outcome_model_data$beta_exposure,
    se_exposure = exposure_outcome_model_data$se_exposure,
    control_function = residuals(exposure_outcome_model_data$exposure_model),
    exposure_powers = 1,
    control_function_powers = 1:min(max(exposure_powers),
                                    max_control_function_power),
    reverse_t_thr = NULL,
    max_exposure_power = 10,
    max_control_function_power = NULL
  ){
  if ((is.null(exposure_outcome_model_data$exposure)  && is.null(exposure)) ||
      (is.null(exposure_outcome_model_data$outcome)   && is.null(outcome))  ||
      (is.null(exposure_outcome_model_data$genotypes) && is.null(genotypes)))
    stop("Exposure, outcome and genotypes must all be provided.")

  if (length(exposure) != length(outcome) || length(exposure) != nrow(genotypes))
    stop("Exposure, outcome and genotype sample sizes must match.")

  structure(
    list(
      exposure = exposure,
      outcome = outcome,
      genotypes = genotypes,
      beta_exposure = beta_exposure,
      se_exposure = se_exposure,
      control_function = control_function,
      exposure_powers = exposure_powers,
      control_function_powers = control_function_powers,
      reverse_t_thr = reverse_t_thr,
      max_exposure_power = max_exposure_power,
      max_control_function_power = max_control_function_power
    ),
    class = c("PolyMRModel", "EOModel")
  ) |>
    preprocess_data() |>
    model()
}


model.PolyMRModel <- function(polymr_model){
  model_exposure(polymr_model) |>
    model_outcome()
}













