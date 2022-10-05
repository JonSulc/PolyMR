model_exposure <- function(polymr_model) {
  polymr_model <- create_exposure_model(polymr_model)

  if (!is.null(polymr_model$reverse_t_thr))
    polymr_model <- filter_reverse_snps(polymr_model)

  polymr_model
}


create_exposure_model <- function(polymr_model) {
  if (is.null(polymr_model$beta_exposure))
    polymr_model$exposure_model <-
      lm(polymr_model$exposure ~ polymr_model$genotypes)

  polymr_model |>
    calculate_control_function()
}


calculate_control_function <- function(polymr_model) {
  if (!is.null(polymr_model$beta_exposure)){
    polymr_model$cf <- with(polymr_model,
                            {
                              exposure - genotypes %*% beta_exposure
                            }) |>
      c()
  } else
    polymr_model$cf <- residuals(polymr_model$exposure_model)

  polymr_model
}


filter_reverse_snps <- function(polymr_model) {
  if (!is.null(polymr_model$beta_exposure) &&
      !is.null(polymr_model$se_exposure)) {
    beta_exposure <- polymr_model$beta_exposure
    se_exposure <- polymr_model$se_exposure
  } else {
    beta_exposure <- coef(polymr_model$exposure_model)[-1]
    se_exposure <- summary(polymr_model$exposure_model)$coef[-1, 2]
  }
  genotype_outcome_stats <-
    summary(lm(polymr_model$outcome ~ polymr_model$genotypes))$coefficients[-1, 1:2]

  to_keep <-
    (abs(beta_exposure) - abs(genotype_outcome_stats[, 1])) /
    sqrt(se_exposure^2 + genotype_outcome_stats[, 2]^2) >
    polymr_model$reverse_t_thr

  if (all(to_keep))
    return(polymr_model)

  polymr_model$genotypes <- polymr_model$genotypes[, to_keep]
  create_exposure_model(polymr_model)
}
