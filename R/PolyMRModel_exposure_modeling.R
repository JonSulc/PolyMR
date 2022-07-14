model_exposure <- function(polymr_model) {
  polymr_model <- create_exposure_model(polymr_model)

  if (!is.null(polymr_model$reverse_t_thr))
    polymr_model <- filter_reverse_snps(polymr_model)

  polymr_model
}


create_exposure_model <- function(polymr_model) {
  polymr_model$exposure_model <-
    lm(polymr_model$exposure ~ polymr_model$genotypes)
  polymr_model
}


filter_reverse_snps <- function(polymr_model) {
  genotype_outcome_stats <-
    summary(lm(polymr_model$outcome ~ polymr_model$genotypes))$coefficients[-1, 1:2]
  genotype_exposure_stats <-
    summary(polymr_model$exposure_model)$coefficients[-1, 1:2]

  to_keep <-
    (abs(genotype_exposure_stats[, 1]) - abs(genotype_outcome_stats[, 1])) /
    sqrt(genotype_exposure_stats[, 2]^2 + genotype_outcome_stats[, 2]^2) >
    polymr_model$reverse_t_thr
  
  if (all(to_keep))
    return(polymr_model)

  polymr_model$genotypes <- polymr_model$genotypes[, to_keep]
  create_exposure_model(polymr_model)
}
