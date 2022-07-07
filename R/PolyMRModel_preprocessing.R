preprocess_data.PolyMRModel <- function(polymr_model){
  NextMethod() |>
    remove_duplicate_snps() |>
    filter_reverse_snps()
}

remove_duplicate_snps <- function(polymr_model){
  polymr_model <- calculate_beta_exposure(polymr_model)
  if (any(is.na(polymr_model$beta_exposure))) {
    warning("Some IVs are not independent of one another, removing duplicates.")
    polymr_model <-
      filter_snps(polymr_model, !is.na(polymr_model$beta_exposure))
  }
  polymr_model
}

filter_reverse_snps <- function(polymr_model){
  if (is.null(polymr_model$reverse_t_thr))
    return(polymr_model)

  polymr_model <- calculate_beta_exposure(polymr_model)
  genotype_outcome_stats <-
    summary(lm(
        polymr_model$outcome ~ polymr_model$genotypes
      ))$coefficients[-1, 1:2]

  to_keep <-
    (abs(polymr_model$beta_exposure) - abs(genotype_outcome_stats[, 1])) /
    sqrt(polymr_model$se_exposure^2 + genotype_outcome_stats[, 2]^2) > polymr_model$reverse_t_thr

  filter_snps(polymr_model, to_keep)
}

filter_snps <- function(polymr_model,
                        to_keep){
  polymr_model$genotypes <-
    polymr_model$genotypes[, to_keep]
  polymr_model$beta_exposure <-
    polymr_model$beta_exposure[to_keep]
  polymr_model$se_exposure[to_keep]

  polymr_model
}
