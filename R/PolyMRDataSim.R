#' Simulate exposure, outcome, and genotype data
#'
#' Simulates exposure, outcome, and genotype data corresponding to the provided
#' causal function.
#'
#' @param sample_size Sample size (default is 10^5).
#' @param n_exposure_snps Number of SNPs explaining the
#'   \code{exposure_heritability} (default is 100). Note that this is not equal
#'   to the number of instruments as it includes SNPs which will be filtered out
#'   upon finalizing the data (see \code{gws_thr}).
#' @param exposure_heritability Heritability of the exposure explained by the
#'   \code{n_exposure_snps} (default is 0.3).
#' @param causal_function Function defining the true relationship between the
#'   exposure and the outcome. It should accept a vector of exposure values and
#'   return a vector of outcome values of the same length. This represents the
#'   pure contribution of the exposure to the outcome and should not include
#'   confounding or noise. Default is \code{function(x) 0.1*x + 0.05*x^2}. See
#'   also [get_polynomial_function()]
#' @param confounders_list A list of objects of class \code{Confounder} (see
#'   [new_Confounder()]) to be added to the data. Default is a single confounder
#'   with linear contributions to both exposure and outcome with coefficients
#'   0.2 and 0.5, respectively.
#' @param finalize Logical indicating whether the data set should be finalized,
#'   i.e. errors added to contribute remaining variance and SNPs filtered based
#'   on genome-wide significance (default is TRUE).
#' @param gws_thr P-value genome-wide significance threshold to filter SNPs for
#'   instrumental variable selection (default is 5e-8).
#'
#' @return A list-like object of class \code{PolyMRDataSim}. It's main
#'   constituents are the \code{exposure} and \code{outcome} vectors, and the
#'   \code{genotypes} matrix. In addition, a number of parameters used in data
#'   generation are kept as named elements, including \code{n_exposure_snps},
#'   \code{exposure_heritability}, and \code{causal_function}. Some intermediate
#'   values are also included, namely the minor allele frequencies (\code{mafs})
#'   and the effects on the exposure (\code{exposure_coefficients}) of the SNPs
#'   remaining after filtering. These represent the ground-truth values used in
#'   the data generation.
#'
#' @examples
#'   simulated_data <- new_PolyMRDataSim(
#'     sample_size = 50000,
#'     n_exposure_snps = 200,
#'     causal_function = function(x) 0.05*exp(x))
#'
new_PolyMRDataSim <- function(
    sample_size = 1e5,
    n_exposure_snps = 100,
    exposure_heritability = 0.3,
    causal_function = get_polynomial_function(c(.1, .05)),
    confounders_list = list(new_Confounder(sample_size)),
    finalize = TRUE,
    gws_thr = 5e-8) {
  stopifnot(is.numeric(sample_size))
  stopifnot(is.numeric(n_exposure_snps))
  stopifnot(is.double(exposure_heritability))
  stopifnot(is.function(causal_function))
  stopifnot(is.list(confounders_list) | "Confounder" %in% class(confounders_list))
  stopifnot(is.logical(finalize))

  polymr_data <- structure(
    list(
      n_exposure_snps = n_exposure_snps,
      exposure_heritability = exposure_heritability,
      causal_function = causal_function
    ),
    class = "PolyMRDataSim"
  )
  polymr_data <- generate_exposure_snp_mafs(polymr_data)
  polymr_data <- generate_exposure_snp_genotypes(polymr_data, sample_size)
  polymr_data <- remove_invariant_snps(polymr_data)
  polymr_data <- generate_exposure_coefficients(polymr_data)
  polymr_data <- generate_genetic_exposure(polymr_data)

  for (confounder in confounders_list)
    polymr_data <- polymr_data + confounder

  if (finalize)
    polymr_data <- finalize_data(polymr_data, gws_thr = gws_thr)
  else
    polymr_data <- recalculate_outcome(polymr_data)

  polymr_data
}

validate_PolyMRDataSim <- function(polymr_data){
  if (!(length(polymr_data$exposure) == length(polymr_data$outcome) &&
        length(polymr_data$exposure) == nrow(polymr_data$genotypes)))
    stop(
      "Sample sizes for 'exposure', 'outcome', and 'genotypes' must match.",
      call = FALSE
    )
  if (length(polymr_data$mafs) != length(polymr_data$exposure_coefficients) ||
      length(polymr_data$mafs) != ncol(polymr_data$genotypes))
    stop(
      "The size of 'mafs' and 'exposure_coefficients' must match the dimensions of the 'genotypes' matrix.",
      call. = FALSE
    )
  if (var(polymr_data$exposure) > 1.01)
    stop(
      "The variance of the 'exposure' (standardized) cannot exceed 1.",
      call. = FALSE
    )
  if (var(polymr_data$outcome) > 1.01)
    stop(
      "The variance of the 'outcome' (standardized) cannot exceed 1.",
      call. = FALSE
    )
}

`[.PolyMRDataSim` <- function(polymr_data, sample_index, snp_index){
  polymr_data$exposure  <- polymr_data$exposure[sample_index]
  polymr_data$outcome   <- polymr_data$outcome[sample_index]
  polymr_data$genotypes <- polymr_data$genotypes[sample_index, snp_index, drop = FALSE]
  polymr_data$mafs <- polymr_data$mafs[snp_index]
  polymr_data$exposure_coefficients <- polymr_data$exposure_coefficients[snp_index]

  polymr_data
}

generate_exposure_snp_mafs <- function(polymr_data){
  polymr_data$mafs <- rbeta(polymr_data$n_exposure_snps, 1, 3)

  while (any(invalid_mafs <- polymr_data$mafs > .5))
    polymr_data$mafs[invalid_mafs] <- rbeta(sum(invalid_mafs), 1, 3)

  polymr_data
}

generate_exposure_snp_genotypes <- function(polymr_data,
                                            sample_size){
  polymr_data$genotypes <-
    apply(matrix(polymr_data$mafs, nrow = 1),
          2,
          function(maf) rbinom(sample_size, 2, maf))

  remove_invariant_snps(polymr_data)
}

remove_invariant_snps <- function(polymr_data){
  genotype_variance <- apply(polymr_data$genotypes, 2, var, na.rm = TRUE)
  invariant_positions <- is.na(genotype_variance) | genotype_variance == 0

  polymr_data[, !invariant_positions]
}

generate_exposure_coefficients <- function(polymr_data){
  polymr_data$exposure_coefficients <- rnorm(
    length(polymr_data$mafs),
    mean = 0,
    sd   = 2 * (polymr_data$mafs * (1-polymr_data$mafs))^(3 / 8)
  )

  polymr_data$exposure_coefficients <-
    polymr_data$exposure_coefficients *
    sqrt(polymr_data$exposure_heritability /
         sum(polymr_data$exposure_coefficients^2))

  polymr_data
}

generate_genetic_exposure <- function(polymr_data){
  polymr_data$exposure <-
    c(scale(polymr_data$genotypes)[,] %*% polymr_data$exposure_coefficients)

  polymr_data
}


`+.PolyMRDataSim` <- function(polymr_data, confounder_obj){
  polymr_data$exposure <- polymr_data$exposure +
    confounder_obj$exposure_confounding_function(
      confounder_obj$confounder_values)
  polymr_data$confounders_list <- c(polymr_data$confounders_list,
                                    list(confounder_obj))

  polymr_data <- recalculate_outcome(polymr_data)
  validate_PolyMRDataSim(polymr_data)

  polymr_data
}


recalculate_outcome <- function(polymr_data){
  polymr_data$outcome <- polymr_data$causal_function(polymr_data$exposure)

  apply_confounders_to_outcome(polymr_data)
}

finalize_data <- function(polymr_data,
                          filter_snps = TRUE,
                          gws_thr = 5e-8){
  polymr_data <- add_error_to(polymr_data, "exposure")
  polymr_data <- recalculate_outcome(polymr_data)
  polymr_data <- add_error_to(polymr_data, "outcome")

  if (filter_snps)
    return(filter_gws_snps(polymr_data, gws_thr))

  polymr_data
}

add_error_to <- function(polymr_data,
                         phenotype = "exposure"){
  error_values <- scale(rnorm(length(polymr_data[[phenotype]])))[,]
  remaining_variance <- 1 - var(polymr_data[[phenotype]])
  if (remaining_variance < 0)
    stop(sprintf("The variance of %s exceeds 1 before adding noise.",
                 phenotype))
  polymr_data[[phenotype]] <- polymr_data[[phenotype]] +
    error_values * sqrt(remaining_variance)

  polymr_data
}

apply_confounders_to_outcome <- function(polymr_data){
  for (confounder in polymr_data$confounders_list)
    polymr_data$outcome <- polymr_data$outcome +
      confounder$outcome_confounding_function(confounder$confounder_values)

  polymr_data
}

filter_gws_snps <- function(polymr_data,
                            gws_thr = 5e-8){
  summary_stats <- summary(lm(
      polymr_data$exposure~polymr_data$genotypes
    ))$coefficients
  polymr_data[ , summary_stats[ -1, "Pr(>|t|)" ] < gws_thr]
}











