#' Mendelian randomization-based approximation of non-linear causal effects
#'
#' This function approximates a non-linear causal effect through a polynomial
#' regression of observational data, correcting for confounding using an
#' instrumental variable-based approach.
#'
#' @param exposure A vector containing the exposure values for each individual.
#' @param outcome A vector containing the outcome values for each individual.
#' @param genotypes The NxM genetic matrix, with a column for each variant and a
#'   row for each individual.
#' @param return_phenotypes_summary Whether to return a data.table containing
#'   the median, mean, and standard deviation of both exposure and outcome
#'   (default is TRUE).
#' @param return_observational_function Whether to return a polynomial
#'   approximation of the observed association between exposure and outcome
#'   (default is TRUE).
#' @param return_binned_observations Whether to return a data.table containing
#'   per-bin summary information, including the median exposure and the median,
#'   mean, and standard deviation of the outcome, binned on exposure (default is
#'   TRUE).
#' @param bins Number of bins for which to return mean and median values
#'   (default is 100).
#' @param starting_exposure_powers A vector containing the exponents for the
#'   exposure terms in the initial model. Default is c(1:10), corresponding to a
#'   10th degree polynomial with all lower terms present.
#' @param max_exposure_power The maximum exponent to use in modeling the
#'   exposure (default is \code{max(starting_exposure_powers)}). If this is
#'   greater greater than \code{max(starting_exposure_powers)}, \code{polymr}
#'   will iteratively increase (by \code{power_step}) the degree of the causal
#'   polynomial function as long as new terms are significant (p <
#'   \code{p_thr_add}).
#' @param max_control_function_power The maximum exponent to use in modeling the
#'   control function component. Default is NULL, in which case the control
#'   function polynomial will include all terms from 1 to the highest degree of
#'   the exposure component.
#' @param power_step The number by which to increment the degree of the exposure
#'   polynomial each iteration until \code{max_exposure_power} is reached or the
#'   new terms are no longer significant (p > \code{p_thr_add}). Default is 2,
#'   as even and odd degree terms have different properties.
#' @param reverse_t Threshold to use for reverse causality filtering (T
#'   statistic), NULL for no filtering (default). A value of 0 represents a
#'   simple filtering out of IVs explaining more variance in the outcome than
#'   the exposure, whereas a value of 1.645 (\code{qnorm(.95)}) would remove
#'   only those where that difference is significant (p < 0.05).
#' @param p_thr_add The p-value threshold determining if newly added exposure
#'   terms should be considered significant enough to further increase the
#'   degree of the polynomial (by \code{power_step}, up to
#'   \code{max_exposure_power}). Default is 0, which will prevent new terms from
#'   being added. A value of NULL will be equivalent to a per-step Bonferroni
#'   correction, i.e. 0.05 / `power_step`.
#' @param p_thr_drop The p-value threshold determining which, if any, exposure
#'   terms should be dropped from the final function. This is done iteratively
#'   and the significance of each term is assessed in the new context before
#'   proceeding again, if necessary, until all remaining terms reach the defined
#'   significance threshold. Default is 1, which will retain all terms. A value
#'   of NULL will use a Bonferroni-corrected threshold at each step.
#' @param drop_higher_control_function_powers Logical indicating whether control
#'   function terms with a higher degree than the highest exposure term should
#'   be dropped. Default is TRUE. Only relevant if \code{p_thr_drop} < 1 or is
#'   NULL.
#'
#'
#' @return Returns a named list of results for PolyMR itself and the other
#'   selected values:
#'   \itemize{
#'   \item \code{phenotypes_summary} is a data.table with the median, mean, and
#'     standard deviation of both exposure and outcome
#'   \item \code{binned_observations} is a data.table with per-bin summary
#'     information, including the median exposure and the median, mean, and
#'     standard deviation of the outcome (binned on the exposure).
#'   \item \code{observational} is a list-like object of class \code{EOModel}
#'     containing:
#'     \itemize{
#'       \item \code{outcome_model}, an object of class \code{lm} containing the
#'         full model. Use \code{summary()} for more details.
#'       \item \code{vcov}, the variance-covariance matrix which can be used to
#'         create the 95% confidence hull for plotting.
#'       \item \code{pval_null_model}, the p-value for the full model
#'         (F-statistic-based).
#'       \item \code{pval_linear_model}, the LRT p-value comparing the full
#'         model to the linear model.
#'       \item \code{r_squared}, the variance explained by the model.
#'     }
#'   \item \code{polymr} is a list-like object of class \code{PolyMRModel}, the
#'     contents of which are similar to those of \code{observational}:
#'     \itemize{
#'       \item \code{outcome_model}, an object of class \code{lm} containing the
#'         full model. Use \code{summary()} for more details.
#'       \item \code{vcov}, the variance-covariance matrix of the coefficient
#'         estimates, which can be used to create the 95% confidence interval
#'         for plotting.
#'       \item \code{pval_null_model}, the LRT p-value comparing the full model
#'         to the model with (all) the control function terms but no exposure
#'         terms.
#'       \item \code{pval_linear_model}, the LRT p-value comparing the full
#'         model to the linear model containing all control function terms but
#'         only the degree 1 (linear) exposure term.
#'       \item \code{r_squared}, the variance of the outcome attributable to the
#'         causal effect of the exposure. This is obtained by comparing the
#'         R-squared of the full model to that of the null model (containing
#'         only the control function terms).
#'     }
#'   }
#' @examples
#' simulated_data <- PolyMR:::new_PolyMRDataSim()
#' polymr(exposure  = simulated_data$exposure,
#'        outcome   = simulated_data$outcome,
#'        genotypes = simulated_data$genotypes,
#'        reverse_t_thr = 0,
#'        p_thr_drop = NULL)
#'
#'
#' @import data.table
#'
#' @export
polymr <- function(exposure,
                   outcome,
                   genotypes,
                   return_phenotypes_summary = TRUE,
                   return_observational_function = TRUE,
                   return_binned_observations = TRUE,
                   bins = 100,
                   starting_exposure_powers = 1:10,
                   max_exposure_power = max(starting_exposure_powers),
                   max_control_function_power = NULL,
                   power_step = 2,
                   reverse_t_thr = NULL,
                   p_thr_add = 0,
                   p_thr_drop = 1,
                   drop_higher_control_function_powers = TRUE){
  if (!is.null(max_control_function_power) &&
      max_control_function_power < max_exposure_power)
    warning("max_control_function_power < max_exposure_power: This is not recommended.")

  results <- list()

  if (return_phenotypes_summary)
    results$phenotypes_summary <- get_phenotypes_summary(exposure, outcome)

  if (return_binned_observations) {
    results$binned_observations <- get_bins(exposure = exposure,
                                            outcome = outcome,
                                            bins = bins)
    results$binned_observations_scaled <- get_bins(exposure = scale(exposure)[,],
                                                   outcome = scale(outcome)[,],
                                                   bins = bins)
  }

  if (return_observational_function)
    results$observational <- new_EOModel(
        exposure = exposure,
        outcome = outcome,
        exposure_powers = starting_exposure_powers,
        max_exposure_power = max_exposure_power
      ) |>
      select_best_outcome_model(
        power_step = power_step,
        p_thr_add = p_thr_add,
        p_thr_drop = p_thr_drop
      ) |>
      return_model()

  results$polymr <-
    new_PolyMRModel(exposure = exposure,
                    outcome = outcome,
                    genotypes = genotypes,
                    exposure_powers = starting_exposure_powers,
                    reverse_t_thr = reverse_t_thr,
                    max_exposure_power = max_exposure_power,
                    max_control_function_power = max_control_function_power) |>
    select_best_outcome_model(
        power_step = power_step,
        p_thr_add = p_thr_add,
        p_thr_drop = p_thr_drop,
        drop_higher_control_function_powers = drop_higher_control_function_powers
      ) |>
      return_model()

  results
}


get_phenotypes_summary <- function(exposure, outcome){
  data.table::data.table(
    phenotype = c("Exposure", "Outcome"),
    median = c(median(exposure), median(outcome)),
    mean   = c(mean(exposure),   mean(outcome)),
    sd     = c(sd(exposure),     sd(outcome))
  )
}


select_best_outcome_model <- function(eo_model,
                                      power_step = 2,
                                      p_thr_add = 0,
                                      p_thr_drop = 1,
                                      drop_higher_control_function_powers = TRUE){
  new_eo_model <- eo_model
  new_terms_are_significant <- any(are_terms_significant(new_eo_model,
                                                         p_thr_add,
                                                         power_step))
  while (new_terms_are_significant){
    eo_model <- new_eo_model
    if (max(eo_model$exposure_powers) < eo_model$max_exposure_power){
      new_eo_model <- increase_order_of_outcome_model(eo_model,
                                                      power_step = power_step)
      new_terms_are_significant <- any(are_terms_significant(new_eo_model,
                                                             p_thr_add,
                                                             power_step))
    } else{
      message("Max power reached")
      break
    }
  }

  while (any(!are_terms_significant(eo_model, p_thr_drop)) && eo_model$exposure_powers > 0)
    eo_model <-
      drop_least_significant_term(
        eo_model,
        drop_higher_control_function_powers = drop_higher_control_function_powers
      )
  eo_model
}


increase_order_of_outcome_model <- function(x, ...){
  UseMethod("increase_order_of_outcome_model")
}

drop_least_significant_term <- function(x, ...){
  UseMethod("drop_least_significant_term")
}

are_terms_significant <- function(x, ...){
  UseMethod("are_terms_significant")
}














