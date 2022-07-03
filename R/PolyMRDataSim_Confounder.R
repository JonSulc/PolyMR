#' Create a confounder for simulated data
#'
#' This function defines the \code{Confounder} class for use with
#' \code{PolyMRDataSim}. It is a simplified way to add a confounder to an
#' exposure-outcome pair, i.e. correlated error between two simulated
#' phenotypes.
#'
#' @param sample_size Sample size of the simulated exposure-outcome data set.
#' @param exposure_confounding_function A function that defines the confounder
#'   contribution to the exposure (default is \code{function(x) 0.2*x}). The
#'   function should accept a vector of confounder values (standardized) and
#'   return a vector of values to add to the un-confounded exposure.
#' @param outcome_confounding_function A function that defines the confounder
#'   contribution to the outcome (default is \code{function(x) 0.5*x}). The
#'   function should accept a vector of confounder values (standardized) and
#'   return a vector of values to add to the un-confounded outcome.
#'
#' @return A list-like object of class \code{Confounder}, containing a vector of
#'   confounder values (\code{confounder_values}) and functions to calculate the
#'   confounder contribution to exposure and outcome (respectively
#'   \code{exposure_confounding_function} and
#'   \code{outcome_confounding_function}). Given an object \code{polymr_data} of
#'   class \code{PolyMRDataSim}, the confounder can be applied through simple
#'   addition, i.e. \code{polymr_data + my_confounder}.
#'
#' @examples
#'   my_confounder <- new_Confounder(1000,
#'     exposure_confounding_function = function(x) 0.3*x,
#'     outcome_confounding_function  = function(x) 0.1*x + 0.05*x^2
#'   )
#' @export
new_Confounder <- function(
    sample_size,
    exposure_confounding_function = get_polynomial_function(0.2),
    outcome_confounding_function  = get_polynomial_function(0.5)
  ){

  structure(
    list(confounder_values = rnorm( sample_size ) |>
           scale() |>
           c(),
         exposure_confounding_function = exposure_confounding_function,
         outcome_confounding_function  = outcome_confounding_function),
    class = "Confounder"
  )
}

