#' Function factory to create polynomial functions
#'
#' This function factory is intended for use with the PolyMRDataSim class. It
#' simplifies the creation of polynomial functions based on its coefficients.
#'
#' @param polynomial_coefficients A vector of coefficients for the polynomial
#'   terms, in ascending order of degree (default is c(0.1)). The length of the
#'   vector determines the degree of the resulting polynomial.
#'
#' @return A polynomial function with the corresponding coefficients.
#'
#' @examples
#'   cubic_polynomial <- get_polynomial_function(c(0.5, 0, 2))
#'
#'   some_values <- rnorm(100)
#'   all(cubic_polynomial(some_values) == 0.5*some_values + 2*some_values^3)
#'
#' @export
get_polynomial_function <- function(polynomial_coefficients = c(.1)){
  # Argument evaluation must be forced to avoid lazy evaluation potentially
  # occurring after the value has changed.
  force(polynomial_coefficients)
  function(values){
    powers <- outer(values,
                    1:length(polynomial_coefficients),
                    "^")
    c(powers %*% polynomial_coefficients)
  }
}
