#' Natural Logarithm Of Multivariate Gamma Function
#'
#' Calculate natural logarithm of multivariate Gamma function for evaluating the
#' marginal log-likelihood of \eqn{x}. This function both adds and subtracts values
#' as in the log-likelihood of \eqn{x}. See \code{\link{logLik}} for more details.
#'
#' @param p for p-adic Gamma function
#' @param val1 first value
#' @param val2 second value
#' @return the natural logarithm of multivariate Gamma function

logGamma <- function(p, val1, val2) {
  base <- 0.0
  for (i in 1:p) {
    gamma_diff <- lgamma(val1 + (1 - i) * 0.5) - lgamma(val2 + (1 - i) * 0.5)
    base <- base + gamma_diff
  }
  return (base)
}
