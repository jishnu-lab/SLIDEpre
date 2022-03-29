#' Transform \eqn{\Phi} to \eqn{\alpha} for log-likelihood calculation
#'
#' @param phi numeric value
#' @param n number of rows/samples in data matrix (integer)
#' @param p number of columns/features in data matrix (integer)
#' @return numeric

toAlpha <- function(phi, n, p) {
  alpha <- (phi - p - 1) / (n + phi - p - 1)
  return (alpha)
}
