#' Transform \eqn{\alpha} to \eqn{\Phi} for log-likelihood calculation
#'
#' @param alpha numeric value
#' @param n number of rows/samples in data matrix (integer)
#' @param p number of columns/features in data matrix (integer)
#' @return numeric

toPhi <- function(alpha, n, p) {
  phi <- ((alpha * n) + ((1 - alpha) * p) + (1 - alpha)) / (1 - alpha)
  phi2 <- (alpha * n) / (1 - alpha) + p + 1
  return (phi)
}
