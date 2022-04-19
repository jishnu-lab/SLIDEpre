#' Adjust Matrix
#'
#' Find posterior mean of \eqn{\Sigma} given data matrix, \eqn{x}, prior correlation matrix,
#' \eqn{\Delta}, and weight, \eqn{\alpha}.
#'
#' @param x standardized data matrix of dimensions \eqn{n \times p}
#' @param alpha weight to use in average between sample correlation matrix and prior correlation matrix
#' @param Delta prior correlation matrix of dimensions \eqn{p \times p} (positive semi-definite)
#' @return An object of class \sQuote{data.frame}
#' @export

adjMat <- function(x, Delta, alpha) {
  n <- nrow(x)
  adj_mat <- (alpha * Delta) + ((1 - alpha) * (stats::cor(x)))
  return (adj_mat)
}
