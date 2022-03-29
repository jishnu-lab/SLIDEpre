#' Calculate the sum of squares of the upper off-diagonal elements of two matrices
#'
#' @param M a matrix
#' @param N a matrix of same dimensions as \code{M}
#' @param weights the weights to use in calculating the sum
#' @return the sum of squares

offSum <- function(M, N, weights) {
  tmp <- (M-N) / weights
  tmp <- t(t(tmp) / weights)
  return(sum((tmp[row(tmp) <= (col(tmp) - 1)])^2))
}
