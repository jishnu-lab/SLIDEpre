#' Calculate the element-wise mean squared error of a matrix.
#'
#' @param A a matrix
#' @param B a matrix of same dimensions as \code{A}
#' @return MSE

Mse <- function(A, B) {
  return(mean((A - B) ** 2))
}
