#' Calculate the mean element-wise absolute deviation between two matrices.
#'
#' @param mat_A a matrix
#' @param mat_B a matrix of same dimensions as \code{mat_A}
#' @return mean absolute deviation

devA <- function(mat_A, mat_B) {
  return (mean(abs(mat_A - mat_B)))
}
