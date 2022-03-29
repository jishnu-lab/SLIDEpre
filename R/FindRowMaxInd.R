#' Calculate indices of each row such that the absolute values of these indices
#' are within \eqn{2\delta} of the maximal absolute value \eqn{M} of this row.
#'
#' @param i the row index
#' @param M the maximal absolute value of row \eqn{i}
#' @param arg_M the indices at which \eqn{M} is achieved in the row
#' @param vector the row
#' @param delta \eqn{\delta}, a numeric constant
#' @param se_est ???
#' @return a vector of indices

FindRowMaxInd <- function(i, M, arg_M, vector, delta, se_est) {
  lbd <- delta * se_est[i] * se_est[arg_M] + delta * se_est[i] * se_est
  indices <- which(M <= lbd + vector)
  return(indices)
}
