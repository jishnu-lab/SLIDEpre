#' Calculate indices of each row such that the absolute values of these indices
#' are within \eqn{2\delta} of the maximal absolute value \eqn{M} of this row.
#' This is an implementation of step 4 of Algorithm 1 in Bing et al. (2020).
#'
#' @param i the row index
#' @param M the maximal absolute value of each row of the covariance/correlation matrix
#' @param arg_M the first index at which each entry in \eqn{M} is achieved
#' @param vector a row of \code{off_Sigma}
#' @param delta \eqn{\delta}, a numeric constant
#' @param se_est vector of standard deviations of features (columns of \code{X})
#' @return a vector of indices

FindRowMaxInd <- function(i, M, arg_M, vector, delta, se_est) {
  ## lbd <- delta * sd(feat i) * sd(feat max M) + delta * sd(feat i) * sd(each feat)
  ## lbd is a vector of values for each column off_Sigma (each feature)
  lbd <- delta * se_est[i] * se_est[arg_M] + delta * se_est[i] * se_est

  ## which indices of M (the max abs values for each row in off_Sigma) are ≤ lbd + value in vector
  ## ALG 1.4 - find set of indices within 2delta of max
  indices <- which(M <= lbd + vector)
  return(indices)
}
