#' For given row, check if it is a pure node by iteratively checking the nodes
#' in \code{Si}. Return TRUE if the given row corresponds to a pure variable.
#' This is an implementation of steps 6-9 of Algorithm 1 in Bing et al. (2020).
#'
#' @param Sigma_row a vector of dimension \eqn{p}
#' @param rowInd the index of the row
#' @param Si a vector of indices that are within \eqn{2\delta} of the maximum value of \code{Sigma_row}
#' @param Ms a vector of the largest absolute values of each of the rows in \code{Sigma}
#' @param arg_Ms a vector of the first index at which each value in Ms is achieved
#' @param delta \eqn{\delta}, a numerical constant
#' @return TRUE or FALSE

TestPure <- function(Sigma_row, rowInd, Si, Ms, arg_Ms, delta, se_est) {
  for (i in 1:length(Si)) {
    #### go row by row through list of indices where abs val is within 2*delta of max
    j <- Si[i] #### j is some index of a row in Sigma that has large enough abs cov/corr
    #### some sort of cutoff value
    delta_j <- (se_est[rowInd] + se_est[arg_Ms[j]]) * se_est[j] * delta
    #### check if abs diff between jth entry in row i and max value in row j is greater than cutoff
    #### Sigma_row[j] = Sigma_{ij} (ith by jth entry)
    #### Ms[j] = max(Sigma_{j.}) (max of jth row)
    if (abs(Sigma_row[j] - Ms[j]) > delta_j) {
      return(FALSE)
    }
  }
  return(TRUE)
}
