#' For given row, check if it is a pure node by iteratively checking the nodes
#' in \code{Si}. Return TRUE if the given row corresponds to a pure variable.
#'
#' @param Sigma a matrix of dimensions \eqn{p \times p}
#' @param rowInd the index of the row
#' @param Si a vector of indices
#' @param Ms a vector of the largest absolute values of each of the rows in \code{Sigma}
#' @param delta \eqn{\delta}, a numerical constant
#' @return TRUE or FALSE

TestPure <- function(Sigma_row, rowInd, Si, Ms, arg_Ms, delta, se_est) {
  for (i in 1:length(Si)) {
    #### go row by row through list of indices where abs val is within 2*delta of max
    j <- Si[i]

    delta_j <- (se_est[rowInd] + se_est[arg_Ms[j]]) * se_est[j] * delta
    if (abs(Sigma_row[j] - Ms[j]) > delta_j) {
      return(FALSE)
    }
  }
  return(TRUE)
}
