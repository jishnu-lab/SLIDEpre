#' Estimate list of pure node indices for given \eqn{\Sigma} and \eqn{\delta}.
#' This code is an implementation of Algorithm 1 from Bing et al. (2020).
#'
#' @param off_Sigma a sample correlation matrix of dimensions \eqn{p \times p}
#' @param delta \eqn{\delta}, a numerical constant
#' @param Ms the largest absolute values of each row of \code{off_Sigma}
#' @param arg_Ms vector of column indices at which the values in \code{Ms} are achieved in \code{off_Sigma}
#' @param se_est standard deviations of features (columns of \code{X})
#' @param merge boolean indicating merge style
#' @return a list including the list of estimated pure node indices and a vector
#' of the estimated pure node indices

FindPureNode <- function(off_Sigma, delta, Ms, arg_Ms, se_est, merge) {
  G <- list() #### groups

  for (i in 1:nrow(off_Sigma)) { #### loop through rows
    row_i <- off_Sigma[i,]
    #### Calculate indices of ith row such that the absolute values of these indices
    #### are within 2 * delta from the maximal absolute value (M) of this row.
    #### these values correspond to columns of off_Sigma
    Si <- FindRowMaxInd(i, Ms[i], arg_Ms[i], row_i, delta, se_est) ## ALG 1.4
    if (length(Si) != 0) {
      #### For given row, check if it is a pure node by iteratively checking the other
      #### indices in Si. Return TRUE if the given row corresponds to a pure variable.
      pureFlag <- TestPure(row_i, i, Si, Ms, arg_Ms, delta, se_est) ## ALG 1.7
      if (pureFlag) { #### if pure
        if (merge) {
          G <- Merge(G, c(Si, i))
        } else {
          G <- Merge_union(G, c(Si, i))
        }
      }
    }
  }
  return(list(pureInd = G, pureVec = unlist(G)))
}
