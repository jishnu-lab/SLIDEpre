#' Calculate the fitted value of \eqn{A_I \cdot C \cdot A_I^\top} for given
#' \eqn{\Sigma} and \eqn{\delta}.
#'
#' @param Sigma a correlation matrix of dimensions \eqn{p \times p}
#' @param delta a threshold parameter
#' @param Ms the calculated maximal values of \eqn{\Sigma} by row
#' @param se_est the estimated standard errors
#' @param diagonal a boolean indicating the diagonal structure
#' @param merge a boolean indicating the merge style
#' @return a list containing a vector of the indices of the estimated pure variables and
#' the fitted value of \eqn{A_I \cdot C \cdot A_I^\top}

CalFittedSigma <- function(Sigma, delta, Ms, arg_Ms, se_est, diagonal, merge) {
  resultPureNode <- FindPureNode(abs(Sigma), delta, Ms, arg_Ms, se_est, merge)

  estPureIndices <- resultPureNode$pureInd
  # lapply(estPureIndices, function(x) cat(x, "\n"))

  if (singleton(estPureIndices))
    return(list(pureVec = NULL, fitted = -1))

  estSignPureIndices <- FindSignPureNode(estPureIndices, Sigma)
  AI <- RecoverAI(estSignPureIndices, length(se_est))
  C <- EstC(Sigma, AI, diagonal)

  if (length(estPureIndices) == 1)
    fitted <- -1
  else {
    subAI <- AI[resultPureNode$pureVec, ]
    fitted <- subAI %*% C %*% t(subAI)
  }
  return(list(pureVec = resultPureNode$pureVec, fitted = fitted))
}
