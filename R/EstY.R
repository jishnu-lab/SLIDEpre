#' Estimate the matrix, \eqn{\Sigma_{\hat{TJ}}}.
#'
#' @param Sigma a sample correlation matrix of dimensions \eqn{p \times p}
#' @param AI a matrix of dimensions \eqn{p \times K}
#' @param pureVec a vector of indices of the pure nodes
#' @return an estimate for \eqn{\Sigma_{\hat{TJ}}} of dimensions \eqn{K \times |J|}

EstY <- function(Sigma, AI, pureVec) {
  SigmaS <- AdjustSign(Sigma, AI)
  SigmaJ <- matrix(SigmaS[ ,-pureVec], nrow = nrow(Sigma))
  SigmaTJ <- matrix(0, ncol(AI), nrow(AI) - length(pureVec))
  for (i in 1:ncol(AI)) {
    groupi <- which(AI[ ,i] != 0)
    SigmaiJ <- as.matrix(SigmaJ[groupi, ])
    SigmaTJ[i, ] <- apply(SigmaiJ, 2, mean) # Average columns along the rows.
  }
  return(SigmaTJ)
}
