#' Use the given \eqn{\delta} to calculate the fitted \eqn{A_I}, pure variable
#' indices in list form and vector form. Also return estimated \eqn{Y} and \eqn{C} for
#' the following Dantzig estimation.
#'
#' @param Sigma a correlation matrix of dimensions \eqn{p \times p}
#' @param optDelta \eqn{\delta}, a numerical constant
#' @param se_est ???
#' @param merge a boolean indicating merge style
#' @return a list including: \eqn{A_I}, a matrix of dimensions \eqn{p \times K},
#' a vector of the indices of estimated pure variables, and a list of the indices of
#' the estimated pure variables


EstAI <- function(Sigma, optDelta, se_est, merge) {
  #### get absolute value of covariance matrix
  off_Sigma <- abs(Sigma)

  #### set entries on main diagonal to 0
  diag(off_Sigma) <- 0

  #### calculate the maximal absolute value for each row of the given matrix
  result_Ms <- FindRowMax(off_Sigma)
  Ms <- result_Ms$M #### maximal abs values
  arg_Ms <- result_Ms$arg_M #### first index where max abs values are achieved

  #### estimate list of pure node indices for given Sigma and delta
  #### a node is pure if
  resultPure <- FindPureNode(off_Sigma, optDelta, Ms, arg_Ms, se_est, merge)
  estPureIndices <- resultPure$pureInd
  print(estPureIndices)
  estPureVec <- resultPure$pureVec

  #### Estimate the sign subpartition of pure node sets. If there is an element
  #### of a list is empty, then a empty list will be put in that position
  estSignPureIndices <- FindSignPureNode(estPureIndices, Sigma)

  #### Recover the estimated submatrix A_I given the pure node group.
  AI <- RecoverAI(estSignPureIndices, nrow(off_Sigma))

  return(list(AI = AI, pureVec = estPureVec, pureSignInd = estSignPureIndices))
}
