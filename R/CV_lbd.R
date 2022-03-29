#' Use cross-validation to select \eqn{\lambda} for estimating \eqn{\Omega}.
#' Split the data into two parts and estimate \eqn{C} on both data sets. Then, for each
#' \eqn{lambda}, calculate \eqn{\Omega} on the first dataset and calculate the loss on the second dataset.
#' Find the lambda which minimizes \eqn{C \cdot \Omega - log(|\Omega|)}
#'
#' @param x a data matrix of dimensions \eqn{n \times p}
#' @param lbdGrids a vector of numerical constants over which to search for the optimal \eqn{\lambda}
#' @param AI the estimated matrix \eqn{A_I} of dimensions \eqn{p \times K}
#' @param pureVec a vector of indices of pure variables
#' @param diagonal a boolean indicating the diagonal structure of \eqn{C}
#' @return the selected optimal \eqn{\lambda}

CV_lbd <- function(x, lbdGrids, AI, pureVec, diagonal) {
  sampInd <- sample(nrow(x), floor(nrow(x) / 2))
  X1 <- x[sampInd, ]
  X2 <- x[-sampInd, ]
  Sigma1 <- crossprod(X1) / nrow(X1)
  Sigma2 <- crossprod(X2) / nrow(X2)
  C1 <- EstC(Sigma1, AI, diagonal)
  C2 <- EstC(Sigma2, AI, diagonal)
  loss <- c()
  for (i in 1:length(lbdGrids)) {
    Omega <- estOmega(lbdGrids[i], C1)
    det_Omega <- det(Omega)
    loss[i] <- ifelse(det_Omega <= 0, Inf, sum(Omega * C2) - log(det_Omega))
  }
  return(lbdGrids[which.min(loss)])
}
