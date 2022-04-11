#' Estimate the matrix \eqn{A_J} using soft thresholding.
#'
#' @param Omega an estimate of \eqn{C^{-1}} of dimensions \eqn{p \times p}
#' @param Y response matrix of dimensions \eqn{K \times |J|}
#' @param lbd \eqn{\lambda}, the chosen constant for the RHS constraint for soft thresholding
#' @return a matrix of dimensions \eqn{|J| \times K}

EstAJInv <- function(Omega, Y, lbd) {
  AJ <- matrix(0, ncol(Y), nrow(Y))
  for (i in 1:ncol(Y)) {
    Atilde <- Omega %*% as.matrix(Y[ ,i])
    AJ[i, ] <- LP(Atilde, lbd)
    if (sum(abs(AJ[i, ])) > 1)
      AJ[i,] <- AJ[i,] / sum(abs(AJ[i, ]))
  }
  return(AJ)
}
