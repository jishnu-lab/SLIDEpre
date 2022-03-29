#' Estimate \eqn{C}. If diagonal=True, estimate only diagonal elements.
#'
#' @param Sigma a sample correlation matrix of dimensions \eqn{p \times p}
#' @param AI a matrix of dimensions \eqn{p \times K}
#' @param diagonal boolean indicating whether to only estimated diagonal elements of \eqn{C}
#' @return an estimate of \eqn{C} of dimensions \eqn{K \times K}

EstC <- function(Sigma, AI, diagonal) {
  K <- ncol(AI)
  C <- diag(0, K, K)
  for (i in 1:K) {
    groupi <- which(AI[ ,i] != 0)
    sigmai <- as.matrix(abs(Sigma[groupi,groupi]))
    tmpEntry <- sum(sigmai) - sum(diag(sigmai))
    C[i,i] <- tmpEntry / (length(groupi) * (length(groupi) - 1))
    if (!diagonal && i < K) {
      for (j in (i+1):K) {
        groupj <- which(AI[ ,j]!=0)
        # adjust the sign for each row
        sigmaij <- AI[groupi,i] * as.matrix(Sigma[groupi, groupj])
        sigmaij <- t(AI[groupj, j] * t(sigmaij))
        C[i,j] <- C[j,i] <- sum(sigmaij) / (length(groupi) * length(groupj))
      }
    }
  }
  return(C)
}
