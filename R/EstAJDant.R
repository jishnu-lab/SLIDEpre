#' Estimate the matrix \eqn{A_J} with a Dantzig-type estimator.
#'
#' @param C_hat a matrix
#' @param Y_hat a matrix
#' @param lbd \eqn{\lambda}, a positive constant
#' @param se_est_J the estimated standard errors of indices in \eqn{J}
#' @return an estimate of matrix \eqn{A_J}

EstAJDant <- function(C_hat, Y_hat, lbd, se_est_J) {
  AJ <- matrix(0, ncol(Y_hat), nrow(Y_hat))
  for (i in 1:ncol(Y_hat)) {
    AJ[i, ] <- Dantzig(C_hat, Y_hat[,i], lbd * se_est_J[i])
    if (sum(abs(AJ[i, ])) > 1)
      AJ[i,] <- AJ[i,] / sum(abs(AJ[i, ]))
  }
  return(AJ)
}
