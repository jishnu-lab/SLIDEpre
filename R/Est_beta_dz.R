#' Estimate the values of the coefficients \eqn{\beta} with a Dantzig-type estimator.
#'
#' @param Y response vector of dimension \eqn{p}
#' @param X data matrix of dimensions \eqn{n \times p}
#' @param A_hat a matrix of dimensions \eqn{p \times K}
#' @param C_hat estimated matrix
#' @param I_hat estimated matrix
#' @param delta_opt numeric value for \eqn{\delta}
#' @param mu ???? \eqn{\mu}
#' @param lbd numeric value for \eqn{\lambda} used in soft thresholding
#' @return Dantzig-type estimator results for \eqn{\beta}

Est_beta_dz <- function(Y, X, A_hat, C_hat, I_hat, delta_opt, mu = 0.5, lbd = 0.5) {
  n <- nrow(X); p <- ncol(X)
  AI <- A_hat[I_hat,]
  h <- solve(crossprod(AI), t(AI) %*% crossprod(X[ ,I_hat], Y) / n)

  C <- C_hat
  mu <- mu * delta_opt
  lbd <- lbd * delta_opt

  K <- nrow(C)
  cvec <- c(1, rep(0, 2*K))
  Amat <- -cvec
  Amat <- rbind(Amat, c(-1, rep(1, 2*K)))
  tmp_constr <- C %x% t(c(1,-1))
  Amat <- rbind(Amat, cbind(-1 * lbd, rbind(tmp_constr, -tmp_constr)))
  bvec <- c(0, 0, mu + h, mu - h)

  lpResult <- linprog::solveLP(cvec, bvec, Amat, lpSolve = T)$solution
  while (length(lpResult) == 0) {
    cat("The penalty lambda =", lbd, "is too small and increased by 0.01...\n")
    lbd <- lbd + 0.01
    Amat[-(1:2), 1] <- lbd
    lpResult <- linprog::solveLP(cvec, bvec, Amat, lpSolve = T)$solution[-1]
  }
  ind <- seq(2, 2*K, 2)
  sol <- lpResult[ind] - lpResult[ind + 1]

  return (as.vector(sol))
}
