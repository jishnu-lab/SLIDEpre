#' Estimate the values of the coefficients \eqn{\beta} with a Dantzig-type estimator.
#'
#' @param y response vector of dimension \eqn{p}
#' @param x data matrix of dimensions \eqn{n \times p}
#' @param A_hat a matrix of dimensions \eqn{p \times K}
#' @param C_hat estimated matrix
#' @param I_hat estimated matrix
#' @param delta numeric value for \eqn{\delta}
#' @param mu ???? \eqn{\mu}
#' @param lambda numeric value for \eqn{\lambda} used in soft thresholding
#' @return Dantzig-type estimator results for \eqn{\beta}

estBetaDant <- function(y, x, A_hat, C_hat, I_hat, delta, mu = 0.5, lambda = 0.5) {
  n <- nrow(x); p <- ncol(x)
  AI <- A_hat[I_hat, ]
  #### h_hat (supplement 2.4)
  h <- solve(crossprod(AI), t(AI) %*% crossprod(x[, I_hat], y) / n)

  C <- C_hat
  mu <- mu * delta
  lbd <- lbd * delta

  #### setting up linear program
  K <- nrow(C)
  c_vec <- c(1, rep(0, 2 * K))
  A_mat <- -c_vec
  A_mat <- rbind(A_mat, c(-1, rep(1, 2 * K)))
  tmp_constr <- C %x% t(c(1, -1))
  A_mat <- rbind(A_mat, cbind(-1 * lambda, rbind(tmp_constr, -tmp_constr)))
  b_vec <- c(0, 0, mu + h, mu - h)
  lp_result <- linprog::solveLP(c_vec, b_vec, A_mat, lpSolve = T)$solution

  while (length(lp_result) == 0) {
    cat("The penalty lambda =", lambda, "is too small and increased by 0.01...\n")
    lambda <- lambda + 0.01
    A_mat[-(1:2), 1] <- lambda
    lp_result <- linprog::solveLP(c_vec, b_vec, A_mat, lpSolve = T)$solution[-1]
  }
  ind <- seq(2, 2 * K, 2)
  sol <- lp_result[ind] - lp_result[ind + 1]

  return (as.vector(sol))
}
