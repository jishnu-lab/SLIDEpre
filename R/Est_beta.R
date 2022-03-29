#' Estimate the values of the coefficients \eqn{\beta} as well as confidence intervals for the estimates
#' if requested.
#'
#' @param Y vector of responses of dimension \eqn{p}
#' @param X data matrix of dimensions \eqn{n \times p}
#' @param sigma sample correlation matrix of the data, \eqn{X}
#' @param A_hat estimated matrix of dimensions \eqn{p \times K}
#' @param C_hat estimated matrix
#' @param Gamma_hat estimated matrix
#' @param I_hat estimated matrix
#' @param I_ind_hat estimated matrix
#' @param CI boolean indicated whether to calculate confidence intervals
#' @param alpha_level value for confidence intervals
#' @param correction type of multiple testing correction to perform
#' @param support whether to do calculations on support of WHAT
#' @return a list including the estimates for \eqn{\beta} and the confidence intervals (if requested)

Est_beta <- function(Y, X, sigma, A_hat, C_hat, Gamma_hat, I_hat, I_ind_hat, CI = T,
                     alpha_level = 0.05, correction = NULL, support = NULL) {
  n <- nrow(X); p <- ncol(X)
  R <- matrix(0, nrow = p, ncol = ncol(A_hat))
  Sigma <- sigma #crossprod(X) / n

  BI <- t(solve(crossprod(A_hat[I_hat,]), t(A_hat[I_hat,])))
  R[I_hat, ] = (Sigma[I_hat, I_hat] - diag(Gamma_hat[I_hat])) %*% BI
  R[-I_hat, ] =  Sigma[-I_hat, I_hat] %*% BI

  if (is.null(support)) {
    beta_est <- solve(crossprod(R), t(R) %*% crossprod(X, Y) / n)
  } else {
    beta_est <- rep(0, ncol(A_hat))
    beta_est[support] <- solve(crossprod(R[, support]), t(R[,support]) %*% crossprod(X, Y) / n)
  }

  sigma2_hat <- Est_sigma2(Y, t(BI) %*% crossprod(X[,I_hat], Y) / n, beta_est, C_hat)
  Omega_hat <- solve(C_hat)

  if (CI) {
    ### (1 - alpha / 2) CI of beta (for each individual component)
    beta_var <- Comp_ASV(sigma2_hat, BI, R, Gamma_hat, beta_est, Omega_hat, I_hat, I_ind_hat)
    beta_var[beta_var < 0] = 0
    if (!is.null(correction)) {
      if (correction == "Bonferroni") {
        alpha_level <- alpha_level / ncol(A_hat)
      }
    }
    CIs_beta <- cbind(lower = beta_est - qnorm(1 - alpha_level / 2) * sqrt(beta_var) / sqrt(n),
                      upper = beta_est + qnorm(1 - alpha_level / 2) * sqrt(beta_var) / sqrt(n))
  } else {
    CIs_beta <- beta_var <- NULL
  }

  return(list(beta = beta_est, CIs = CIs_beta, beta_var = beta_var))
}
