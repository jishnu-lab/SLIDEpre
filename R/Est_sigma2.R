#' I am not sure what this does.
#'
#' @param Y a response vector of dimension \eqn{p}
#' @param h Not sure
#' @param beta_est estimated values of \eqn{\beta}
#' @param C_hat estimated matrix
#' @return ???

Est_sigma2 <- function(Y, h, beta_est, C_hat) {
  n <- length(Y)
  sigma2_hat <- crossprod(Y) / n - 2 * t(beta_est) %*% h + t(beta_est) %*% C_hat %*% beta_est
  return (ifelse(sigma2_hat < 0, 0, sigma2_hat))
}
