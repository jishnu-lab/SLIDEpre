#' Additional prediction function using the support of \eqn{\beta}.
#'
#' @param X data matrix of dimensions \eqn{n \times p}
#' @param A_hat a matrix of dimensions \eqn{p \times K}
#' @param C_hat a matrix of dimensions
#' @param est_Gamma a matrix of dimensions
#' @param S_beta ????
#' @return ????

Pred_Z_BLP <- function(X, A_hat, C_hat, est_Gamma, S_beta) {
  est_Gamma_inv <- diag(est_Gamma ** (-1))
  G_hat <- crossprod(A_hat, est_Gamma_inv) %*% A_hat + solve(C_hat)
  Z_hat <- X %*% est_Gamma_inv %*% A_hat %*% ginv(G_hat)
  return (Z_hat[,S_beta,drop = F])
}
