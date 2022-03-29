#' An additional prediction function that uses the support of \eqn{\beta}.
#'
#' @param X a data matrix of dimensions \eqn{n \times p}
#' @param A_hat a matrix of dimensions \eqn{p \times K}
#' @param C_hat a matrix of dimensions
#' @param S_beta ???
#' @return ????

Pred_Z_BLP_avg <- function(X, A_hat, C_hat, S_beta) {
  G_hat <- crossprod(A_hat) + solve(C_hat)
  Z_hat <- X %*% A_hat %*% MASS::ginv(G_hat)
  return (Z_hat[,S_beta,drop = F])
}
