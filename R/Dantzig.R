#' I am not really sure what this does.
#'
#' @param C_hat a matrix
#' @param y response vector of dimension ???
#' @param lbd \eqn{\lambda}, a positive constant
#' @return ???

Dantzig <- function(C_hat, y, lbd) {
  K <- length(y) ## number of clusters (columns in Y_hat)
  cvec <- rep(1, 2 * K) ## just 1s because we want to minimize x itself
  print(cvec)
  bvec <- c(lbd + y, lbd - y, rep(0, 2 * K))
  print(bvec)
  new_C_hat <- matrix(0, K, 2 * K)
  for (i in 1:K) {
    new_C_hat[i, ] <- c(C_hat[i,], -C_hat[i,]) ## two copies of C_hat, but second copy is opposite signed
  }
  print(new_C_hat)
  Amat <- rbind(new_C_hat, -new_C_hat, diag(-1, nrow = 2 * K))
  print(Amat)
  ## minimize cvec'x subject to Amatx ≤ bvec and x ≥ 0
  LPsol <- linprog::solveLP(cvec, bvec, Amat, lpSolve = T)$solution
  beta <- LPsol[1:K] - LPsol[(K + 1):(2 * K)] ## split up positive and negative portions of beta
  return(beta)
}


