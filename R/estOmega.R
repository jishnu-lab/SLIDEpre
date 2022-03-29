#' For a given \eqn{\lambda} and \eqn{C}, find \eqn{C^{-1}}.
#'
#' @param lbd \eqn{\lambda}
#' @param C a square, symmetric matrix
#' @return \eqn{\Omega = C^{-1}}
#' @export

estOmega <- function(lbd, C) {
  K <- nrow(C)
  omega <- matrix(0, K, K)
  for (i in 1:K) {
    omega[,i] <- solve_row(i, C, lbd)
  }
  return(omega)
}
