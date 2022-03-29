#' Recover the estimated submatrix \eqn{A_I} by given the pure node group.
#'
#' @param estGroupList a list of group indices of the pure nodes with sign
#' @param p the number of features
#' @return a matrix of dimensions \eqn{p \times K}

RecoverAI <- function(estGroupList, p) {
  K <- length(estGroupList)
  A <- matrix(0, p, K)
  for (i in 1:K) {
    groupi <- estGroupList[[i]]
    A[groupi[[1]],i] <- 1
    groupi2 <- groupi[[2]]
    if (length(groupi2) != 0)
      A[groupi2, i] <- -1
  }
  return(A)
}
