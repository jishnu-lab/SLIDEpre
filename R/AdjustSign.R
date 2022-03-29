#' Perform a sign operation on matrix \eqn{\Sigma} according to the sign of \eqn{A_I}.
#'
#' @param Sigma a sample correlation matrix of dimensions \eqn{p \times p}
#' @param AI a matrix of dimensions \eqn{p \times K}
#' @return a matrix of dimensions \eqn{p \times p}

AdjustSign <- function(Sigma, AI) {
  SigmaS <- matrix(0, nrow(AI), nrow(AI))
  for (i in 1:nrow(AI)) {
    index <- which(AI[i, ] != 0)
    if (length(index) != 0) {
      SigmaS[i, ] = sign(AI[i,index]) * Sigma[i, ]
    }
  }
  return(SigmaS)
}
