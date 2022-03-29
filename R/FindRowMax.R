#' Calculate the maximal absolute value for each row of the given matrix.
#'
#' @param Sigma a matrix of dimensions \eqn{p \times p}
#' @return a list including the indices of the maximal absolute values of the rows of
#' \eqn{\Sigma} as well as the values

FindRowMax <- function(Sigma) {
  p <- nrow(Sigma)
  M <- arg_M <- rep(0, p)
  for (i in 1:p) {
    row_i <- Sigma[i,]
    arg_M[i] <- which.max(row_i)
    M[i] <- row_i[arg_M[i]]
  }
  return(list(arg_M = arg_M, M = M))
}
