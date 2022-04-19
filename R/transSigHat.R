#' Transform \eqn{\hat{\Sigma}} NOT USED ANYMORE
#'
#' Create the transformed version of \eqn{\hat{\Sigma}} by replaced rows/columns of \eqn{\hat{\Sigma}}, the
#' sample correlation matrix of the data, \eqn{x}, with a specified value, \eqn{r}
#'
#' @param x data matrix of dimensions \eqn{n \times p}
#' @param imp vector of indices of important features
#' @param r replacement value
#' @return a transformed version of \eqn{\hat{\Sigma}}, the expert knowledge matrix for Essential Regression
#' @export

transSigHat <- function(x, imp, r) {
  p <- ncol(x)
  n <- nrow(x)
  sig_hat <- crossprod(x) / n
  for (i in 1:length(imp)) {
    ind <- imp[i]
    sig_hat[, ind] <- r
    sig_hat[ind, ] <- r
  }
  diag(sig_hat) <- 1
  return (sig_hat)
}
