#' Compute the corresponding p-value for a given correlation, \eqn{r}, and sample size, \eqn{n}.
#'
#' @param r numeric, correlation
#' @param n numeric, sample size
#' @return a list included the correlation, \eqn{r}, the sample size, \eqn{n},
#' the t-test statistic, \eqn{t}, the p-value, \eqn{p}, and the standard error, \eqn{se}
#' @export

corrToP <- function(r, n) {
  #### t-test statistic
  t <- (r * sqrt(n - 2)) / sqrt(1 - r^2)
  #### df = n - 2; 2 tailed p-value
  p <- 2 * (1 - pt(abs(t), (n - 2)))
  #### standard error
  se <- sqrt((1 - r^2) / (n - 2))
  out <- list(r, n, t, p, se)
  names(out) <- c("r", "n", "t", "p", "se")
  return(out)
}
