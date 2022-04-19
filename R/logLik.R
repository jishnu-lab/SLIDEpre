#' Log-Likelihood of Data
#'
#' Evaluate the marginal log-likelihood of the data, \eqn{x}, given weight \eqn{\alpha} and expert matrix \eqn{\Delta}:
#' ADD EQUATION
#'
#' @param alpha weight for the intensity of the balance between the expert knowledge matrix, \eqn{\Delta}, and the data
#' @param x data matrix of dimensions \eqn{n \times p}
#' @param Delta the correlation matrix of dimensions \eqn{p \times p} created with regard to expert knowledge
#' @return numeric
#' @export

logLik <- function(alpha, Delta, x) {
  n <- nrow(x)
  p <- ncol(x)
  xtx <- crossprod(x)
  run_loglik <- -(n * p * 0.5) * log(pi)
  phi <- toPhi(alpha, n, p)
  val1 <- (phi + n) * 0.5
  val2 <- phi * 0.5
  run_loglik <- run_loglik + logGamma(p, val1, val2)
  eig_Delta <- eigen(Delta)
  eig_Delta <- eig_Delta$values
  ## determinant = product of eig vals
  log_det_Delta <- sum(log(eig_Delta))
  addition <- 0.5 * phi * ((p * log(phi - p - 1)) + log_det_Delta)
  run_loglik <- run_loglik + addition
  ## lower determinant
  phip1 <- phi - p - 1
  low_det <- (phip1 * Delta) + xtx
  low_det_eigs <- eigen(low_det)
  low_det_eigs <- low_det_eigs$values
  log_low_det <- sum(log(low_det_eigs))
  addition <- (0.5 * (phi + n)) * log_low_det
  run_loglik <- run_loglik - addition
  return (as.numeric(run_loglik))
}
