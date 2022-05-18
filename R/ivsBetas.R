#' IVS \eqn{\beta} Estimation
#'
#' Find new estimates for \eqn{\beta}s using IVS and Essential Regression
#' estimation method.
#'
#' @importFrom magrittr '%>%'
#' @importFrom foreach '%dopar%'
#' @param x a data matrix of dimensions \eqn{n \times p}
#' @param y a response vector of dimension \eqn{n}
#' @param er_res the result of \link{plainER} or \link{priorER}
#' @param imps a vector of important variables
#' @param estim the type of estimator used for Bayesian Model Averaging. options include
#' "BMA" = Bayesian Model Averaging Model, "HPM" = Highest Probability Model,
#' "MPM" = Median Probability Model
#' @return a vector of \eqn{\beta} estimates
#' @export

ivsBetas <- function() {
return()
}


