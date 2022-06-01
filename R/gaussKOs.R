#' Gaussian Knockoffs
#'
#' A feature selection method using Model X Knockoffs
#'
#' @importFrom magrittr '%>%'
#' @importFrom foreach '%dopar%'
#' @param y a response vector of dimension \eqn{n}, standardized
#' @param z a matrix of dimensions \eqn{p \times K}
#' @param statistic a statistic to use for variable selection
#' @param fdr the target false discovery rate
#' @param iters the number of iterations to perform
#' @return a vector of selected variable indices
#' @export

gaussKOs <- function(z, y, statistic = knockoff::stat.glmnet_lambdasmax, fdr = 0.05, iters = 20) {
  covar_mat <- t(z) %*% z / dim(z)[1]
  covar_mat <- makePosDef(covar_mat)
  mus <- rep(0, dim(z)[2])
  gauss_kos <- function(x) {
    knockoff::create.gaussian(x, mu = mus, Sigma = covar_mat)
  }
  selected_list <- list()
  knockoff_list <- list()

  selected_list <- foreach::foreach(i = 1:iters) %dopar% {
    result <- knockoff::knockoff.filter(z, y, knockoffs = gauss_kos, statistic = statistic, offset = 0, fdr = fdr)
    result$selected
  }

  len_list <- lapply(selected_list, function(x) { length(x) })
  mm <- which.max(unlist(len_list)) ## use iteration with longest selected_list
  selected_vars <- selected_list[[mm]]
  return(selected_vars)
}
