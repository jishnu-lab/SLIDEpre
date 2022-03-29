#' Find the absolute minimum correlation of a cluster by using the nodes in cluster as a list of
#' rows and searching over these rows for the minimum value in the sample correlation matrix.
#'
#' @param sigma a correlation matrix of dimensions \eqn{p \times p}
#' @param cluster a cluster found by Essential Regression (list format)
#' @param clust_pure a vector of the pure nodes in the cluster
#' @return the absolute minimum correlation found out of all of the rows provided


findClustMin <- function(sigma, cluster, clust_pure) {
  clust <- unlist(clust)
  sub_sigma <- samp_corr
  ## set diagonal to infinity because we never want to select a value on the diagonal
  ## this should never happen anyway, but it is good to be extra cautious
  diag(sub_sigma) <- Inf
  ## subset the correlation matrix to just the nodes in the provided cluster
  sub_sigma <- sub_sigma[clust, clust_pure]
  ## correlation can be negative, so find minimum of the absolute value
  sub_sig_abs <- abs(sub_sigma)
  min_abs <- which(sub_sig_abs == min(sub_sig_abs), arr.ind=TRUE)
  ## get correlation from not absolute value matrix
  min_corr <- sub_sigma[min_abs[1], min_abs[2]]
  return (min_corr)
}
