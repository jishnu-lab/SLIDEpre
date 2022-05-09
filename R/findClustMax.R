#' Find Clusters Maximums
#'
#' Find the absolute minimum of the maximum correlations of each of the target nodes in the
#' given cluster.
#'
#' @param sigma a correlation matrix of dimensions \eqn{p \times p}
#' @param cluster a cluster found by Essential Regression (list format)
#' @param clust_targ a vector of the target nodes in the cluster
#' @return the absolute maximum correlation found out of all of the rows provided
#' @export


findClustMax <- function(sigma, cluster, clust_targ) {
  clust <- unlist(cluster)
  sub_sigma <- sigma
  ## set diagonal to 0 because we never want to select a value on the diagonal
  ## this should never happen anyway, but it is good to be extra cautious
  diag(sub_sigma) <- 0
  ## subset the correlation matrix to just the nodes in the provided cluster
  sub_sigma <- as.data.frame(sub_sigma[clust, clust_targ])
  colnames(sub_sigma) <- clust_targ
  ## correlation can be negative, so find minimum of the absolute value
  sub_sig_abs <- abs(sub_sigma)
  colmax <- apply(sub_sig_abs, 2, max)
  where_absmax <- which(sub_sig_abs == min(colmax), arr.ind=TRUE)
  if (length(where_absmax) > 2) {
    where_absmax <- where_absmax[1, ]
  }
  ## get correlation from not absolute value matrix
  max_corr <- sub_sigma[where_absmax[1], where_absmax[2]]
  return (max_corr)
}
