#' Find the median correlation of a cluster by using the nodes in cluster as a list of
#' rows and searching over these rows for the median value in the sample correlation matrix.
#'
#' @param sigma a correlation matrix of dimensions \eqn{p \times p}
#' @param cluster a cluster found by Essential Regression (list format)
#' @param clust_pure a vector of the pure nodes in the cluster
#' @return the median correlation found out of all of the rows provided


findClustMed <- function(sigma, cluster, clust_pure) {
  clust <- unlist(cluster)
  sub_sigma <- sigma
  ## set diagonal NA so we can remove them for median calculation
  diag(sub_sigma) <- NA
  ## subset the correlation matrix to just the nodes in the provided cluster
  sub_sigma <- sub_sigma[clust, clust_pure]
  sub_sigma <- as.vector(sub_sigma) %>% na.omit()
  med_corr <- median(sub_sigma)
  return (med_corr)
}
