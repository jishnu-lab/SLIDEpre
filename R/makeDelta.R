#' Construct the expert knowledge matrix used as the target for the sample correlation matrix
#' reweighting. This is the first step of including prior knowledge in Essential Regression.
#'
#' @param y a response vector of dimension \eqn{n}
#' @param x a data matrix of dimensions \eqn{n \times p}
#' @param imp a vector of the indices of the features determined to be important by expert knowledge
#' @param delta_range a vector of values to search over when picking the optimal \eqn{\delta} by cross-validation
#' @param lbd_range a vector of values to search over when picking the optimal \eqn{\lambda} by cross-validation
#' @return
#' @export

makeDelta <- function(y, x, imp, delta_range, lbd_range, beta_est = "LS", diagonal = F, merge = F,
                      equal_var = F, alpha_level = 0.05, support = F, correction = "Bonferroni") {
  ## must do an initial run of plain Essential Regression in order to get information about
  ## the data structure. use this for running ER-priors
  n <- nrow(x)
  p <- ncol(x)
  samp_corr <- crossprod(x) / n
  er_res <- ER(Y = y, X = x, sigma = samp_corr, delta = delta_range, lbd = lbd_range,
               beta_est = beta_est, CI = T, pred = F, rep_CV = 100, diagonal = diagonal,
               merge = merge, equal_var = equal_var, alpha_level = alpha_level, support = support,
               correction = correction, verbose = F)

  results <- readER(er_res)
  for (i in 1:length(results$clusters)) {
    clust <- results$clusters[[i]]
    pure <- results$pure_vars
    clust_pure <- intersect(unlist(clust), pure)
    clust_min <- findClustMin(samp_corr, clust, clust_pure)
    samp_corr[imp, clust_pure] <- clust_min
    samp_corr[clust_pure, imp] <- clust_min
  }

  Delta <- Matrix::nearPD(samp_corr)
  Delta <- Delta$mat

  return (Delta)
}
