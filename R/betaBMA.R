#' PriorER - Bayesian Model Averaging For \eqn{\beta}s
#'
#' Find new estimates for \eqn{\beta}s using Bayesian Adaptive
#' Sampling for Bayesian Model Averaging as in the \code{BAS} package by Clyde et al.
#'
#' @importFrom magrittr '%>%'
#' @importFrom foreach '%dopar%'
#' @param x a data matrix of dimensions \eqn{n \times p}
#' @param y a response vector of dimension \eqn{n}
#' @param er_res the result of \link{plainER} or \link{priorER}
#' @param imps a vector of important variables
#' @return nothing is returned, saves boxplot of cross-validation results for user to use
#' in selecting optimal \eqn{\delta}
#' @export

betaBMA <- function(x, y, er_res, imps) {
  #### find important features in each cluster
  er_read <- readER(er_res)
  clust_feats <- list()
  imp_clusts <- NULL
  for (i in 1:er_res$K) {
    cluster <- unlist(er_read$clusters[[i]])
    imp_feats_cluster <- intersect(cluster, imps)
    clust_feats[[length(clust_feats) + 1]] <- imp_feats_cluster
    if (length(imp_feats_cluster) > 0) {
      imp_clusts <- c(imp_clusts, i)
    }
  }

  ## get Zs
  prior_z <- predZ(x = scale(x, T, T), er_res = er_res)
  ## scale y
  scale_y <- scale(y, T, T)
  ## create initial probabilities vector
  ## (0.5 for uncertain clusters, 1 for important clusters to force inclusion)
  imp_probs <- rep(0.5, ncol(prior_z))
  imp_probs[imp_clusts] <- 1


  ## do BAS for BMA
  cat("Running bas.lm . . .")
  new_betas <- BAS::bas.lm(scale_y ~ prior_z,
                           prior = "g-prior",
                           alpha = 1,
                           modelprior = BAS::beta.binomial(),
                           force.heredity = FALSE,
                           pivot = TRUE,
                           method = "MCMC",
                           include.always = '~ ',
                           MCMC.iterations = 100000)

  cat("Finding median probability model . . .")
  mpm <- coef(new_betas, estimator = "MPM")
  cat("Finding highest probability model . . . ")
  hpm <- coef(new_betas, estimator = "HPM")
  cat("Finding Bayesian model average model . . .")
  bma <- coef(new_betas, estimator = "BMA")

  rm(new_betas)
  gc()

  return (list("MPM" = mpm$postmean,
               "HPM" = hpm$postmean,
               "BMA" = bma$postmean))
}
