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
#' @param estim the type of estimator used for Bayesian Model Averaging. options include
#' "BMA" = Bayesian Model Averaging Model, "HPM" = Highest Probability Model,
#' "MPM" = Median Probability Model
#' @return a vector of \eqn{\beta} estimates
#' @export

betaBMA <- function(x, y, er_res, imps, estim = "HPM") {
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

  ## scale y
  scale_y <- scale(y, T, T)
  ## scale x
  scale_x <- scale(x, T, T)
  ## get Zs
  prior_z <- predZ(x = scale_x, er_res = er_res)

  ## Step 1: Important Clusters
  ## create initial probabilities vector
  z_imp <- prior_z[, imp_clusts]
  loadings <- er_res$A[, imp_clusts]
  z_imp_probs <- rep(0, ncol(z_imp))
  ## make initial inclusion probabilities - weighted average of loadings in A matrix
  for (i in 1:ncol(loadings)) {
    column <- loadings[, i]
    abs_col <- abs(column)
    col_sum <- sum(abs(column))
    n_nonzero <- sum(abs_col != 0)
    z_imp_probs[i] <- col_sum / n_nonzero
  }

  p <- ncol(z_imp)
  cat("Running BAS.lm . . . ")
  imp_betas <- BAS::bas.lm(scale_y ~ z_imp,
                           prior = "g-prior",
                           alpha = 1,
                           modelprior = BAS::beta.binomial(alpha = 1, beta = 1),
                           force.heredity = FALSE,
                           pivot = TRUE,
                           method = "MCMC",
                           initprobs = z_imp_probs,
                           MCMC.iterations = p * 2^15,
                           thin = p) ## thin every p iterations

  imp_betas <- coef(imp_betas, estimator = estim)
  imp_betas <- imp_betas$postmean ## get posterior mean
  imp_betas <- imp_betas[-1] ## remove intercept estimate

  ## Step 2: Non-Important Clusters
  new_y <- scale_y - z_imp %*% imp_betas
  new_A <- er_res$A[, -imp_clusts]
  new_C <- er_res$C[-imp_clusts, -imp_clusts]
  new_I_clust <- er_res$I_clust[-imp_clusts]
  new_I <- unlist(new_I_clust)
  ## re-estimate betas for the non-important features
  nonimp_betas <- estBeta(y = new_y,
                          x = scale_x,
                          sigma = cor(scale_x),
                          A_hat = new_A,
                          C_hat = new_C,
                          Gamma_hat = er_res$Gamma,
                          I_hat = new_I,
                          I_hat_list = new_I_clust,
                          conf_int = T,
                          alpha = 0.05,
                          correction = TRUE,
                          support = NULL)
  ## get beta estimates
  nonimp_beta <- nonimp_betas$beta_hat

  ## concatenate the nonimportant and important beta estimates
  all_betas <- rep(0, ncol(prior_z))
  all_betas[imp_clusts] <- unlist(imp_betas)
  all_betas[-imp_clusts] <- unlist(nonimp_beta)
  return(list("beta_est" = all_betas,
              "imp_clusters" = imp_clusts))
}
