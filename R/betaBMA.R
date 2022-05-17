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

betaBMA <- function(x, y, er_res, imps, method = "split", estim = "HPM") {
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
  ## make initial inclusion probabilities
  for (i in 1:ncol(loadings)) {
    column <- loadings[, i] ## get one column of loadings matrix A
    abs_col <- abs(column) ## take absolute value
    abs_col <- abs_col[imps] ## get just rows of important features
    num_nonzero <- length(which(abs_col > 0)) ## get number of nonzero entries
    z_imp_probs[i] <- num_nonzero / length(which(column != 0)) ## divide by total number of nonzero entries
  }
  ## rescale probabilities to be [0.5, 1.0]
  z_imp_probs <- scales::rescale(z_imp_probs, to = c(0.5, 1.0))

  ## split zs by correlation similarity
  z_corr <- cor(z_imp)
  corr_dist <- stats::as.dist(1 - z_corr)
  corr_tree <- stats::hclust(corr_dist, method = "complete")
  corr_dend <- stats::as.dendrogram(corr_tree)

  ## find location of cutting
  k <- 1
  k_unknown <- TRUE
  while (k_unknown) {
    k <- k + 1
    clusters <- dendextend::cutree(corr_dend, k = k)
    cl_tab <- table(clusters)
    biggest_group <- 25
    for (i in 1:length(cl_tab)) {
      if (cl_tab[i] > biggest_group) {
        biggest_group <- cl_tab[i]
      }
    }
    if (biggest_group <= 25) {
      k_unknown <- FALSE
    }
  }

  cat("Running BAS.lm . . . \n")
  all_imp_betas <- rep(0, ncol(z_imp))
  for (i in 1:k) { ## loop through dendrogram clusters
    cat("    Running BMA for Dendrogram Cluster #", i, ". . . \n")
    ## get one dendrogram cluster
    cluster <- clusters[which(clusters == i)]
    cluster_inds <- names(cluster) %>%
      as.numeric()
    z_clust <- z_imp[, cluster_inds]
    z_clust_imp_probs <- z_imp_probs[cluster_inds]

    ## do BMA
    imp_betas_bas <- BAS::bas.lm(scale_y ~ z_clust,
                                 prior = "g-prior",
                                 alpha = 1,
                                 modelprior = BAS::beta.binomial(alpha = 1, beta = 1),
                                 force.heredity = FALSE,
                                 pivot = TRUE,
                                 method = "BAS",
                                 initprobs = z_clust_imp_probs)

    imp_betas <- coef(imp_betas_bas, estimator = estim)
    imp_betas <- imp_betas$postmean ## get posterior means
    imp_betas <- imp_betas[-1] ## remove intercept estimate

    all_imp_betas[which(clusters == i)] <- imp_betas
  }

  ## Step 2: Non-Important Clusters
  new_y <- scale_y - z_imp %*% all_imp_betas
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
  all_betas[imp_clusts] <- unlist(all_imp_betas)
  all_betas[-imp_clusts] <- unlist(nonimp_beta)
  return(list("beta_est" = all_betas,
              "imp_clusters" = imp_clusts))
}
