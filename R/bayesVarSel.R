#' An implementation of the variational algorithm for Bayesian
#' variable selection paper by Huang, Wang, Liang (2016) \link{https://arxiv.org/pdf/1602.07640.pdf}.
#'
#' @param er_res a list output by the function \code{\link{plainER}} or \code{\link{priorER}}
#' @param y a vector of responses
#' @param alpha a numeric constant for the p-value threshold
#' @return a list including the clusters, pure variables, mixed variables, the significant
#' cluster indices, and the features found in the significant clusters
#' found by ER

bayesVarSel <- function(y, x, imps, er_res) {
  #### get data dimensions
  n <- nrow(x)
  p <- ncol(x)

  #### get number of clusters
  K <- er_res$K

  #### get which clusters have important variables
  er_feats <- readER(er_res = er_res, y = y)
  imp_clusts <- c()
  for (i in 1:length(er_feats$clusters)) {
    clust <- unlist(er_feats$clusters[[i]])
    has_imp <- intersect(imps, clust)
    if (length(has_imp) > 0) {
      imp_clusts <- c(imp_clusts, i)
    }
  }

  #### get betas and their variances from ER
  betas <- er_res$beta
  beta_vars <- er_res$beta_var

  #### initialize vectors of phis, mus, sig_sqs, sigma_hat_sq
  phis <- rep(0, K)
  phis[imp_clusts] <- 1
  mus <- betas
  sig_sqs <- beta_vars
  sigma_hat_sq <- 1
  beta_hat <- 1

  #### calculate initial maximum entropy of phis


  while (change >= threshold) {
    for (j in 1:K) {
      #### for now, fix v1 = 1
      v1 <- 1
      #### update mu_j
      phi_x_mu <- phis[-j] * mus[-j]
      mu_sum <- x[, -j] * phi_x_mu
      mus[j] <- (t(y - mu_sum) %*% x[, j]) / (n + (1 / v1))

      #### update sig_sq_j
      sig_sqs[j] <- (sigma_hat_sq) / (n + (1 / v1))

      #### update phi_j
      phis[j] <- gtools::invlogit(log(beta_hat / (1 - beta_hat)) - 0.5 * log((v1 * sigma_hat_sq) / sigma_hat_sq[j]) + (mu[j] * mu[j]) / (2 * sig_hat_sq[j]))
    }
  }







}
