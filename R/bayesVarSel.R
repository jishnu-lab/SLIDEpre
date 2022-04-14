#' An implementation of the variational algorithm for Bayesian
#' variable selection paper by Huang, Wang, Liang (2016) \link{https://arxiv.org/pdf/1602.07640.pdf}.
#'
#' @param er_res a list output by the function \code{\link{plainER}} or \code{\link{priorER}}
#' @param y a vector of responses
#' @param alpha a numeric constant for the p-value threshold
#' @return a list including the clusters, pure variables, mixed variables, the significant
#' cluster indices, and the features found in the significant clusters
#' found by ER

bayesVarSel <- function(x, y, imps, er_res, thresh = 0.01) {
  z <- predZ(x, er_res)
  #### get data dimensions
  n <- nrow(z)

  #### get number of clusters
  K <- er_res$K

  #### get which clusters have important variables
  er_feats <- readER(er_res = er_res)
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
  phis <- rep(0.5, K)
  phis[imp_clusts] <- 0.99
  mus <- betas
  sig_sqs <- beta_vars
  sig_hat_sq <- 0.1
  theta_hat <- 0.1

  #### for now, fix v1, a0, b0, nu, lambda = 1
  v1 <- 0.5
  a0 <- 1
  b0 <- 1
  nu <- 1
  lambda <- 1

  #### calculate initial change in betas
  old_phis <- phis
  old_mus <- mus
  max_change <- Inf

  while (max_change >= threshold) {
    for (j in 1:K) {
      #### update mu_j
      phi_x_mu <- phis[-j] * mus[-j]
      mu_sum <- sweep(z[, -j], 2, phi_x_mu, "*")
      mu_sum <- rowSums(mu_sum)
      mus[j] <- (t(y - mu_sum) %*% z[, j]) / (n + (1 / v1))

      #### update sig_sq_j
      sig_sqs[j] <- (sig_hat_sq) / (n + (1 / v1))

      #### update phi_j
      phis[j] <- gtools::inv.logit(log(theta_hat / (1 - theta_hat)) - 0.5 * log((v1 * sig_hat_sq) / sig_sqs[j]) + (mus[j] * mus[j]) / (2 * sig_sqs[j]))
    }

    #### update theta_hat
    theta_hat <- (sum(phis) + K * a0 - K) / (K + a0 + b0 - 2)

    #### update sigma_hat_sq
    beta_bar <- mus * phis
    l2_sq <- (norm((y - z %*% beta_bar), type = "2"))^2
    summand <- (n * (1 - phis) + 1 / v1) * phis * mus * mus
    summand <- summand + (n + 1 / v1) * phis * sig_sqs
    numerator <- l2_sq + sum(summand) + nu * lambda
    denominator <- n + prod((phis + nu + 2))
    sig_hat_sq <- numerator / denominator

    #### calculate change in betas
    changes <- old_phis * old_mus - phis * mus
    max_change <- max(changes)

    #### update old phis, mus
    old_phis <- phis
    old_mus <- mus
  }

  return (phis * mus)
}
