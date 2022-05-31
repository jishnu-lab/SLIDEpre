#' Variational Bayes Variable Selection
#'
#' An implementation of the variational algorithm for Bayesian
#' variable selection paper by Huang, Wang, Liang (2016).
#'
#' @param x a data matrix of dimensions \eqn{n \times p}
#' @param y a response vector of dimension \eqn{n}
#' @param prior_clusts a vector of cluster indices to keep
#' @param er_res a list output by the function link[EssReg]{plainER} or \link[EssReg]{priorER}
#' @param thresh a numeric constant used in determining convergence
#' @return a list including the clusters, pure variables, mixed variables, the significant
#' cluster indices, and the features found in the significant clusters
#' found by ER
#' @export

bayesVarSel <- function(z, y, prior_clusts, er_res, thresh = 0.01) {
  n <- length(y)
  z <- scale(z, center = T, scale = T)
  #### get number of clusters
  K <- er_res$K

  #### get which clusters have important variables
  er_feats <- readER(er_res = er_res)

  #### get betas and their variances from ER
  betas <- er_res$beta
  beta_vars <- er_res$beta_var

  #### initialize vectors of phis, mus, sig_sqs, sigma_hat_sq
  phis <- rep(0.5, K)
  phis[prior_clusts] <- 1
  mus <- betas / phis
  sig_sqs <- beta_vars
  sig_hat_sq <- 0.1
  theta_hat <- 0.1

  #### for now, fix v1, a0, b0, nu, lambda = 1
  v1 <- 1
  a0 <- 0.5
  b0 <- 0.5
  nu <- 1
  lambda <- 1

  #### calculate initial change in betas
  old_phis <- phis
  old_mus <- mus
  max_change <- Inf

  while (max_change > thresh) {
    for (j in 1:K) {
      #### update mu_j
      beta <- phis * mus
      beta_j <- phis[j] * mus[j]
      mu_temp <- z %*% beta - z[, j] * beta_j
      mus[j] <- (t(y - mu_temp) %*% z[, j]) / (n + (1 / v1))

      #### update sig_sq_j
      sig_sqs[j] <- (sig_hat_sq) / (n + (1 / v1))

      #### update phi_j
      if (theta_hat <= 0) {
        theta_hat <- 1e-10
      }
      phis[j] <- gtools::inv.logit(log(theta_hat / (1 - theta_hat)) -
                                     0.5 * log((v1 * sig_hat_sq) / sig_sqs[j]) +
                                     (mus[j] * mus[j]) / (2 * sig_sqs[j]))
    }

    #### update theta_hat
    theta_hat <- (sum(phis) + a0 - 1) / (K + a0 + b0 - 2)
    #### update sigma_hat_sq
    beta_bar <- mus * phis
    l2_sq <- (norm((y - z %*% beta_bar), type = "2"))^2
    summand <- (n * (1 - phis) + 1 / v1) * phis * mus * mus
    summand <- summand + (n + 1 / v1) * phis * sig_sqs
    numerator <- l2_sq + sum(summand) + nu * lambda
    denominator <- n + prod(phis) + nu + 2
    sig_hat_sq <- numerator / denominator

    #### calculate change in betas
    changes <- old_phis * old_mus - phis * mus
    max_change <- max(abs(changes))

    #### update old phis, mus
    old_phis <- phis
    old_mus <- mus
  }

  return (phis * mus)
}
