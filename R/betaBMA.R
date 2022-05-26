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
#' @param imps_z a vector of important \eqn{Z}s
#' @param estim the type of estimator used for Bayesian Model Averaging. options include
#' "BMA" = Bayesian Model Averaging Model, "HPM" = Highest Probability Model,
#' "MPM" = Median Probability Model
#' @return a vector of \eqn{\beta} estimates
#' @export

betaBMA <- function(x, y, er_res, imps, imps_z, estim = "HPM") {
  ## scale y
  scale_y <- scale(y, T, T)
  ## scale x
  scale_x <- scale(x, T, T)
  ## get Zs
  prior_z <- predZ(x = scale_x, er_res = er_res)

  ## Step 1: Important Clusters
  ## create initial probabilities vector
  z_imp <- prior_z[, imps_z] %>%
    as.matrix()
  loadings <- er_res$A[, imps_z] %>%
    as.matrix()
  z_imp_probs <- rep(0, ncol(z_imp))

  ## make initial inclusion probabilities
  if (!is.null(imps)) {
    for (i in 1:ncol(loadings)) {
      column <- loadings[, i] ## get one column of loadings matrix A
      abs_col <- abs(column) ## take absolute value
      abs_col <- abs_col[imps] ## get just rows of important features
      num_nonzero <- length(which(abs_col > 0)) ## get number of nonzero entries
      z_imp_probs[i] <- num_nonzero / length(which(column != 0)) ## divide by total number of nonzero entries
    }
    ## rescale probabilities to be [0.5, 1.0]
    z_imp_probs <- scales::rescale(z_imp_probs, to = c(0.5, 1.0))
  } else {
    z_imp_probs <- rep(0.5, ncol(z_imp))
  }

  ## do BMA
  imp_betas_bas <- BAS::bas.lm(scale_y ~ z_imp, ## INTERCEPT
                               prior = "g-prior",
                               alpha = 1,
                               modelprior = BAS::beta.binomial(alpha = 1, beta = 1),
                               force.heredity = FALSE,
                               pivot = TRUE,
                               method = "BAS",
                               initprobs = z_imp_probs)

  imp_betas <- coef(imp_betas_bas, estimator = estim)
  imp_betas <- imp_betas$postmean ## get posterior means

  ## Step 2: Non-Important Clusters
  z_imp_inter <- cbind(1, z_imp) ## add column for intercept
  new_y <- scale_y - z_imp_inter %*% imp_betas
  new_A <- er_res$A[, -imps_z]
  new_C <- er_res$C[-imps_z, -imps_z]
  new_I_clust <- er_res$I_clust[-imps_z]
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
  all_betas <- rep(0, ncol(as.matrix(prior_z)))
  imp_betas_nointer <- imp_betas[-1]
  all_betas[imps_z] <- unlist(imp_betas_nointer)
  all_betas[-imps_z] <- unlist(nonimp_beta)
  return(c(all_betas, imp_betas[1]))
}
