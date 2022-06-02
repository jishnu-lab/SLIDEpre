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
#' @param priors a vector of important variables
#' @param priors_z a vector of important \eqn{Z}s
#' @param estim the type of estimator used for Bayesian Model Averaging. options include
#' "BMA" = Bayesian Model Averaging Model, "HPM" = Highest Probability Model,
#' "MPM" = Median Probability Model
#' @return a list containing: beta, a vector of \eqn{\beta} estimates;
#' sig_z, the top 10 zs according to posterior probability; BMA_results,
#' an object returned by \code{BAS::bas.lm()}
#' @export

betaBMA <- function(x, y, er_res, priors, priors_z, estim = "HPM") {
  ## scale y
  scale_y <- scale(y, T, T)
  ## scale x
  scale_x <- scale(x, T, T)
  ## get Zs
  prior_z <- predZ(x = scale_x, er_res = er_res)

  ## Step 1: Important Clusters
  ## create initial probabilities vector
  z_imp <- prior_z[, priors_z] %>%
    as.matrix()
  loadings <- er_res$A[, priors_z] %>%
    as.matrix()
  z_imp_probs <- rep(0, ncol(z_imp))

  ## make initial inclusion probabilities
  if (!is.null(priors)) {
    for (i in 1:ncol(loadings)) {
      column <- loadings[, i] ## get one column of loadings matrix A
      abs_col <- abs(column) ## take absolute value
      abs_col <- abs_col[priors] ## get just rows of important features
      num_nonzero <- length(which(abs_col > 0)) ## get number of nonzero entries
      z_imp_probs[i] <- num_nonzero / length(which(column != 0)) ## divide by total number of nonzero entries
    }
    ## rescale probabilities to be [0.5, 1.0]
    z_imp_probs <- scales::rescale(z_imp_probs, to = c(0.5, 1.0))
  } else {
    z_imp_probs <- rep(0.5, ncol(z_imp))
  }

  ## do BMA
  imp_betas_bas <- BAS::bas.lm(scale_y ~ z_imp,
                               prior = "g-prior",
                               alpha = 1,
                               modelprior = BAS::beta.binomial(alpha = 1, beta = 1),
                               force.heredity = FALSE,
                               pivot = TRUE,
                               method = "MCMC",
                               MCMC.iterations = 100000,
                               thin = 1000,
                               initprobs = z_imp_probs)

  imp_betas <- coef(imp_betas_bas, estimator = estim)
  imp_betas <- imp_betas$postmean ## get posterior means as coefficient estimates
  sig_betas <- sort(imp_betas_bas$probne0, index.return = TRUE, decreasing = TRUE) ## sort posterior probs
  sig_betas <- sig_betas$ix[1:11] ## select top 11 by posterior probs
  if (1 %in% sig_betas) { ## if intercept is selected, get rid of it
    sig_betas <- sig_betas[-which(sig_betas == 1)]
  } else {
    sig_betas <- sig_betas[1:10] ## if intercept is not selected, get rid of lowest prob
  }
  sig_zs <- priors_z[(sig_betas - 1)] ## get indices of Zs corresponding to top 10 post probs (adjust for intercept index)

  return(list("beta" = imp_betas, ## beta estimates
              "sig_beta" = sig_betas, ## significant betas by post probs
              "sig_z" = sig_zs, ## significant Zs by post probs
              "BMA_results" = imp_betas_bas)) ## bma results from BAS::bas.lm()
}
