#' Prior Information Essential Regression
#'
#' Run Essential Regression with prior information.
#'
#' @importFrom magrittr '%>%'
#' @param y a response vector of dimension \eqn{n}
#' @param x a data matrix of dimensions \eqn{n \times p}
#' @param sigma a sample correlation matrix of dimensions \eqn{p \times p}
#' @param imps a vector of important feature indices/names
#' @param delta \eqn{\delta}, a numerical constant used for thresholding
#' @param thresh_fdr a numerical constant used for thresholding the correlation matrix to
#' control the false discovery rate, default is 0.2
#' @param beta_est a string indicating the type of estimation to use for \eqn{\beta}
#' @param conf_int a boolean indicating whether to calculate confidence intervals for the \eqn{\beta} estimates
#' @param pred a boolean indicating whether to do prediction
#' @param lambda \eqn{\lambda}, a numerical constant used in thresholding
#' @param rep_cv number of replicates for cross-validation
#' @param diagonal a boolean indicating the diagonal structure of the data ???
#' @param merge a boolean indicating the merge type
#' @param equal_var a boolean indicating whether there is equal variance ??
#' @param alpha_level \eqn{\alpha}, a numerical constant used in confidence interval calculation
#' @param thresh a numerical constant used as the threshold for convergence in Variational Bayes
#' @param support a boolean ???
#' @param correction a boolean flag indicating whether to perform Bonferroni multiple testing correction
#' @param change_all a boolean indicating whether to change all entries in \eqn{\hat{\Sigma}}
#' for an important feature (T) or to just change to more extreme values (F)
#' @param verbose a boolean indicating whether to include printing
#' @return a list of results from the Essential Regression framework including: \eqn{K = } number of clusters,
#' \eqn{\hat{A}}, \eqn{\hat{C}}, \eqn{\hat{I}}, the indices of the pure variables, \eqn{\hat{\Gamma}},
#' \eqn{\hat{\beta}}, \eqn{\alpha}-level confidence intervals (if requested), prediction results (if requested),
#' the optimal value of \eqn{\lambda} determined by cross-validation, the optimal value of \eqn{\delta}
#' determined by cross-validation, \eqn{Q}, and the variances of \eqn{\hat{\beta}}
#' @export

priorER <- function(y, x, imps, sigma = NULL, delta, thresh_fdr = 0.2, beta_est = "NULL",
                    conf_int = F, pred = T, lambda = 0.1, rep_cv = 50, diagonal = F,
                    merge = F, equal_var = F, alpha_level = 0.05, thresh = 0.001,
                    support = NULL, correction = T, change_all = F, verbose = F) {
  #### run plainER() first
  plain_er <- plainER(y = y,
                      x = x,
                      sigma = NULL,
                      thresh_fdr = thresh_fdr,
                      delta = delta,
                      beta_est = beta_est,
                      conf_int = conf_int,
                      pred = pred,
                      lambda = lambda,
                      rep_cv = rep_cv,
                      diagonal = diagonal,
                      merge = merge,
                      equal_var = equal_var,
                      alpha_level = alpha_level,
                      support = support,
                      correction = correction,
                      verbose = verbose)
  opt_delta <- plain_er$opt_delta
  opt_lambda <- plain_er$opt_lambda
  sigma <- plain_er$thresh_sigma
  plain_betas <- sigBetas(betas = plain_er$beta, cutoff = alpha_level * 2)

  #### at this point, we have fully run LOVE/done cross-validation and can now
  #### begin to incorporate the prior information
  #### adjust sigma matrix to reflect prior information found in imps
  feats <- readER(plain_er)
  incl_feats <- unlist(feats$clusters) %>% unique()
  excl_feats <- setdiff(seq(1, ncol(sigma)), incl_feats)
  imp_inds <- imps
  if (typeof(imps) == "character") {
    imp_inds <- indName(imps, colnames(x), to_ind = T)
  }
  excl_imp <- intersect(imp_inds, excl_feats)

  #### Essential Regression with Prior Information Part I ######################
  if (length(excl_imp) > 0) {
    #### construct the expert knowledge matrix
    prior_sigma <- makeDelta(x = x,
                             sigma = sigma,
                             imps = excl_imp,
                             er_res = plain_er,
                             change_all = change_all,
                             equal_var = equal_var)
    #### balance the prior_sigma with the sample correlation matrix
    bal_sigma <- findAlpha(x = x,
                           Delta = prior_sigma,
                           alpha_range = seq(0.01, 0.99, 0.01))
    #### run plainER with sigma = bal_sigma
    prior_er <- plainER(y = y,
                        x = x,
                        sigma = bal_sigma$adj_mat,
                        delta = opt_delta,
                        thresh_fdr = NULL,
                        beta_est = beta_est,
                        conf_int = conf_int,
                        pred = pred,
                        lambda = opt_lambda,
                        rep_cv = rep_cv,
                        diagonal = diagonal,
                        merge = merge,
                        equal_var = equal_var,
                        alpha_level = alpha_level,
                        support = support,
                        correction = correction,
                        verbose = verbose)
  } else { ## if all of the important features were in clusters, then don't need to do Part I
    prior_er <- plain_er
  }

  #### Essential Regression with Prior Information Part II #####################
  #### find significant betas
  betas <- sigBetas(betas = prior_er$beta, cutoff = alpha_level * 2)
  sig_betas <- unlist(c(betas$pos_sig, betas$neg_sig))

  #### find important features in each cluster
  clust_feats <- list()
  imp_clusts <- NULL
  prior_er_res <- readER(prior_er)
  for (i in 1:prior_er$K) {
    cluster <- unlist(prior_er_res$clusters[[i]])
    imp_feats_cluster <- intersect(cluster, imps)
    clust_feats[[length(clust_feats) + 1]] <- imp_feats_cluster
    if (length(imp_feats_cluster) > 0) {
      imp_clusts <- c(imp_clusts, i)
    }
  }

  #### attempt to make the clusters with important features significant
  prior_z <- predZ(x = x, er_res = prior_er)
  new_betas <- bayesVarSel(z = prior_z, y = y, imp_clusts = imp_clusts, er_res = prior_er, thresh = thresh)
  new_betas <- sigBetas(betas = new_betas, cutoff = alpha_level * 2)

  #### Essential Regression with Prior Information Part III ####################
  #### compile results
  results <- list("plainER_results" = plain_er, ## plainER - no prior info
                  "priorER_results" = prior_er, ## priorER - prior info
                  "sample_sigma" = sigma, ## sigma used in plainER
                  "prior_knowledge_sigma" = prior_sigma, ## Delta used to represent prior knowledge
                  "priorER_sigma" = bal_sigma, ## sigma used in priorER
                  "plainER_betas" = plain_betas, ## plainER betas
                  "priorER_betas" = betas, ## priorER betas - no bayes
                  "prior_adj_betas" = new_betas) ## variational bayes adjusted betas

  return (results)
}
