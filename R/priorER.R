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
#' @param support a boolean ???
#' @param correction a string indicating the type of multi-testing correction to perform
#' @param change_all a boolean indicating whether to change all entries in \eqn{\hat{\Sigma}}
#' for an important feature (T) or to just change to more extreme values (F)
#' @param verbose a boolean indicating whether to include printing
#' @return a list of results from the Essential Regression framework including: \eqn{K = } number of clusters,
#' \eqn{\hat{A}}, \eqn{\hat{C}}, \eqn{\hat{I}}, the indices of the pure variables, \eqn{\hat{\Gamma}},
#' \eqn{\hat{\beta}}, \eqn{\alpha}-level confidence intervals (if requested), prediction results (if requested),
#' the optimal value of \eqn{\lambda} determined by cross-validation, the optimal value of \eqn{\delta}
#' determined by cross-validation, \eqn{Q}, and the variances of \eqn{\hat{\beta}}
#' @export

priorER <- function(y, x, sigma, imps, delta, thresh_fdr = 0.2, beta_est = "NULL",
                    conf_int = F, pred = T, lambda = 0.1, rep_cv = 50, diagonal = F,
                    merge = F, equal_var = F, alpha_level = 0.05, support = NULL,
                    correction = "Bonferroni", change_all = F, verbose = F) {
  #### run plainER() first
  plain_er <- plainER(y = y,
                      x = x,
                      sigma = sigma,
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

  if (length(excl_imp) == 0) {
    print("All Important Features Already In ER")
    return (plain_ER)
  }

  prior_sigma <- makeDelta(y = y,
                           x = x,
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
  prior_er[["prior_sigma"]] <- list("Delta" = prior_sigma,
                                    "adj_sigma" = bal_sigma)
  return (prior_er)
}
