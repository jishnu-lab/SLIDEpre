#' Hierarchical Essential Regression.
#'
#' Essential Regression with interaction terms among the components of \eqn{\mathbf{Z}}.
#'
#' @importFrom magrittr '%>%'
#' @param y a response vector of dimension \eqn{n}, must be continous
#' @param x a data matrix of dimensions \eqn{n \times p}
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
#' @param correction a boolean flag indicating whether to perform Bonferroni multiple testing correction
#' @param sig_betas either NULL or the proportion of \eqn{\beta}s to use when determining significant components
#' of \eqn{\mathbf{Z}}
#' @param niter the number of iterations to perform in the variable selection step
#' @return a list containing the results of the intial run of `plainER()`, the newly estimated \eqn{\beta},
#' the predicted values using the new \eqn{\beta}, and a list containing the results of variable
#' selection via Model X knockoffs
#' @export

hierER <- function(y, x, delta, thresh_fdr = 0.2, beta_est = "LS",
                   conf_int = T, pred = T, lambda = 0.1, rep_cv = 50, diagonal = F,
                   merge = F, equal_var = F, alpha_level = 0.05, support = NULL,
                   correction = T, sig_betas = NULL, niter = 20) {
  #### Data Housekeeping #######################################################
  raw_y <- y
  raw_x <- x

  #### run plainER() first #####################################################
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
                      correction = correction)

  x_std <- scale(x, T, T)
  y_std <- scale(y, T, T)
  zs <- predZ(x = x_std,
              er_res = plain_er)
  covar_mat <- crossprod(zs) / dim(zs)[2]
  covar_mat <- makePosDef(covar_mat)
  mus <- rep(0, dim(zs)[2])

  if (!is.null(sig_betas)) { ## top 5%
    sig_betas <- sigBetas(betas = plain_er$beta,
                          cutoff = sig_betas)
    marginal_vars <- c(unlist(sig_betas$pos_sig), unlist(sig_betas$neg_sig))
  } else { ## model x knockoffs
    marginal_vars <- gaussKOs(z = zs,
                              y = y_std,
                              statistic = knockoff::stat.glmnet_lambdasmax,
                              fdr = 0.1,
                              iters = niter)
  }

  ## Correct_it() ##############################################################
  if (is.null(dim(zs[, marginal_vars]))) { ## case with only one marg var
    n <- length(zs)
    covar_mat_marg <- crossprod(zs) / n
    ## just use ^-1 since one number
    y_hat <- (diag(1, nrow = length(y_std)) - zs %*% ((1 / covar_mat_marg) %*% t(zs))) %*% y_std
  } else {
    n <- nrow(zs)
    covar_mat_marg <- crossprod(zs) / n
    covar_mat_marg <- makePosDef(covar_mat_marg)
    y_hat <- (diag(1, nrow = length(y_std)) - zs %*% (solve(covar_mat_marg) %*% t(zs))) %*% y_std
  }

  ## make_interact_union() #####################################################
  colnames(zs) <- paste0('Z', seq(1, ncol(zs)))
  z_marg <- zs[, marginal_vars] %>% as.matrix()
  z_int <- diag(z_marg[, 1]) %*% as.matrix(zs)
  colnames(z_int) <- unlist(strsplit(paste0(colnames(z_marg), ".", colnames(zs)), " "))
  inter_colnames <- NULL
  for (i in 1:length(marginal_vars)) {
    inter_colnames <- c(inter_colnames, paste0(marginal_vars[i], ".", colnames(zs)))
  }
  for (i in 2:dim(z_marg)[2]) {
    z_int <- cbind(z_int, diag(z_int[, i]) %*% as.matrix(zs))
  }
  colnames(z_int) <- inter_colnames
  z_interact_list <- list(z_marginal = z_marg, z_interaction = z_int)

  ## Do Final Knockoffs ########################################################
  z_int <- z_interact_list$z_interaction
  sig_vars <- gaussKOs(z = z_int,
                       y = y_hat,
                       statistic = knockoff::stat.glmnet_lambdasmax,
                       fdr = 0.1,
                       iters = niter)
  sig_interact <- colnames(z_int)[sig_vars]

  hier_results <- list("marg_vars_sig" = marginal_vars,
                       "int_vars_sig" = sig_interact,
                       "int_vars" = z_int)

  marg_vars <- hier_results$marg_vars_sig ## marginal vars
  int_vars <- hier_results$int_vars_sig ## interaction vars
  z_int <- hier_results$int_vars ## interaction Zs

  #### Re-estimate Betas #######################################################
  ## first term is the significant marginal variables, selected from the original Zs
  ## second term is the significant interaction terms, selected from the interaction Zs
  new_betas <- coef(lm(y_std ~ zs[, marg_vars] + z_int[, int_vars]))
  pred_vals <- cbind(1, zs[, marg_vars], z_int[, int_vars]) %*% new_betas ## prediction


  return(list("plainER_results" = plain_er,
              "new_betas" = new_betas,
              "pred_vals" = pred_vals,
              "hier_results" = hier_results))
}
