#' Function to read all, pure, and mixed variables from ER results.
#'
#' @param er_res a list output by the function \code{\link{plainER}} or \code{\link{priorER}}
#' @param y a vector of responses
#' @param alpha a numeric constant for the p-value threshold
#' @return a list including the clusters, pure variables, mixed variables, the significant
#' cluster indices, and the features found in the significant clusters
#' found by ER

readER <- function(er_res, y, alpha = 0.05) {
  #### get clusters with structure
  clusters <- recoverGroup(er_res$A)
  #### get all features included in clusters
  samp_feats <- unlist(clusters) %>% unique()
  #### get pure variables
  pure_vars <- pureRowInd(er_res$A)
  #### get mixed variables
  mix_vars <- setdiff(samp_feats, pure_vars)
  #### get significant clusters
  #p_vals <- 2 * pnorm(abs(er_res$beta) / sqrt(er_res$beta_var / length(y)), lower.tail = F)
  #Z_ind <- which(p_vals <= alpha)
  #A_hat <- er_res$A
  #sig_vars <- which(rowSums(abs(A_hat[,Z_ind, drop = F])) != 0)
  return(list("clusters" = clusters,
              "pure_vars" = pure_vars,
              "mix_vars" = mix_vars)) #,
              #"sig_clusters" = Z_ind,
              #"sig_vars" = sig_vars))
}
