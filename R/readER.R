#' Function to read all, pure, and mixed variables from ER results.
#'
#' @param er_res a list output by the function \code{\link{plainER}} or \code{\link{priorER}}
#' @return a list including the clusters, pure variables, and mixed varlues
#' found by ER

readER <- function(er_res) {
  clusters <- recoverGroup(er_res$A)
  samp_feats <- unlist(clusters) %>% unique()
  pure_vars <- pureRowInd(er_res$A)
  mix_vars <- setdiff(samp_feats, pure_vars)
  return(list("clusters" = clusters,
              "pure_vars" = pure_vars,
              "mix_vars" = mix_vars))
}
