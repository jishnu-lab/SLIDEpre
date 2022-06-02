#' Read Essential Regression Results.
#'
#' Function to read all, pure, and mixed variables from ER results.
#'
#' @importFrom magrittr '%>%'
#' @param er_res a list output by the function \link[EssReg]{plainER} or \link[EssReg]{priorER}
#' @param y a vector of responses
#' @return a list including the clusters, pure variables, mixed variables, the significant
#' cluster indices, and the features found in the significant clusters
#' found by ER
#' @export

readER <- function(er_res) {
  #### get clusters with structure
  clusters <- recoverGroup(er_res$A)
  #### get all features included in clusters
  samp_feats <- unlist(clusters) %>% unique()
  #### get pure variables
  pure_vars <- pureRowInd(er_res$A)
  #### get mixed variables
  mix_vars <- setdiff(samp_feats, pure_vars)
  return(list("clusters" = clusters,
              "pure_vars" = pure_vars,
              "mix_vars" = mix_vars))
}
