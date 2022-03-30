#' Construct the expert knowledge matrix used as the target for the sample correlation matrix
#' reweighting. This is the first step of including prior knowledge in Essential Regression.
#'
#' @param y a response vector of dimension \eqn{n}
#' @param x a data matrix of dimensions \eqn{n \times p}
#' @param imp a vector of the features determined to be important by expert knowledge
#' @param er_res the output of a run of \code{\link{ER}}
#' @param type string indicating what type of replacement calculation to do
#' @param as_names boolean indicating whether the important features are reported as a vector of
#' indices (integers) or a vector of feature names (strings)
#' @return a matrix resulting from replacing the low correlations with the median or minimal value
#' in the given row/column. returns NULL if all important features in \code{imp} are already present
#' in the ER results
#' @export

makeDelta <- function(y, x, imp, er_res, type = "med", as_names = T) {
  ## must do an initial run of plain Essential Regression in order to get information about
  ## the data structure. use this for running ER with prior knowledge
  n <- nrow(x)
  p <- ncol(x)
  feat_names <- colnames(x)
  samp_corr <- crossprod(x) / n

  results <- readER(er_res)
  feats <- unlist(results$clusters) %>% unique()
  if (as_names) { ## translate to indices
    feats <- indName(feats, feat_names, F)
  }
  imp_excl <- setdiff(imp, feats)
  if (length(imp_excl) == 0) {
    cat("All important features are already included in ER results.")
    return (NULL)
  }

  ## get frequencies of mixed variables
  freq_tab <- getFreq(er_res)
  ## number of clusters
  nclust <- length(results$clusters)
  ## list of target vars
  targ_vars <- c()
  ## need to keep track of which clusters have been taken care of in targ_vars
  clusts_hit <- c()
  for (i in 1:nclust) {
    if (sum(clusts_hit[, i]) == 0 || is.na(sum(clusts_hit[, i]))) {
      ## get just the vars that appear in cluster i
      clust_col <- freq_tab[, c(1:2, i+2)] %>%
        dplyr::filter(.[[3]] == 1)
      ## add most frequent variable in that cluster to list of target vars
      most_freq <- as.numeric(clust_col[1, 1])
      targ_vars <- c(targ_vars, most_freq)
      ## get the corresponding row from freq_tab
      freq_row <- freq_tab[which(freq_tab$var == most_freq), ][3:ncol(freq_tab)]
      clusts_hit_row <- ifelse(freq_row == 1, most_freq, 0)
      clusts_hit <- rbind(clusts_hit, clusts_hit_row)
    }
  }
  ## get rid of NA variables - this cluster will be dropped
  targ_vars <- na.omit(targ_vars)
  ## translate to feature names
  targ_vars <- indName(targ_vars, feat_names, F)

  ## subset the sample correlation matrix to get just the submatrix of the
  ## variables that will be changed
  change_vars <- c(targ_vars, imp_excl)
  change_ind <- indName(change_vars, feat_names)
  sub_sc <- samp_corr[change_ind, change_ind]

  ## find maximum correlation of the target variables
  ## change the values of the entries of the important variables to these values
  for (i in 1:length(clusts_hit)) {
    cluster <- results$clusters[[i]]
    clust_targ <- clusts_hit[, i]
    clust_targ <- clust_targ[which(clust_targ > 0)]
    if (length(clust_targ) > 0) {
      if (type == "med") {
        clust_val <- findClustMed(samp_corr, cluster, clust_targ)
      } else if (type == "min") {
        clust_val <- findClustMin(samp_corr, cluster, clust_targ)
      } else {
        clust_val <- findClustMax(samp_corr, cluster, clust_targ)
      }
      for (j in 1:length(imp_excl)) {
        imp_var <- imp_excl[j]
        for (k in 1:length(clust_targ)) {
          targ_var <- indName(clust_targ[k], feat_names, F)
          sub_sc[imp_var, targ_var] <- ifelse(abs(clust_val) > abs(sub_sc[imp_var, targ_var]),
                                                 clust_val, sub_sc[imp_var, targ_var])
          sub_sc[targ_var, imp_var] <- ifelse(abs(clust_val) > abs(sub_sc[targ_var, imp_var]),
                                                 clust_val, sub_sc[targ_var, imp_var])
        }
      }
    }
  }

  ## make positive definite
  sub_sc_pd <- Matrix::nearPD(sub_sc, corr = T)
  sub_sc_pd <- sub_sc_pd$mat %>% as.matrix()

  ## make all other columns positive definite
  comp_sc_names <- setdiff(feat_names, rownames(sub_sc))
  comp_sub_sc <- samp_corr[comp_sc_names, comp_sc_names]
  comp_pd <- Matrix::nearPD(comp_sub_sc, corr = T)
  comp_pd <- comp_pd$mat %>% as.matrix()

  ## replace columns and rows of sample correlation matrix with the new
  ## submatrices constructed above
  Delta <- samp_corr
  Delta[rownames(sub_sc_pd), colnames(sub_sc_pd)] <- sub_sc_pd
  Delta[rownames(comp_pd), colnames(comp_pd)] <- comp_pd

  return (Delta)
}
