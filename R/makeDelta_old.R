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

makeDelta_old <- function(y, x, imp, er_res, type = "med", as_names = T) {
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
      ## add least frequent variable in that cluster to list of target vars
      least_freq <- as.numeric(clust_col[1, 1])
      targ_vars <- c(targ_vars, least_freq)
      ## get the corresponding row from freq_tab
      freq_row <- freq_tab[which(freq_tab$var == least_freq), ][3:ncol(freq_tab)]
      clusts_hit_row <- ifelse(freq_row == 1, least_freq, 0)
      clusts_hit <- rbind(clusts_hit, clusts_hit_row)
    }
  }
  ## get rid of NA variables - this cluster will be dropped
  targ_vars <- na.omit(targ_vars)
  ## save indices for later
  targ_ind <- targ_vars
  ## translate to feature names
  targ_vars <- indName(targ_vars, feat_names, F)

  ## find within cluster values
  clust_vals <- matrix(Inf, nrow = nrow(clusts_hit), ncol = ncol(clusts_hit))
  clust_vals <- data.frame(targ_vars, clust_vals)
  colnames(clust_vals) <- c("var", paste0("cluster", seq(1, nclust)))
  for (i in 1:ncol(clusts_hit)) { ## loop through the clusters
    cluster <- results$clusters[[i]]
    clust_targ <- clusts_hit[, i]
    clust_targ <- clust_targ[which(clust_targ > 0)]
    if (length(clust_targ) > 0) { ## if the cluster has mixed variables
      if (type == "med") { ## median
        clust_val <- findClustMed(samp_corr, cluster, clust_targ)
      } else if (type == "min") { ## min correlation
        clust_val <- findClustMin(samp_corr, cluster, clust_targ)
      } else { ## max correlation
        clust_val <- findClustMax(samp_corr, cluster, clust_targ)
      }
      for (j in 1:nrow(clust_vals)) { ## loop through the mixed variables of interest
        feature <- clust_vals$var[j]
        feat_ind <- indName(feature, feat_names, T) ## get feature index
        if (freq_tab[freq_tab$var == feat_ind, i + 2] == 1) { ## if that mixed var j appears in cluster i
          clust_vals[clust_vals$var == feature, i + 1] <- clust_val ## change the value for var j, clust i
        }
      }
    }
  }

  ## since we can't have corr(imp_var, mix_var) be two different values, we select the absolute
  ## minimum correlation found for mix_var across all clusters it appears in
  min_func <- function(x) {
    min_val <- min(abs(x))
    which_min <- which(abs(x) == min_val)
    which_min <- which_min[1]
    return (x[which_min])
  }
  clust_vals$min <- apply(clust_vals[, 2:ncol(clust_vals)], 1, min_func)

  ## subset data to just the mixed variables in clust_vals
  targ_data <- x[, targ_ind]
  change_corrs <- clust_vals$min
  change_data <- x
  for (k in 1:length(imp_excl)) { ## loop through excluded vars
    ## only change the target correlation if the value found above is more extreme
    ## than the value already in the sample correlation matrix
    targ_corrs <- ifelse(abs(samp_corr[imp_excl[k], targ_ind]) < abs(change_corrs), change_corrs, samp_corr[imp_excl[k], targ_ind])
    change_data[, imp_excl[k]] <- faux::rnorm_pre(targ_data, mu = 0, sd = 1, r = targ_corrs, empirical = T)
  }

  ## new sample correlation matrix
  new_samp_corr <- crossprod(change_data) / n

  ## make positive definite
  Delta <- makePosDef(new_samp_corr) %>% as.data.frame()
  colnames(Delta) <- feat_names
  rownames(Delta) <- feat_names

  return (as.matrix(Delta))
}
