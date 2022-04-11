#' Construct the expert knowledge matrix used as the target for the sample correlation matrix
#' reweighting. This is the first step of including prior knowledge in Essential Regression.
#' The column and row in the sample correlation matrix corresponding to the important feature
#' are adjusted in the following way:
#'
#' First, find the absolute maximum entry in the sample correlation matrix (excluding diagonal elements)
#' Then calculate the cutoff/threshold value used in determining node purity by \code{\link{FindRowMaxInd()}}
#' Next, calculate the values used in the replacement vector by subtracting the cutoff/threshold value from the
#' absolute maximum entry for each row. We then get a vector of values used to replace each of \eqn{\hat{\Sigma}_{imp, i}}
#' and \eqn{\hat{\Sigma}_{i, imp}}. If the \eqn{i}th replacement value in the replacement vector is greater in absolute
#' value than the value already found at \eqn{\hat{\Sigma}_{i, imp}} and \eqn{\hat{\Sigma}_{imp, i}}, then replace
#' that entry in \eqn{\hat{\Sigma}} with the replacement value. Otherwise, leave the entry alone.
#' Finally, we find the nearest positive definite matrix to our new \eqn{\hat{\Sigma}} with a diagonal of 1s.
#'
#' @param y a response vector of dimension \eqn{n}
#' @param x a data matrix of dimensions \eqn{n \times p}
#' @param imp the index of the important feature
#' @param er_res the output of a run of \code{\link{ER}}
#' @param equal_var a boolean indicating whether the columns of \code{x} are assumed to have equal variance
#' @return the expert knowledge matrix, \eqn{\Delta}, of dimensions \eqn{p \times p}
#' @export

makeDelta <- function(y, x, imp, er_res, change_all = F, equal_var = F) {
  ## must do an initial run of plain Essential Regression in order to get information about
  ## the data structure. use this for running ER with prior knowledge
  n <- nrow(x)
  p <- ncol(x)
  feat_names <- colnames(x)
  samp_corr <- crossprod(x) / n
  delta <- er_res$opt_delta

  ## get row maxes without last col/row
  abs_sc <- abs(samp_corr)
  diag(abs_sc) <- 0
  sc_ms <- findRowMax(abs_sc)
  max_inds <- sc_ms$max_inds
  max_vals <- sc_ms$max_vals

  if (equal_var) {
    se_est <- rep(1, p)
  } else {
    se_est <- apply(x, 2, sd) #### get sd of columns for feature matrix
  }

  ## read er_res
  er_results <- readER(er_res)
  pure_vars <- er_results$pure_vars
  mix_vars <- er_results$mix_vars

  new_imp <- c()
  for (i in 1:nrow(abs_sc)) {
    row_i <- abs_sc[i,]
    max_ind <- max_inds[i]
    cutoff <- (delta * se_est[i] * se_est[max_ind] + delta * se_est[i] * se_est)[imp]
    replacement <- sign(samp_corr[i, max_ind]) * (max_vals[i] - cutoff - 1e-10)
    if (change_all) {
      new_imp <- c(new_imp, replacement)
    } else {
      new_imp <- c(new_imp, ifelse(abs(replacement) > abs(samp_corr[i, imp]), replacement, samp_corr[i, imp]))
    }
  }

  if (change_all) {
    min_new_imp <- min(new_imp)
    new_imp <- rep(min_new_imp, length(new_imp))
  }

  Delta <- samp_corr
  Delta[imp, ] <- new_imp
  Delta[, imp] <- new_imp
  Delta[imp, imp] <- samp_corr[1, 1]

  if (!matrixcalc::is.positive.definite(Delta)) {
    ## use Matrix::nearPD so that diagonal is 1
    #Delta <- Matrix::nearPD(Delta, corr = T)$mat %>% as.matrix()
    Delta <- makePosDef(Delta) %>% as.matrix()
  }

  colnames(Delta) <- feat_names
  rownames(Delta) <- feat_names

  return (Delta)
}
