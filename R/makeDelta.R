#' Construct \eqn{\Delta}
#'
#' Construct the expert knowledge matrix used as the target for the sample correlation matrix
#' reweighting. This is the first step of including prior knowledge in Essential Regression.
#' The column and row in the sample correlation matrix corresponding to the important feature
#' are adjusted in the following way:
#'
#' First, find the absolute maximum entry in the sample correlation matrix (excluding diagonal elements)
#' Then calculate the cutoff/threshold value used in determining node purity by \link[EssReg]{findRowMaxInd}
#' Next, calculate the values used in the replacement vector by subtracting the cutoff/threshold value from the
#' absolute maximum entry for each row. We then get a vector of values used to replace each of \eqn{\hat{\Sigma}_{imp, i}}
#' and \eqn{\hat{\Sigma}_{i, imp}}. If the \eqn{i}th replacement value in the replacement vector is greater in absolute
#' value than the value already found at \eqn{\hat{\Sigma}_{i, imp}} and \eqn{\hat{\Sigma}_{imp, i}}, then replace
#' that entry in \eqn{\hat{\Sigma}} with the replacement value. Otherwise, leave the entry alone.
#' Finally, we find the nearest positive definite matrix to our new \eqn{\hat{\Sigma}} with a diagonal of 1s.
#'
#' @param x a data matrix of dimensions \eqn{n \times p}
#' @param sigma a sample correlation matrix of dimensions \eqn{p \times p}
#' @param imps a vector of indices of the important feature
#' @param er_res the output of a run of \link[EssReg]{plainER}
#' @param change_all a boolean flag indicating whether to adjust all correlations or just those that are smaller
#' in absolute value than the value to change to
#' @param equal_var a boolean indicating whether the columns of \code{x} are assumed to have equal variance
#' @return the expert knowledge matrix, \eqn{\Delta}, of dimensions \eqn{p \times p}
#' @export

makeDelta <- function(x, sigma, imps, er_res, change_all = F, equal_var = F) {
  ## must do an initial run of plain Essential Regression in order to get information about
  ## the data structure. use this for running ER with prior knowledge
  n <- nrow(x)
  p <- ncol(x)
  feat_names <- colnames(x)
  delta <- er_res$opt_delta * sqrt(log(max(p, n)) / n)

  ## get row maxes without last col/row
  abs_sc <- abs(sigma)
  diag(abs_sc) <- 0
  sc_ms <- findRowMax(abs_sc)
  max_inds <- sc_ms$max_inds
  max_vals <- sc_ms$max_vals

  if (equal_var) {
    se_est <- rep(1, p)
  } else {
    se_est <- apply(x, 2, stats::sd) #### get sd of columns for feature matrix
  }

  ## read er_res
  er_results <- readER(er_res)
  pure_vars <- er_results$pure_vars
  mix_vars <- er_results$mix_vars

  new_imps <- matrix(0, nrow = length(imps), ncol = p)
  for (i in 1:nrow(abs_sc)) {
    row_i <- abs_sc[i,]
    max_ind <- max_inds[i]
    cutoffs <- rep(2 * delta, length(imps)) #(delta * se_est[i] * se_est[max_ind] + delta * se_est[i] * se_est)[imps]
    replacements <- max_vals[i] - cutoffs - 1e-10
    ## adjust sign to match original correlation, or set to 0 if cutoffs > max_vals[i]
    replacements <- ifelse(replacements < 0, 0, sign(sigma[i, max_ind]) * replacements)
    if (change_all) { ## if changing entire row/column
      new_imps[, i] <- replacements
    } else { ## if only changing row/column entries that are smaller than replacements
      replacements <- ifelse(abs(replacements) > abs(sigma[i, imps]), replacements, sigma[i, imps])
      new_imps[, i] <- replacements
    }
  }

  if (change_all) { ## if changing entire row/column
    min_new_imps <- apply(new_imps, 1, min)
    new_imps <- matrix(min_new_imps, nrow = length(min_new_imps), ncol = p, byrow = FALSE)
  }

  Delta <- sigma
  for (i in 1:nrow(new_imps)) { ## replace values
    replace_row <- new_imps[i, ]
    imp <- imps[i]
    Delta[imp, ] <- replace_row
    Delta[, imp] <- replace_row
    Delta[imp, imp] <- sigma[1, 1] ## replace diagonal elements
  }

  if (!matrixcalc::is.positive.definite(Delta)) { ## make positive definite
    Delta <- makePosDef(Delta)
    Delta <- as.matrix(Delta)
  }

  colnames(Delta) <- feat_names
  rownames(Delta) <- feat_names

  return (Delta)
}
