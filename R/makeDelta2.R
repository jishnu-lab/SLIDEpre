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

makeDelta2 <- function(y, x, imp, er_res, equal_var = F) {
  ## must do an initial run of plain Essential Regression in order to get information about
  ## the data structure. use this for running ER with prior knowledge
  n <- nrow(x)
  p <- ncol(x)
  feat_names <- colnames(x)
  samp_corr <- crossprod(x) / n
  delta <- er_res$optDelta

  ## get row maxes without last col/row
  abs_sc <- abs(samp_corr)
  diag(abs_sc) <- 0
  sc_ms <- FindRowMax(abs_sc)
  arg_Ms <- sc_ms$arg_M
  Ms <- sc_ms$M

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
    arg_M <- arg_Ms[i]
    cutoff <- (delta * se_est[i] * se_est[arg_M] + delta * se_est[i] * se_est)[imp]
    replacement <- sign(samp_corr[i, arg_M]) * (Ms[i] - cutoff)
    new_imp <- c(new_imp, ifelse(abs(replacement) > abs(samp_corr[i, imp]), replacement, samp_corr[i, imp]))
  }

  Delta <- samp_corr
  Delta[imp, ] <- new_imp
  Delta[, imp] <- new_imp

  if (!matrixcalc::is.positive.definite(Delta)) {
    ## use Matrix::nearPD so that diagonal is 1
    Delta <- Matrix::nearPD(Delta, corr = T)$mat %>% as.matrix()
  }

  return (Delta)
}
