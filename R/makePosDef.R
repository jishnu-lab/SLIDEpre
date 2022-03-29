#' Transform the given matrix so that it is positive definite.
#'
#' @param samp_corr a sample correlation matrix of dimension \eqn{p \times p}
#' @return a positive definite matrix of the same dimension as \code{samp_corr}
#' @export

makePosDef <- function(samp_corr) {
  if (!matrixcalc::is.positive.definite(samp_corr)) {
    ## make positive definite
    sc_eig2 <- eigen(samp_corr)
    sc_eval2 <- ifelse(sc_eig2$values < 1e-10, 1e-4, sc_eig2$values)
    sc_adj2 <- sc_eig2$vectors %*% diag(sc_eval2) %*% t(sc_eig2$vectors)
    if (!matrixcalc::is.symmetric.matrix(sc_adj2)) {
      ## make symmetric
      sc_adj2_sym <- (sc_adj2 + t(sc_adj2)) / 2
      return (sc_adj2_sym)
    } else {
      return (sc_adj2)
    }
  } else {
    return (samp_corr)
  }
}

