#' Find \eqn{\rho}, NOT USED ANYMORE
#'
#' Calculate the optimal value of \eqn{\rho} to use in \eqn{\Delta} that maximizes the marginal log-likelihood of the data, \eqn{x}, for fixed \eqn{\alpha}.
#'
#' @param x data matrix of dimensions \eqn{n \times p}
#' @param imp a vector of indices indicating important features
#' @param alpha \eqn{\alpha}, a positive value between 0 and 1
#' @param rho_range a range of values to try for \eqn{\rho} (should be some subset of \[-1, 1\])
#' @return a list containing the marginal log-likelihoods, the optimal value of \eqn{\rho},
#' the adjusted correlation matrix, and a plot of the marginal log-likelihoods
#' @importFrom foreach %dopar%
#' @export

findRho <- function(x, imp, alpha, rho_range) {
  i <- NULL
  temp <- foreach::foreach (i = 1:length(rho_range), .combine = rbind) %dopar% {
    rho <- rho_range[i]
    Delta <- transSigHat(x, imp, rho)
    Delta <- makePosDef(Delta)
    #colnames(Delta) <- seq(1, 21)
    #makeHeatmap(Delta, "", T, T)
    loglik <- logLik(alpha, Delta, x)
    c("rho" = rho, "loglik" = loglik)
  }
  ll_df <- as.data.frame(temp)
  rhos <- ll_df$rho
  logliks <- ll_df$loglik
  max_vals <- ll_df[which.max(ll_df$loglik),]
  out_plot <- ggplot2::ggplot(ll_df, ggplot2::aes(x = rho, y = loglik)) +
    ggplot2::geom_point() +
    ggplot2::geom_vline(xintercept = max_vals$rho, col = "red") +
    ggplot2::labs(title = "Marginal Log-Likelihoods", x = "rho", y = "marginal log-likelihood")

  Delta <- transSigHat(x, imp, max_vals$rho)
  Delta <- makePosDef(Delta)
  adj_mat <- adjMat(x, Delta, alpha)
  adapt_shrink_res <- list("marg_logliks" = ll_df,
                           "optim_rho" = max_vals$rho,
                           "adj_mat" = adj_mat,
                           "marg_plot" = out_plot)
  return (adapt_shrink_res)
}
