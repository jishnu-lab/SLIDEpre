#' Calculate the optimal value of \eqn{\alpha} that maximizes the marginal log-likelihood of the data, \eqn{x}, for a fixed \eqn{\Delta}.
#'
#' @param x data matrix of dimensions \eqn{n \times p}
#' @param Delta the expert knowledge matrix of dimensions \eqn{p \times p}, positive definite
#' @param alpha_range range of values of \eqn{\alpha} to conduct search for MLE over (should be some subset of \eqn{(0, 1)})
#' @return a list containing the marginal log-likelihoods, the optimal value of \eqn{\alpha}, the optimal value of \eqn{\rho} (if requested),
#' the adjusted correlation matrix, and a plot of the marginal log-likelihoods
#' @importFrom foreach %dopar%
#' @export

findAlpha <- function(x, Delta = NULL, alpha_range) {
    temp <- foreach::foreach (i = 1:length(alpha_range), .combine = rbind) %dopar% {
      alpha <- alpha_range[i]
      loglik <- logLik(alpha, Delta, x)
      c("alpha" = alpha, "loglik" = loglik)
    }

    ll_df <- as.data.frame(temp)
    alphas <- ll_df$alpha
    logliks <- ll_df$loglik
    max_vals <- ll_df[which.max(ll_df$loglik),]

    out_plot <- ggplot2::ggplot(ll_df, ggplot2::aes(x = alpha, y = loglik)) +
      ggplot2::geom_point() +
      ggplot2::geom_vline(xintercept = max_vals$alpha, col = "red") +
      ggplot2::labs(title = "Marginal Log-Likelihoods - Alpha", x = "alpha", y = "marginal log-likelihood")

    adj_mat <- adjMat(x, Delta, max_vals$alpha)

    adapt_shrink_res <- list("marg_logliks" = ll_df,
                             "optim_alpha" = max_vals$alpha,
                             "adj_mat" = adj_mat,
                             "marg_plot" = out_plot)

  return (adapt_shrink_res)
}
