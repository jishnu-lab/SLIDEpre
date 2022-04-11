#' Cross validation for choosing \eqn{\delta}. For each delta from the given grids,
#' first split the data into two data sets. Obtain \eqn{I}, \eqn{A_I} and \eqn{C} from data set 1.
#' Then calculate \eqn{A_I \cdot C \cdot A_I^\top} and choose \eqn{\delta} which minimizes the criterion
#' \eqn{A_I \cdot C \cdot A_I^\top - \Sigma(\text{data set 2})}
#'
#' @param x data matrix of dimensions \eqn{n \times p}
#' @param delta_scaled a vector of numerical constants over which to perform the search for the optimal \eqn{\delta}
#' @param diagonal a boolean indicating the diagonal structure of \eqn{C}
#' @param se_est estimated standard errors
#' @param merge a boolean indicating merge style
#' @return the selected optimal \eqn{\delta}

cvDelta <- function(x, deltas_scaled, diagonal, se_est, merge) {
  n <- nrow(x); p <- ncol(x)
  samp_ind <- sample(n, floor(n / 2))
  x_train <- x[samp_ind, ]
  x_val <- x[-samp_ind, ]
  sigma_train <- crossprod(x_train) / nrow(x_train);
  diag(sigma_train) <- 0
  sigma_val <- crossprod(x_val) / nrow(x_val)

  result_max <- findRowMax(abs(sigma_train))
  max_vals <- result_max$max_vals
  max_inds <- result_max$max_inds

  loss <- c()
  for (i in 1:length(deltas_scaled)) {
    result_fitted <- calFittedSigma(sigma = sigma_train, delta = deltas_scaled[i],
                                   max_vals = max_vals, max_inds = max_inds,
                                   se_est = se_est, diagonal = diagonal, merge = merge)
    fit_sigma <- result_fitted$fit_sigma
    pure_vec <- result_fitted$pure_vec

    if (is.null(dim(fit_sigma)) && fit_sigma == -1) {
      loss[i] <- Inf
    } else {
      denom <- length(pure_vec) * (length(pure_vec) - 1)
      loss[i] <- 2 * offSum(sigma_val[pure_vec, pure_vec], fit_sigma, se_est[pure_vec]) / denom
    }
  }
  return(deltas_scaled[which.min(loss)])
}
