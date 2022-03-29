#' Perform prediction on the data using the ER predictor.
#'
#' @param Y a response vector of dimension \eqn{n}
#' @param X a data matrix of dimensions \eqn{n \times p}
#' @param sigma a correlation matrix of dimensions \eqn{p \times p}
#' @param A_hat a matrix of dimensions \eqn{p \times K}
#' @param Gamma_hat a matrix of dimensions
#' @param I_hat a matrix of dimensions
#' @return a list including \eqn{\hat{\theta}}, the predicted values, and \eqn{\Theta}

Prediction <- function(Y, X, sigma, A_hat, Gamma_hat, I_hat) {
  K <- ncol(A_hat)
  p <- ncol(X); n <- nrow(X)

  R <- matrix(0, nrow = K, ncol = p)
  Sigma <- sigma #crossprod(X) / n

  BI <- solve(t(A_hat[I_hat,]) %*% A_hat[I_hat,], t(A_hat[I_hat,]))
  R[, I_hat] = BI %*% (Sigma[I_hat, I_hat] - diag(Gamma_hat[I_hat]))
  R[, -I_hat] = BI %*% Sigma[I_hat, -I_hat]


  Q <- X %*% t(R)
  theta_hat <- try(t(R) %*% solve(crossprod(Q), crossprod(Q, Y)), silent = T)
  if (class(theta_hat)[1] == "try-error")
    theta_hat <- t(R) %*% MASS::ginv(crossprod(Q)) %*% crossprod(Q, Y)

  pred_val <- X %*% theta_hat
  return(list(theta = theta_hat, fit = pred_val, Theta = t(R)))
}
