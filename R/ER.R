#' Run Essential Regression.
#'
#' @param Y a response vector of dimension \eqn{n}
#' @param X a data matrix of dimensions \eqn{n \times p}
#' @param sigma a sample correlation matrix of dimensions \eqn{p \times p}
#' @param delta \eqn{\delta}, a numerical constant used for thresholding
#' @param beta_est a string indicating the type of estimation to use for \eqn{\beta}
#' @param CI a boolean indicating whether to calculate confidence intervals for the \eqn{\beta} estimates
#' @param pred a boolean indicating whether to do prediction
#' @param lbd \eqn{\lambda}, a numerical constant used in thresholding
#' @param rep_CV number of replicates for cross-validation
#' @param diagonal a boolean indicating the diagonal structure of the data ???
#' @param merge a boolean indicating the merge type
#' @param equal_var a boolean indicating whether there is equal variance ??
#' @param alpha_level \eqn{\alpha}, a numerical constant used in confidence interval calculation
#' @param support a boolean ???
#' @param correction a string indicating the type of multi-testing correction to perform
#' @param verbose a boolean indicating whether to include printing
#' @return a list of results from the Essential Regression framework including: \eqn{K = } number of clusters,
#' \eqn{\hat{A}}, \eqn{\hat{C}}, \eqn{\hat{I}}, the indices of the pure variables, \eqn{\hat{\Gamma}},
#' \eqn{\hat{\beta}}, \eqn{\alpha}-level confidence intervals (if requested), prediction results (if requested),
#' the optimal value of \eqn{\lambda} determined by cross-validation, the optimal value of \eqn{\delta}
#' determined by cross-validation, \eqn{Q}, and the variances of \eqn{\hat{\beta}}
#' @export

ER <- function(Y, X, sigma, delta, beta_est = "NULL", CI = F, pred = T, lbd = 1,
               rep_CV = 50, diagonal = F, merge = F, equal_var = F,
               alpha_level = 0.05, support = NULL, correction = "Bonferroni",
               verbose = F) {
  Sigma <- sigma
  n <- nrow(X);  p <- ncol(X) #### feature matrix dimensions
  if (equal_var) {
    se_est <- rep(1, p)
  } else {
    se_est <- apply(X, 2, sd) #### get sd of columns for feature matrix
  }
  #### scale delta by this value to satisfy some requirements so that the
  #### statistical guarantees in the paper hold
  deltaGrids <- delta * sqrt(log(max(p, n)) / n)

  #### if deltaGrids has more than 1 element, then do rep_CV # of replicates
  #### of CV_Delta and select median of replicates
  optDelta <- ifelse(length(deltaGrids) > 1,
                     median(replicate(rep_CV, CV_Delta(X, deltaGrids, diagonal, se_est, merge))),
                     deltaGrids)
  if (verbose) { #### select first delta in deltaGrids that is >= optDelta
    cat("Selecting the delta =", delta[min(which(deltaGrids >= optDelta))], "\n")
  }

  #### estimate membership matrix Ai
  #### also returns a vector of the indices of estimated pure variables
  #### and a list of the indices of estimated pure variables
  resultAI <- EstAI(Sigma, optDelta, se_est, merge)
  pure_numb <- sapply(resultAI$pureSignInd, FUN = function(x) {length(c(x$pos, x$neg))})
  if (sum(pure_numb == 1) > 0) {
    cat("Changing ``merge'' to ``union'' and reselect delta ... \n")
    optDelta <- ifelse(length(deltaGrids) > 1,
                       median(replicate(rep_CV, CV_Delta(X, deltaGrids, diagonal, se_est, merge = F))),
                       deltaGrids)
    resultAI <- EstAI(Sigma, optDelta, se_est, merge = F)
  }

  A_hat <- resultAI$AI
  I_hat <- resultAI$pureVec
  I_hat_ind <- resultAI$pureSignInd

  if (is.null(I_hat)) {
    cat("Algorithm fails due to non-existence of pure variable.\n")
    stop()
  }

  C_hat <- EstC(Sigma, A_hat, diagonal)
  Gamma_hat <- rep(0, p)
  Gamma_hat[I_hat] <- diag(Sigma[I_hat, I_hat]) - diag(A_hat[I_hat,] %*% C_hat %*% t(A_hat[I_hat,]))
  Gamma_hat[Gamma_hat < 0] <- 1e-2

  if (pred) {
    pred_result <- Prediction(Y, X, Sigma, A_hat, Gamma_hat, I_hat)

    # the matrix to predict Z
    Theta_hat <- pred_result$Theta
    Q <- try(Theta_hat %*% solve(crossprod(X %*% Theta_hat) / n, crossprod(Theta_hat)), silent = T)
    if (class(Q)[1] == "try-error") {
      Q <- Theta_hat %*% ginv(crossprod(X %*% Theta_hat) / n) %*% crossprod(Theta_hat)
    }
  } else {
    pred_result <- Q <- NULL
  }



  if (beta_est == "NULL" || beta_est == "Dantzig") {
    beta_hat <- beta_CIs <- beta_var <- NULL
    if (beta_est == "Dantzig") {
      beta_hat <- Est_beta_dz(Y, X, A_hat, C_hat, I_hat, optDelta, 0.5, 0.5)
    }
  } else {
    if (CI) {
      if (length(resultAI$pureVec) != nrow(Sigma)) {
        Y_hat <- EstY(Sigma, A_hat, resultAI$pureVec)
        optLbd <- ifelse(length(lbd) > 1,
                         median(replicate(rep_CV, CV_lbd(X, lbd, resultAI$AI, resultAI$pureVec, diagonal))),
                         lbd)
        lbd <- optLbd
        if (lbd > 0) {
          AI_hat <- abs(A_hat[I_hat, ]) ## just rows of pure variables
          sigma_bar_sup <- max(solve(crossprod(AI_hat), t(AI_hat)) %*% se_est[I_hat])
          AJ <- EstAJDant(C_hat, Y_hat, lbd * optDelta * sigma_bar_sup, sigma_bar_sup + se_est[-I_hat])
        } else {
          AJ <- t(solve(C_hat, Y_hat))
        }
        A_hat[-resultAI$pureVec, ] <- AJ
      }

      Gamma_hat[-I_hat] <- diag(Sigma[-I_hat, -I_hat]) - diag(A_hat[-I_hat,] %*% C_hat %*% t(A_hat[-I_hat,]))
      Gamma_hat[Gamma_hat < 0] <- 1e2


    }
    res_beta <- Est_beta(Y, X, Sigma, A_hat, C_hat, Gamma_hat, I_hat, I_hat_ind, CI = CI,
                         alpha_level = alpha_level, correction = correction)
    beta_hat <- res_beta$beta
    beta_CIs <- res_beta$CIs
    beta_var <- res_beta$beta_var
  }
  return(list(K = ncol(A_hat), A = A_hat, C = C_hat, I = I_hat, I_ind = I_hat_ind, Gamma = Gamma_hat,
              beta = beta_hat, beta_CIs = beta_CIs, pred = pred_result, lbd = lbd,
              optDelta = delta[min(which(deltaGrids >= optDelta))],
              Q = Q, beta_var = beta_var))
}
