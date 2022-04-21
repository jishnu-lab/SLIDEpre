#' Vanilla Essential Regression
#'
#' Perform Essential Regression without prior knowledge.
#'
#' @param y a response vector of dimension \eqn{n}
#' @param x a data matrix of dimensions \eqn{n \times p}
#' @param sigma a sample correlation matrix of dimensions \eqn{p \times p}
#' @param kept_sigma a matrix of dimensions \eqn{p \times p} with values 1/0 indicating
#' whether a given entry was removed via thresholding (1) or not (0) with function \link[EssReg]{threshSigma}
#' @param delta \eqn{\delta}, a numerical constant used for thresholding
#' @param thresh_fdr a numerical constant used for thresholding the correlation matrix to
#' control the false discovery rate, default is 0.2
#' @param beta_est a string indicating the type of estimation to use for \eqn{\beta}
#' @param conf_int a boolean indicating whether to calculate confidence intervals for the \eqn{\beta} estimates
#' @param pred a boolean indicating whether to do prediction
#' @param lambda \eqn{\lambda}, a numerical constant used in thresholding
#' @param rep_cv number of replicates for cross-validation
#' @param diagonal a boolean indicating the diagonal structure of the data ???
#' @param merge a boolean indicating the merge type
#' @param equal_var a boolean indicating whether there is equal variance ??
#' @param alpha_level \eqn{\alpha}, a numerical constant used in confidence interval calculation
#' @param support a boolean ???
#' @param correction a boolean flag indicating whether to perform Bonferroni multiple testing correction
#' @param verbose a boolean indicating whether to include printing
#' @param out_path a string path to where to save output
#' @return a list of results from the Essential Regression framework including: \eqn{K} = number of clusters,
#' \eqn{\hat{A}}, \eqn{\hat{C}}, \eqn{\hat{I}}, the indices of the pure variables, \eqn{\hat{\Gamma}},
#' \eqn{\hat{\beta}}, \eqn{\alpha}-level confidence intervals (if requested), prediction results (if requested),
#' the optimal value of \eqn{\lambda} determined by cross-validation, the optimal value of \eqn{\delta}
#' determined by cross-validation, \eqn{Q}, and the variances of \eqn{\hat{\beta}}
#' @export

plainER <- function(y, x, sigma, delta, thresh_fdr = 0.2, beta_est = "NULL",
                    conf_int = F, pred = T, lambda = 0.1, rep_cv = 50, diagonal = F,
                    merge = F, equal_var = F, alpha_level = 0.05, support = NULL,
                    correction = TRUE, verbose = F, out_path = NULL) {
  n <- nrow(x);  p <- ncol(x) #### feature matrix dimensions
  if (equal_var) {
    se_est <- rep(1, p)
  } else {
    se_est <- apply(x, 2, stats::sd) #### get sd of columns for feature matrix
  }
  #### scale delta by this value to satisfy some requirements so that the
  #### statistical guarantees in the paper hold
  delta_scaled <- delta * sqrt(log(max(p, n)) / n)

  #### save correlation matrix heatmap
  if (!is.null(out_path)) {
    pdf_file <- paste0(out_path, "/corr_mat_heatmap.pdf")
    dir.create(file.path(dirname(pdf_file)), showWarnings = F)
    grDevices::pdf(file = pdf_file)
    makeHeatmap(sigma, "Correlation Matrix Heatmap", T, T)
    grDevices::dev.off()
  }

  #### threshold sigma to control for FDR
  if (!is.null(thresh_fdr)) {
    control_fdr <- threshSigma(x = x,
                               sigma = sigma,
                               thresh = thresh_fdr)
    sigma <- control_fdr$thresh_sigma
    kept_entries <- control_fdr$kept_entries
  }

  #### save thresholding correlation matrix heatmap
  if (!is.null(out_path)) {
    pdf_file <- paste0(out_path, "/thresh_corr_mat_heatmap.pdf")
    dir.create(file.path(dirname(pdf_file)), showWarnings = F)
    grDevices::pdf(file = pdf_file)
    makeHeatmap(sigma, "FDR Thresholded Correlation Matrix Heatmap", T, T)
    grDevices::dev.off()
  }

  #### if delta has more than 1 element, then do rep_CV # of replicates
  #### of CV_Delta and select median of replicates
  opt_delta <- ifelse(length(delta_scaled) > 1,
                      stats::median(replicate(rep_cv, cvDelta(x = x, fdr_entries = kept_entries,
                                                              deltas_scaled = delta_scaled,
                                                              diagonal = diagonal, se_est = se_est,
                                                              merge = merge))),
                      delta_scaled)

  #### estimate membership matrix Ai
  #### also returns a vector of the indices of estimated pure variables
  #### and a list of the indices of estimated pure variables
  result_AI <- estAI(sigma = sigma, delta = opt_delta, se_est = se_est, merge = merge)
  pure_numb <- sapply(result_AI$pure_list, FUN = function(x) {length(c(x$pos, x$neg))})
  if (sum(pure_numb == 1) > 0) {
    cat("Changing ``merge'' to ``union'' and reselect delta ... \n")
    opt_delta <- ifelse(length(delta_scaled) > 1,
                        stats::median(replicate(rep_cv, cvDelta(x = x, fdr_entries = kept_entries,
                                                                deltas_scaled = delta_scaled,
                                                                diagonal = diagonal, se_est = se_est,
                                                                merge = F))),
                        delta_scaled)
    result_AI <- estAI(sigma = sigma, delta = opt_delta, se_est = se_est, merge = F)
  }

  A_hat <- result_AI$AI
  I_hat <- result_AI$pure_vec
  I_hat_list <- result_AI$pure_list

  if (is.null(I_hat)) {
    cat("Algorithm fails due to non-existence of pure variable.\n")
    stop()
  }

  C_hat <- estC(sigma = sigma, AI = A_hat, diagonal = diagonal)

  #### ER Supplement (2.1)
  #### Gamma_hat_{ii} = Sigma_hat_{ii} - A_hat_{i.}^T Sigma_hat_{Z} A_hat_{i.} for i in I_hat
  #### Gamma_hat_{ji} = 0 for all j ≠ i
  Gamma_hat <- rep(0, p)
  Gamma_hat[I_hat] <- diag(sigma[I_hat, I_hat]) - diag(A_hat[I_hat, ] %*% C_hat %*% t(A_hat[I_hat, ]))
  Gamma_hat[Gamma_hat < 0] <- 1e-2 #### replace negative values with something very close to 0

  #### Prediction --- what is Q??
  if (pred) {
    pred_result <- prediction(y = y, x = x, sigma = sigma, A_hat = A_hat,
                              Gamma_hat = Gamma_hat, I_hat = I_hat)

    #### theta_hat (supplement 2.2)
    theta_hat <- pred_result$theta_hat
    #### Inference In Latent Factor Regression With Clusterable Features
    #### Z_tilde = Q*X
    Q <- try(theta_hat %*% solve(crossprod(x %*% theta_hat) / n, crossprod(theta_hat)), silent = T)
    if (class(Q)[1] == "try-error") {
      Q <- theta_hat %*% MASS::ginv(crossprod(x %*% theta_hat) / n) %*% crossprod(theta_hat)
    }
  } else {
    pred_result <- Q <- NULL
  }

  #### Beta Estimation
  if (beta_est == "NULL" || beta_est == "Dantzig") {
    beta_hat <- beta_cis <- beta_var <- NULL
    if (beta_est == "Dantzig") {
      beta_hat <- estBetaDant(y = y, x = x, A_hat = A_hat, C_hat = C_hat,
                              I_hat = I_hat, delta = opt_delta, mu = 0.5,
                              lambda = 0.5)
    }
  } else {
    if (conf_int) {
      if (length(result_AI$pure_vec) != nrow(sigma)) {
        sigma_TJ <- estSigmaTJ(sigma = sigma, AI = A_hat, pure_vec = result_AI$pure_vec)
        # opt_lambda <- ifelse(length(lambda) > 1,
        #                  stats::median(replicate(rep_cv, cvLambda(x = x,
        #                                                           fdr_entries = kept_entries,
        #                                                           lambdas = lambda,
        #                                                           AI = result_AI$AI,
        #                                                           pure_vec = result_AI$pure_ec,
        #                                                           diagonal = diagonal))),
        #                  lambda)
        opt_lambda <- lambda
        if (opt_lambda > 0) {
          AI_hat <- abs(A_hat[I_hat, ]) ## just rows of pure variables
          sigma_bar_sup <- max(solve(crossprod(AI_hat), t(AI_hat)) %*% se_est[I_hat]) ## not sure what this does
          AJ <- estAJDant(C_hat = C_hat, sigma_TJ = sigma_TJ,
                          lambda = opt_lambda * opt_delta * sigma_bar_sup,
                          se_est_J = sigma_bar_sup + se_est[-I_hat])
        } else {
          AJ <- t(solve(C_hat, sigma_TJ))
        }
        A_hat[-result_AI$pure_vec, ] <- AJ
      }

      Gamma_hat[-I_hat] <- diag(sigma[-I_hat, -I_hat]) - diag(A_hat[-I_hat,] %*% C_hat %*% t(A_hat[-I_hat, ]))
      Gamma_hat[Gamma_hat < 0] <- 1e2 #### replace negative values with 100
    }
    res_beta <- estBeta(y = y, x = x, sigma = sigma, A_hat = A_hat,
                        C_hat = C_hat, Gamma_hat = Gamma_hat, I_hat = I_hat,
                        I_hat_list = I_hat_list, conf_int = conf_int,
                        alpha_level = alpha_level, correction = correction)
    beta_hat <- res_beta$beta_hat
    beta_conf_int <- res_beta$conf_int
    beta_var <- res_beta$beta_var
    return(list(K = ncol(A_hat),
                A = A_hat,
                C = C_hat,
                I = I_hat,
                I_clust = I_hat_list, ## original is I_ind
                Gamma = Gamma_hat,
                beta = beta_hat,
                beta_conf_int = beta_conf_int, ## original is beta_CIs
                beta_var = beta_var,
                pred = pred_result,
                opt_lambda = opt_lambda,
                opt_delta = opt_delta / sqrt(log(max(p, n)) / n),
                Q = Q,
                thresh_sigma = sigma))
  }
  return(list(K = ncol(A_hat),
              A = A_hat,
              C = C_hat,
              I = I_hat,
              I_clust = I_hat_list,
              Gamma = Gamma_hat,
              beta = beta_hat,
              pred = pred_result,
              opt_lambda = opt_lambda,
              opt_delta = opt_delta / sqrt(log(max(p, n)) / n),
              Q = Q,
              thresh_sigma = sigma))
}
