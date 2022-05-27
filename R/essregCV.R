#' Cross-Validation For Essential Regression
#'
#' Perform k-fold cross-validation for essential regression to select
#'
#' @importFrom magrittr '%>%'
#' @importFrom foreach '%dopar%'
#' @param k an integer for number of folds to use in cross-validation
#' @param y a response vector of dimension \eqn{n}
#' @param x a data matrix of dimensions \eqn{n \times p}
#' @param priors a vector of names or indices of important features in data matrix if
#' using \link{priorER} or \code{NULL} if using \link{plainER}
#' @param delta \eqn{\delta}, a numerical constant used for thresholding
#' @param thresh_fdr a numerical constant used for thresholding the correlation matrix to
#' control the false discovery rate, default is 0.2
#' @param perm_option a string indicating the type of permutation type do perform
#' (can be NULL, "x", "x_y", or "y_before_split")
#' @param y_factor a boolean flag indicating whether \eqn{y} is categorical (\code{T}) or not (\code{F})
#' @param beta_est a string indicating the type of estimation to use for \eqn{\beta}
#' @param sel_corr a boolean flag indicating whether to perform cross-validation by evaluating the correlation
#' between the predicted and true values of \eqn{y} (\code{T}) or by evaluating the prediction error via mse or auc (\code{F})
#' @param lambda \eqn{\lambda}, a numerical constant used in thresholding
#' @param rep_cv number of replicates for cross-validation
#' @param diagonal a boolean indicating the diagonal structure of the data ???
#' @param merge a boolean indicating the merge type
#' @param equal_var a boolean indicating whether there is equal variance ??
#' @param alpha_level \eqn{\alpha}, a numerical constant used in confidence interval calculation
#' @param thresh a numerical constant used as the threshold for convergence in Variational Bayes
#' @param support a boolean ???
#' @param svm a boolean flag indicating whether to fit svm/svr and use that for prediction
#' @param correction a boolean flag indicating whether to perform Bonferroni multiple testing correction
#' @param change_all a boolean indicating whether to change all entries in \eqn{\hat{\Sigma}}
#' for an important feature (T) or to just change to more extreme values (F)
#' @param rep an integer indicating the replicate number
#' @param out_path path for saving output
#' @return An object of class \sQuote{data.frame}
#' @export

essregCV <- function(k = 5, y, x, priors = NULL, delta, thresh_fdr = 0.2, lambda = 0.1,
                     rep_cv = 50, alpha_level = 0.05, thresh = 0.001, perm_option = NULL,
                     beta_est = "NULL", sel_corr = T, svm = F, y_factor = F, out_path,
                     diagonal = F, merge = F, equal_var = F, support = NULL, correction = T,
                     change_all = F, delta_grid = NULL, rep) {

  if (y_factor) {
    eval_type <- "auc"
  } else {
    eval_type <- "mse"
  }

  ## create output directory
  new_dir <- paste0(out_path, "replicate", rep, "/")
  dir.create(file.path(new_dir), showWarnings = F, recursive = T)

  ## divide into folds
  ## first part is partition()
  total_num <- nrow(x)
  num_group <- k
  remainder <- total_num %% num_group # get the remainder
  num_per_group <- total_num %/% num_group
  partition <- rep(num_per_group, num_group) + c(rep(1, remainder), rep(0, num_group - remainder))
  ## second part is extract()
  pre_vec <- sample(1:nrow(x))
  extract <- vector("list", length(partition))
  extract[[1]] <- pre_vec[1:partition[1]]
  for (i in 2:length(partition)) {
    temp_ind <- sum(partition[1:(i - 1)]) + 1
    extract[[i]] <- pre_vec[temp_ind:(temp_ind + partition[i] - 1)]
  }
  group_inds <- extract

  ## methods list
  methods <- c("plainER", "plainER_allZs", "plainER_IVS", "lasso", "IVS")
  if (!is.null(priors)) {
    methods <- c(methods, "priorER", "priorER_allZs", "priorER_IVS")
  }

  if (!is.null(perm_option)) {
    if (perm_option == "x") { ## permute columns of x
      perm_col_ind <- sample(1:ncol(x))
      methods <- c(methods, paste0(methods, "_x"))
    } else if (perm_option == "y_before_split") { ## permute y before splitting into train/valid
      perm_row_ind <- sample(1:nrow(x))
      methods <- c(methods, paste0(methods, "_ybs"))
    } else if (perm_option == "x_y") {
      perm_col_ind <- sample(1:ncol(x))
      perm_row_ind <- sample(1:nrow(x))
      methods <- c(methods, paste0(methods, "_xy"))
    } else {
      methods <- c(methods, paste0(methods, "_y"))
    }
  }

  ## initialization
  results <- NULL
  if (sel_corr) {
    for (h in 1:length(methods)) {
      method <- methods[h]
      results[[method]] <- list()
    }
  }

  for (i in 1:k) { ## loop through folds
    cat("FOLD ", i, ". . . . \n")
    valid_ind <- group_inds[[i]] ## validation indices
    train_y <- y[-valid_ind] ## training y's
    valid_y <- y[valid_ind] ## validation y's
    train_x <- x[-valid_ind, ] ## training x's
    valid_x <- matrix(x[valid_ind, ], ncol = ncol(x)) ## validation x's

    for (j in 1:length(methods)) { ## loop through methods
      method_j <- methods[j]
      cat("CV for ", method_j, ". . . \n")

      ## permute and standardize sets
      if (grepl(x = method_j, pattern = "ybs", fixed = TRUE)) {
        train_y <- y[perm_row_ind][-valid_ind]
        stands <- standCV(train_y = train_y,
                          train_x = train_x,
                          valid_y = valid_y,
                          valid_x = valid_x)
      } else if (grepl(x = method_j, pattern = "x_y", fixed = TRUE)) {
        train_x <- train_x[, perm_col_ind]
        train_y <- y[perm_row_ind][-valid_ind]
        stands <- standCV(train_y = train_y,
                          train_x = train_x,
                          valid_y = valid_y,
                          valid_x = valid_x)
      } else if (grepl(x = method_j, pattern = "x", fixed = TRUE)) {
        train_x <- train_x[, perm_col_ind]
        stands <- standCV(train_y = train_y,
                          train_x = train_x,
                          valid_y = valid_y,
                          valid_x = valid_x)
      } else if (grepl(x = method_j, pattern = "y", fixed = TRUE)) { ## just permute y
        perm_ind <- sample(1:nrow(train_x))
        train_y <- train_y[perm_ind]
        stands <- standCV(train_y = train_y,
                          train_x = train_x,
                          valid_y = valid_y,
                          valid_x = valid_x)
      } else {
        stands <- standCV(train_y = train_y,
                          train_x = train_x,
                          valid_y = valid_y,
                          valid_x = valid_x)
      }

      train_x_std <- stands$train_x
      train_y_std <- stands$train_y
      valid_x_std <- stands$valid_x
      valid_y_std <- stands$valid_y

      if (grepl(x = method_j, pattern = "plainER_IVS", fixed = TRUE)) { ## plain essential regression with IVS
        load(paste0(new_dir, methods[1], "_fold", i, ".rda")) ## load in the plainER results from other loop

        ## calculate predicted values
        train_z <- predZ(x = train_x_std, er_res = res)
        ivs <- IVS(y = train_y_std, z = train_z)
        valid_z <- predZ(x = valid_x_std, er_res = res)

        new_betas <- betaBMA(x = train_x,
                             y = train_y,
                             er_res = res,
                             imps = NULL,
                             imps_z = ivs,
                             estim = "HPM")
        beta_train <- cbind(train_z, 1) %*% new_betas
        #beta_valid <- valid_z[, ivs] %*% new_betas
        beta_valid <- cbind(valid_z, 1) %*% new_betas
        pred_vals <- beta_valid
      } else if (grepl(x = method_j, pattern = "priorER_IVS", fixed = TRUE)) { ## prior essential regression with IVS
        load(paste0(new_dir, methods[2], "_fold", i, ".rda")) ## load in the priorER results from other loop

        ## calculate predicted values
        train_z <- predZ(x = train_x_std, er_res = res)
        ivs <- IVS(y = train_y_std, z = train_z)
        valid_z <- predZ(x = valid_x_std, er_res = res)


        new_betas <- betaBMA(x = train_x,
                             y = train_y,
                             er_res = res,
                             imps = priors,
                             imps_z = ivs,
                             estim = "HPM")
        beta_train <- cbind(train_z, 1) %*% new_betas
        #beta_valid <- valid_z[, ivs] %*% new_betas
        beta_valid <- cbind(valid_z, 1) %*% new_betas
        pred_vals <- beta_valid
      } else if (grepl(x = method_j, pattern = "plainER_allZs", fixed = TRUE)) { ## plain essential regression, predict with all Zs
        load(paste0(new_dir, methods[1], "_fold", i, ".rda")) ## load in the plainER results from other loop

        ## get things for svm
        pred_all_betas <- res$pred$er_predictor
        beta_train <- train_x_std %*% pred_all_betas
        beta_valid <- valid_x_std %*% pred_all_betas
        pred_vals <- beta_valid
      } else if (grepl(x = method_j, pattern = "priorER_allZs", fixed = TRUE)) { ## prior essential regression with all Zs
        load(paste0(new_dir, methods[2], "_fold", i, ".rda")) ## load in the priorER results from other loop

        ## get things for svm
        pred_all_betas <- res$priorER_results$pred$er_predictor
        beta_train <- train_x_std %*% pred_all_betas
        beta_valid <- valid_x_std %*% pred_all_betas
        pred_vals <- beta_valid
      } else if (grepl(x = method_j, pattern = "plainER", fixed = TRUE)) { ## plain essential regression
        res <- plainER(y = train_y,
                       x = train_x,
                       sigma = NULL,
                       delta = delta,
                       lambda = lambda,
                       thresh_fdr = thresh_fdr,
                       rep_cv = rep_cv,
                       alpha_level = alpha_level,
                       beta_est = beta_est,
                       conf_int = T,
                       pred = T,
                       diagonal = diagonal,
                       merge = merge,
                       equal_var = equal_var,
                       support = support,
                       correction = correction)

        ## calculate predicted values
        train_z <- predZ(x = train_x_std, er_res = res)
        sig_betas <- sigBetas(betas = res$beta, cutoff = 0.1) ## Top 5% positive/negative
        sig_betas <- c(unlist(sig_betas$pos_sig), unlist(sig_betas$neg_sig))
        valid_z <- predZ(x = valid_x_std, er_res = res)

        new_betas <- betaBMA(x = train_x,
                             y = train_y,
                             er_res = res,
                             imps = NULL,
                             imps_z = sig_betas,
                             estim = "HPM")
        beta_train <- cbind(train_z, 1) %*% new_betas
        #beta_valid <- valid_z[, ivs] %*% new_betas
        beta_valid <- cbind(valid_z, 1) %*% new_betas
        pred_vals <- beta_valid
      } else if (grepl(x = method_j, pattern = "priorER", fixed = TRUE)) { ## prior essential regression
        res <- priorER(y = train_y,
                       x = train_x,
                       imps = priors,
                       sigma = NULL,
                       delta = delta,
                       lambda = lambda,
                       thresh_fdr = thresh_fdr,
                       thresh = thresh,
                       rep_cv = rep_cv,
                       alpha_level = alpha_level,
                       beta_est = beta_est,
                       conf_int = T,
                       pred = T,
                       diagonal = diagonal,
                       merge = merge,
                       equal_var = equal_var,
                       support = support,
                       correction = correction,
                       change_all = change_all)

        ## calculate predicted values
        train_z <- predZ(x = train_x_std, er_res = res$priorER_results)
        sig_betas <- sigBetas(betas = res$priorER_results$beta, cutoff = 0.1) ## Top 5% positive/negative
        sig_betas <- c(unlist(sig_betas$pos_sig), unlist(sig_betas$neg_sig))
        valid_z <- predZ(x = valid_x_std, er_res = res$priorER_results)

        new_betas <- betaBMA(x = train_x,
                             y = train_y,
                             er_res = res,
                             imps = priors,
                             imps_z = sig_betas,
                             estim = "HPM")
        beta_train <- cbind(train_z, 1) %*% new_betas
        #beta_valid <- valid_z[, ivs] %*% new_betas
        beta_valid <- cbind(valid_z, 1) %*% new_betas
        pred_vals <- beta_valid
      } else if (grepl(x = method_j, pattern = "IVS", fixed = TRUE)) { ## just IVS
        ivs <- IVS(train_y_std, train_x_std)
        ivs_x <- train_x_std[, ivs]
        new_betas <- coef(lm(train_y_std ~ ivs_x))

        beta_train <- cbind(1, ivs_x) %*% new_betas
        beta_valid <- cbind(1, valid_x_std[, ivs]) %*% new_betas
        pred_vals <- beta_valid
      } else { ## lasso for comparison
        if ((nrow(train_x_std) / 10) < 3) { ## sample size too small
          res <- glmnet::cv.glmnet(train_x_std, train_y_std, alpha = 1, nfolds = 5, standardize = F, grouped = F)
        } else {
          res <- glmnet::cv.glmnet(train_x_std, train_y_std, alpha = 1, nfolds = 10, standardize = F, grouped = F)
        }
        beta_hat <- coef(res, s = res$lambda.min)[-1]
        sub_beta_hat <- which(beta_hat != 0)
        if (length(sub_beta_hat) == 0) { ## if lasso selects no variable, randomly pick 5 features instead
          cat("Lasso selects no features - Randomly selecting 5 features. . . \n")
          sub_beta_hat <- sample(1:ncol(train_x_std), 5)
        }
        ## get things for svm
        beta_train <- train_x_std[, sub_beta_hat]
        beta_valid <- valid_x_std[, sub_beta_hat, drop = F]
        pred_vals <- glmnet::predict.glmnet(res$glmnet.fit, valid_x_std, s = res$lambda.min)
      }

      save(res, train_x, train_y, valid_x, valid_y, stands, valid_ind, file = paste0(new_dir, method_j, "_fold", i, ".rda"))

      if (svm) { ## if svm flag == TRUE, use svm/svr to get predicted values for validation set
        ## fit svm
        res_tune <- e1071::tune.svm(x = beta_train, y = train_y_std, cost = 2 ^ seq(-3, 3, 1)) ## linear kernel
        ## get predicted values on validation set
        pred_vals <- stats::predict(object = res_tune$best.model, newdata = beta_valid)
      }

      if (sel_corr) { ## if using correlation to evaluate model fit
        method_res <- results[[method_j]]
        method_res[[length(method_res) + 1]] <- cbind(valid_y_std, pred_vals)
        results[[method_j]] <- method_res
      } else {
        if (eval_type == "auc") {
          pred_obj <- prediction(pred_vals, valid_y_std)
          perf_rates <- performance(pred_obj, "tpr", "fpr")
          perf_auc <- performance(pred_vals, "auc")
          auc <- as.numeric(perf_auc@y.values)
          iter_res <- c(i, method_j, auc, perf@y.values[[1]], perf@xvalues[[1]])
          results <- rbind(results, iter_res)
        } else {
          mse <- mean((valid_y_std - pred_vals)^2)
          iter_res <- c(i, method_j, mse)
          results <- rbind(results, iter_res)
        }
      }
    }
  }
  ## set results data frame column names
  if (sel_corr) {
    res_df <- NULL
    for (l in 1:length(methods)) {
      method <- methods[l]
      method_res <- results[[method]]
      unlist_method_res <- as.data.frame(do.call(rbind, method_res))
      spear_corr <- cor(unlist_method_res[, 1], unlist_method_res[, 2], method = "spearman") ## calc spearman corr
      res_df <- rbind(res_df, c(method, spear_corr))
    }
    colnames(res_df) <- c("method", "spearman_corr")
    res_df <- res_df %>%
      as.data.frame() %>%
      dplyr::mutate(spearman_corr = as.numeric(as.character(spearman_corr)))
    saveRDS(res_df, file = paste0(new_dir, "output_table.rds"))
    return (res_df)
  } else {
    if (eval_type == "mse") {
      colnames(results) <- c("fold", "method", "mse")
      results <- as.data.frame(results)
      results <- results %>%
        dplyr::group_by(method) %>%
        dplyr::summarise(mean_mse = mean(as.numeric(mse)))
    } else {
      colnames(results) <- c("fold", "method", "auc", "tpr", "fpr")
      results <- as.data.frame(results)
      results <- results %>%
        dplyr::group_by(method) %>%
        dplyr::summarise(mean_auc = mean(as.numeric(auc)),
                         mean_tpr = mean(as.numeric(tpr)),
                         mean_fpr = mean(as.numeric(fpr)))
    }
    saveRDS(results, file = paste0(new_dir, "output_table.rds"))
    return (results)
  }
}
