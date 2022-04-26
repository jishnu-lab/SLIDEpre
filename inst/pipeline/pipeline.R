## load libraries
library(EssReg)
library(dplyr)
library(doParallel)
library(foreach)

## FILL IN #####
yaml_path = "/Users/aerosengart/Documents/Das Lab/test.yaml"
priors = NULL
cores = 6
################

## set up parallelization
doParallel::registerDoParallel(cores = cores)

## process arguments
er_input <- yaml::yaml.load_file(yaml_path)
x <- as.matrix(utils::read.csv(er_input$x_path, header = F))
y <- as.matrix(utils::read.csv(er_input$y_path, row.names = 1))

deltas <- list(seq(0.0001, 0.001, 0.0001),
               seq(0.001, 0.01, 0.001),
               seq(0.01, 0.1, 0.01),
               seq(0.1, 1, 0.1))

## Step 1: Coarse Delta Search
foreach (i = 1:length(deltas)) %dopar% {
  plainER(y = y,
          x = x,
          sigma = cor(x),
          delta = deltas[[i]],
          beta_est = er_input$beta_est,
          lambda = 0.5,
          rep_cv = er_input$rep_cv,
          diagonal = er_input$diagonal,
          merge = er_input$merge,
          equal_var = er_input$equal_var,
          alpha_level = er_input$alpha_level,
          support = er_input$support,
          correction = er_input$correction,
          verbose = er_input$verbose,
          thresh_fdr = er_input$thresh_fdr)
} -> coarse_res

## Step 2: K-Fold Cross-Validation To Find Delta Magnitude
#coarse_search <- NULL
corr_bp_data <- NULL
for (i in 1:length(coarse_res)) {
  mag_delta <- coarse_res[[i]]$opt_delta
  cat("DELTA = ", mag_delta, " . . . \n")
  foreach (j = 1:er_input$nreps, .combine = rbind) %dopar% {
    temp <- essregCV(k = er_input$k,
                     x = x,
                     y = y,
                     y_factor = er_input$y_factor,
                     delta = mag_delta,
                     perm_option = "y",
                     beta_est = er_input$beta_est,
                     sel_corr = er_input$sel_corr,
                     lambda = 0.5,
                     svm = F,
                     rep_cv = er_input$rep_cv,
                     diagonal = er_input$diagonal,
                     merge = er_input$merge,
                     equal_var = er_input$equal_var,
                     alpha_level = er_input$alpha_level,
                     support = er_input$support,
                     change_all = er_input$change_all,
                     correction = er_input$correction,
                     verbose = er_input$verbose,
                     thresh_fdr = er_input$thresh_fdr)
  } -> delta_rep
  # delta_rep <- delta_rep %>%
  #   dplyr::group_by(method)
  # if (er_input$y_factor) {
  #   delta_rep <- delta_rep %>%
  #     dplyr::summarise(mean_auc = mean(as.numeric(mean_auc)),
  #                      mean_tpr = mean(as.numeric(mean_tpr)),
  #                      mean_fpr = mean(as.numeric(mean_fpr)))
  # } else if (er_input$sel_corr) {
  #   #corr_bp_data[[length(corr_bp_data) + 1]] <- delta_rep
  #   delta_rep <- delta_rep %>%
  #     dplyr::summarise(mean_corr = mean(as.numeric(spearman_corr)))
  # } else {
  #   delta_rep <- delta_rep %>%
  #     dplyr::summarise(mean_mse = mean(as.numeric(mean_mse)))
  # }
  corr_bp_data[[length(corr_bp_data) + 1]] <- delta_rep
  # coarse_search[[length(coarse_search) + 1]] <- list("delta" = mag_delta,
  #                                                    "eval" = delta_rep)
}

## select optimal delta from k-fold cross-validation
#fin_delta <- data.frame()
# if (er_input$y_factor) {
#   for (i in 1:length(coarse_search)) {
#     coarse_res <- coarse_search[[i]]
#     res <-  c(coarse_res$delta,
#               coarse_res$eval[coarse_res$eval$method == "plainER", ]$mean_auc)
#     fin_delta <- rbind(fin_delta, res)
#   }
#   colnames(fin_delta) <- c("delta", "mean_auc")
#   best_delta <- which(fin_delta$mean_auc == max(fin_delta$mean_auc))
#   best_delta <- fin_delta$delta[best_delta]
# } else if (er_input$sel_corr) {
#   ## create boxplot of replicate correlations
#   sel_corr_res <- NULL
#   for (i in 1:length(corr_bp_data)) {
#     bp_data <- corr_bp_data[[i]] %>%
#       as.data.frame() %>%
#       dplyr::filter(method == "plainER") %>%
#       dplyr::mutate(delta = coarse_search[[i]]$delta)
#     sel_corr_res <- rbind(sel_corr_res, bp_data)
#   }
#   sel_corr_res <- sel_corr_res %>%
#     dplyr::mutate(delta = as.factor(delta))
#   pdf_file <- paste0(er_input$out_path, "/delta_selection_boxplot.pdf")
#   dir.create(file.path(dirname(pdf_file)), showWarnings = F)
#   delta_boxplot <- ggplot2::ggplot(data = sel_corr_res, ggplot2::aes(x = delta, y = spearman_corr)) +
#     ggplot2::geom_boxplot()
#   ggplot2::ggsave(pdf_file, delta_boxplot)
#   best_delta <- readline(prompt = "Enter Delta: ") %>%
#     as.numeric()
# } else {
#   for (i in 1:length(coarse_search)) {
#     coarse_res <- coarse_search[[i]]
#     res <-  c(coarse_res$delta,
#               coarse_res$eval[coarse_res$eval$method == "plainER", ]$mean_mse)
#     fin_delta <- rbind(fin_delta, res)
#   }
#   colnames(fin_delta) <- c("delta", "mean_mse")
#   best_delta <- which(fin_delta$mean_mse == min(fin_delta$mean_mse))
#   best_delta <- fin_delta$delta[best_delta]
# }

## create boxplot of replicate correlations
sel_corr_res <- NULL
for (i in 1:length(corr_bp_data)) {
  bp_data <- corr_bp_data[[i]] %>%
    as.data.frame() %>%
    dplyr::filter(method == "plainER" || method == "plainER_y") %>%
    dplyr::mutate(delta = coarse_search[[i]]$delta)
  sel_corr_res <- rbind(sel_corr_res, bp_data)
}
sel_corr_res <- sel_corr_res %>%
  dplyr::mutate(delta = as.factor(delta))
pdf_file <- paste0(er_input$out_path, "/delta_selection_boxplot.pdf")
dir.create(file.path(dirname(pdf_file)), showWarnings = F)
delta_boxplot <- ggplot2::ggplot(data = sel_corr_res, ggplot2::aes(x = delta, y = spearman_corr)) +
  ggplot2::geom_boxplot()
ggplot2::ggsave(pdf_file, delta_boxplot)

## Step 3: Fine Delta Search ###################################################
d_lbd <- best_delta - best_delta / 2
d_ubd <- best_delta + best_delta / 2
delta_grid <- seq(d_lbd, d_ubd, best_delta / 100)

fine_delta_er <- plainER(y = y,
                         x = x,
                         sigma = cor(x),
                         std_y = T,
                         std_x = T,
                         delta = delta_grid,
                         beta_est = er_input$beta_est,
                         lambda = 0.1,
                         rep_cv = er_input$rep_cv,
                         diagonal = er_input$diagonal,
                         merge = er_input$merge,
                         equal_var = er_input$equal_var,
                         alpha_level = er_input$alpha_level,
                         support = er_input$support,
                         correction = er_input$correction,
                         verbose = er_input$verbose,
                         thresh_fdr = er_input$thresh_fdr)

best_delta <- fine_delta_er$opt_delta

## Step 4: Lambda Search #######################################################
lambda_search <- NULL
corr_bp_data <- NULL
for (i in 1:length(er_input$lambda)) {
  lambda <- er_input$lambda[[i]]
  cat("LAMBDA = ", lambda, " . . . \n")
  foreach (j = 1:er_input$nreps, .combine = rbind) %dopar% {
    temp <- essregCV(k = er_input$k,
                     x = x,
                     y = y,
                     y_factor = er_input$y_factor,
                     delta = best_delta,
                     perm_option = "x",
                     beta_est = er_input$beta_est,
                     sel_corr = er_input$sel_corr,
                     lambda = lambda,
                     rep_cv = er_input$rep_cv,
                     diagonal = er_input$diagonal,
                     merge = er_input$merge,
                     equal_var = er_input$equal_var,
                     alpha_level = er_input$alpha_level,
                     support = er_input$support,
                     change_all = er_input$change_all,
                     correction = er_input$correction,
                     verbose = er_input$verbose,
                     thresh_fdr = er_input$thresh_fdr)
  } -> lambda_rep
  lambda_rep <- lambda_rep %>%
    dplyr::group_by(method)
  if (er_input$y_factor) {
    lambda_rep <- lambda_rep %>%
      dplyr::summarise(mean_auc = mean(as.numeric(mean_auc)),
                       mean_tpr = mean(as.numeric(mean_tpr)),
                       mean_fpr = mean(as.numeric(mean_fpr)))
  } else if (er_input$sel_corr) {
    corr_bp_data[[length(corr_bp_data) + 1]] <- lambda_rep
    lambda_rep <- lambda_rep %>%
      dplyr::summarise(mean_corr = mean(as.numeric(spearman_corr)))
  } else {
    lambda_rep <- lambda_rep %>%
      dplyr::summarise(mean_mse = mean(as.numeric(mean_mse)))
  }
  lambda_search[[length(lambda_search) + 1]] <- list("lambda" = lambda,
                                                     "eval" = delta_rep)
}

## select optimal lambda from k-fold cross-validation
fin_lambda <- data.frame()
if (er_input$y_factor) {
  for (i in 1:length(lambda_search)) {
    lambda_res <- lambda_search[[i]]
    res <-  c(lambda_res$lambda,
              lambda_res$eval[lambda_res$eval$method == "plainER", ]$mean_auc)
    fin_lambda <- rbind(fin_lambda, res)
  }
  colnames(fin_lambda) <- c("lambda", "mean_auc")
  best_lambda <- which(fin_lambda$mean_auc == max(fin_lambda$mean_auc))
  best_lambda <- fin_lambda$lambda[best_lambda]
} else if (er_input$sel_corr) {
  sel_corr_res <- NULL
  for (i in 1:length(corr_bp_data)) {
    bp_data <- corr_bp_data[[i]] %>%
      as.data.frame() %>%
      dplyr::filter(method == "plainER") %>%
      dplyr::mutate(lambda = coarse_search[[i]]$lambda)
    sel_corr_res <- rbind(sel_corr_res, bp_data)
  }
  sel_corr_res <- sel_corr_res %>%
    dplyr::mutate(lambda = as.factor(lambda))
  pdf_file <- paste0(er_input$out_path, "/lambda_selection_boxplot.pdf")
  dir.create(file.path(dirname(pdf_file)), showWarnings = F)
  lambda_boxplot <- ggplot2::ggplot(data = sel_corr_res, ggplot2::aes(x = lambda, y = spearman_corr)) +
    ggplot2::geom_boxplot()
  ggplot2::ggsave(pdf_file, lambda_boxplot)
  best_lambda <- readline(prompt = "Enter Lambda: ") %>%
    as.numeric()
} else {
  for (i in 1:length(lambda_search)) {
    lambda_res <- lambda_search[[i]]
    res <-  c(lambda_res$lambda,
              lambda_res$eval[lambda_res$eval$method == "plainER", ]$mean_mse)
    fin_lambda <- rbind(fin_lambda, res)
  }
  colnames(fin_lambda) <- c("lambda", "mean_mse")
  best_lambda <- which(fin_lambda$mean_mse == min(fin_lambda$mean_mse))
  best_lambda <- fin_lambda$lambda[best_lambda]
}

final_results <- c("opt_delta" = best_delta[1],
                   "opt_lambda" = best_lambda[1])

