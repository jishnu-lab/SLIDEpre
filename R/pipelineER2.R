#' Essential Regression Pipeline - Steps 3, 4
#'
#' Run Essential Regression Pipeline with K-Fold Cross-Validation
#'     Step 3 - Fine Grid Search for \eqn{\delta}
#'     Step 4 - K-Fold Cross-Validation to find \eqn{\lambda}
#'
#' @importFrom magrittr '%>%'
#' @importFrom foreach '%dopar%'
#' @param yaml_path the path to a .yaml file containing all necessary parameters/arguments
#' for Essential Regression
#' @param steps an integer or string indicating which steps of the pipeline to perform: 3, 4, "all"
#' @return nothing is returned, saves boxplot of cross-validation results for user to use
#' in selecting optimal \eqn{\lambda}
#' @export

pipelineER2 <- function(yaml_path, steps = "all") {
  ## process arguments
  er_input <- yaml::yaml.load_file(yaml_path)
  x <- as.matrix(utils::read.csv(er_input$x_path, row.names = 1)) ## not standardized
  y <- as.matrix(utils::read.csv(er_input$y_path, row.names = 1)) ## not standardized

  dir.create(file.path(er_input$out_path), showWarnings = F, recursive = T)

  if (er_input$y_factor) {
    y <- toCont(y)
    saveRDS(y, file = paste0(er_input$out_path, "pipeline2_y_mapping.rds"))
    orig_y <- y$cat_y
    y <- y$cont_y
  }

  if (steps == 3) {
    ## Step 3: Fine Delta Search ###############################################
    if (length(er_input$delta) == 1) {
      d_lbd <- er_input$delta - er_input$delta / 2
      d_ubd <- er_input$delta + er_input$delta / 2
      delta_grid <- seq(d_lbd, d_ubd, er_input$delta / 100)
    } else {
      delta_grid <- er_input$delta
    }

    fine_delta_er <- plainER(y = y,
                             x = x,
                             sigma = cor(x),
                             delta = delta_grid,
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
                             thresh_fdr = er_input$thresh_fdr,
                             out_path = er_input$out_path)

    best_delta <- fine_delta_er$opt_delta
    saveRDS(fine_delta_er, file = paste0(er_input$out_path, "delta", best_delta, "_pipeline_step3.rds"))
  } else if (steps == 4) {
    ## Step 4: Lambda Search ###################################################
    corr_bp_data <- NULL
    delta <- er_input$delta
    for (i in 1:length(er_input$lambda)) {
      lambda <- er_input$lambda[[i]]
      out_path <- paste0(er_input$out_path, "lambda_", lambda, "/")
      cat("LAMBDA = ", lambda, " . . . \n")
      foreach::foreach (j = 1:er_input$nreps, .combine = rbind) %dopar% {
        if (file.exists(file = paste0(out_path, "replicate", j, "/output_table.rds"))) {
          temp <- readRDS(paste0(out_path, "replicate", j, "/output_table.rds"))
        } else {
          temp <- essregCV(k = er_input$k,
                           x = x,
                           y = y,
                           delta = delta,
                           perm_option = er_input$perm_option,
                           beta_est = er_input$beta_est,
                           sel_corr = er_input$sel_corr,
                           y_factor = er_input$y_factor,
                           lambda = lambda,
                           out_path = er_input$out_path,
                           rep_cv = er_input$rep_cv,
                           diagonal = er_input$diagonal,
                           merge = er_input$merge,
                           equal_var = er_input$equal_var,
                           alpha_level = er_input$alpha_level,
                           support = er_input$support,
                           change_all = er_input$change_all,
                           correction = er_input$correction,
                           thresh_fdr = er_input$thresh_fdr,
                           rep = j)
        }
      } -> lambda_rep
      corr_bp_data[[length(corr_bp_data) + 1]] <- list("lambda" = lambda,
                                                       "result" = lambda_rep)
    }
    saveRDS(corr_bp_data, file = paste0(er_input$out_path, "pipeline_step4.rds"))
  } else {
    ## Step 3: Fine Delta Search ###############################################
    if (length(er_input$delta) == 1) {
      d_lbd <- er_input$delta - er_input$delta / 2
      d_ubd <- er_input$delta + er_input$delta / 2
      delta_grid <- seq(d_lbd, d_ubd, er_input$delta / 100)
    } else {
      delta_grid <- er_input$delta
    }

    fine_delta_er <- plainER(y = y,
                             x = x,
                             sigma = cor(x),
                             delta = delta_grid,
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
                             thresh_fdr = er_input$thresh_fdr,
                             out_path = er_input$out_path)

    best_delta <- fine_delta_er$opt_delta

    saveRDS(fine_delta_er, file = paste0(er_input$out_path, "pipeline_step3.rds"))

    ## Step 4: Lambda Search ###################################################
    corr_bp_data <- NULL
    for (i in 1:length(er_input$lambda)) {
      lambda <- er_input$lambda[[i]]
      cat("LAMBDA = ", lambda, " . . . \n")
      out_path <- paste0(er_input$out_path, "lambda_", lambda, "/")
      foreach::foreach (j = 1:er_input$nreps, .combine = rbind) %dopar% {
        if (file.exists(file = paste0(out_path, "replicate", j, "/output_table.rds"))) {
          temp <- readRDS(paste0(out_path, "replicate", j, "/output_table.rds"))
        } else {
          temp <- essregCV(k = er_input$k,
                           x = x,
                           y = y,
                           delta = best_delta,
                           perm_option = er_input$perm_option,
                           beta_est = er_input$beta_est,
                           sel_corr = er_input$sel_corr,
                           priors = er_input$priors,
                           y_factor = er_input$y_factor,
                           lambda = lambda,
                           rep_cv = er_input$rep_cv,
                           diagonal = er_input$diagonal,
                           merge = er_input$merge,
                           equal_var = er_input$equal_var,
                           alpha_level = er_input$alpha_level,
                           support = er_input$support,
                           change_all = er_input$change_all,
                           correction = er_input$correction,
                           thresh_fdr = er_input$thresh_fdr,
                           out_path = out_path,
                           rep = j)
        }
      } -> lambda_rep
      saveRDS(lambda_rep, file = paste0(er_input$out_path, "essregCV_lambda_", lambda, ".rds"))

      ## make CV plot
      sel_corr_res <- lambda_rep %>%
        dplyr::mutate(method = as.factor(method))
      pdf_file <- paste0(er_input$out_path, "lambda_", lambda, "_boxplot.pdf")
      dir.create(file.path(dirname(pdf_file)), showWarnings = F, recursive = T)
      delta_boxplot <- ggplot2::ggplot(data = sel_corr_res,
                                       ggplot2::aes(x = method, y = spearman_corr, fill = method)) +
        ggplot2::geom_boxplot()
      ggplot2::ggsave(pdf_file, delta_boxplot)

      corr_bp_data[[length(corr_bp_data) + 1]] <- list("lambda" = lambda,
                                                       "result" = lambda_rep)
    }
    saveRDS(corr_bp_data, file = paste0(er_input$out_path, "pipeline_step4.rds"))

    ## create boxplot of replicate correlations ################################
    sel_corr_res <- NULL
    for (i in 1:length(corr_bp_data)) {
      bp_data <- corr_bp_data[[i]]
      bp_lambda <- bp_data$lambda
      bp_df <- bp_data$result %>%
        dplyr::filter(method == "plainER" | method == "plainER_y") %>%
        dplyr::mutate(lambda = bp_lambda)
      sel_corr_res <- rbind(sel_corr_res, bp_df)
    }
    sel_corr_res <- sel_corr_res %>%
      dplyr::mutate(lambda = as.factor(lambda),
                    method = as.factor(method))
    pdf_file <- paste0(er_input$out_path, "lambda_selection_boxplot.pdf")
    dir.create(file.path(dirname(pdf_file)), showWarnings = F, recursive = T)
    lambda_boxplot <- ggplot2::ggplot(data = sel_corr_res,
                                      ggplot2::aes(x = lambda, y = spearman_corr, fill = method)) +
      ggplot2::geom_boxplot()
    ggplot2::ggsave(pdf_file, lambda_boxplot)
  }
}
