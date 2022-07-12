#' Essential Regression Pipeline - Step 5
#'
#' Run Essential Regression Pipeline with K-Fold Cross-Validation
#'     Step 5 - Run final replicates of K-Fold CV using optimal \eqn{\delta} and optimal \eqn{\lambda}.
#'              Create a boxplot and run \code{plainER()} one last time for final results.
#'
#' @importFrom magrittr '%>%'
#' @importFrom foreach '%dopar%'
#' @param yaml_path the path to a .yaml file containing all necessary parameters/arguments
#' for Essential Regression
#' @return nothing is returned, saves boxplot of cross-validation results and final ER output
#' @export

pipelineER3 <- function(yaml_path) {
  ## process arguments
  er_input <- yaml::yaml.load_file(yaml_path)
  x <- as.matrix(utils::read.csv(er_input$x_path, row.names = 1)) ## not standardized
  y <- as.matrix(utils::read.csv(er_input$y_path, row.names = 1)) ## not standardized

  dir.create(file.path(er_input$out_path), showWarnings = F, recursive = T)

  if (er_input$y_factor) {
    y <- toCont(y, er_input$y_order)
    saveRDS(y, file = paste0(er_input$out_path, "pipeline3_y_mapping.rds"))
    orig_y <- y$cat_y
    y <- y$cont_y
  }

  ##  Step 5: K-Fold Cross-Validation With Optimal Delta and Lambda  ###########
  foreach::foreach (j = 1:er_input$nreps, .combine = rbind) %dopar% {
    temp <- NULL
    while (is.null(temp)) {
      temp <- essregCV(k = er_input$k,
                       x = x,
                       y = y,
                       delta = er_input$delta,
                       perm_option = er_input$perm_option,
                       sel_corr = er_input$sel_corr,
                       y_factor = er_input$y_factor,
                       lambda = er_input$lambda,
                       out_path = er_input$out_path,
                       rep_cv = er_input$rep_cv,
                       alpha_level = er_input$alpha_level,
                       lasso = er_input$lasso,
                       pcr = er_input$pcr,
                       plsr = er_input$plsr,
                       thresh_fdr = er_input$thresh_fdr,
                       rep = j)
    }
    temp
  } -> lambda_rep
  saveRDS(lambda_rep, paste0(er_input$out_path, "pipeline_step5.RDS"))

  ## create boxplot of replicate correlations ##################################
  bp_df <- lambda_rep %>%
    dplyr::mutate(method = as.factor(method))

  if (!is.null(er_input$perm_option)) {
    bp_df <- bp_df %>%
      dplyr::mutate(perm = sub(".*_", "", method)) %>%
      dplyr::mutate(perm = ifelse(perm == method, "no_perm", paste0(perm, "_perm"))) %>%
      dplyr::mutate(method_perm = sub("*_.", "", method)) %>%
      dplyr::mutate(method = as.factor(method),
                    perm = as.factor(perm)) %>%
      dplyr::mutate(alpha = ifelse(perm == "no_perm", 1, 0.9))
  } else {
    bp_df <- bp_df %>%
      dplyr::mutate(method = as.factor(method),
                    method_perm = as.factor(method),
                    alpha = 1)
  }

  pdf_file <- paste0(er_input$out_path, "/opt_delta_lambda_violinplot.pdf")
  dir.create(file.path(dirname(pdf_file)), showWarnings = F, recursive = T)

  if (er_input$sel_corr) {
    lambda_boxplot <- ggplot2::ggplot(data = bp_df,
                                      ggplot2::aes(x = method,
                                                   y = spear_corr,
                                                   fill = method_perm,
                                                   alpha = alpha)) +
      ggplot2::geom_violin() +
      ggplot2::labs(fill = "Method") +
      ggplot2::scale_alpha(guide = 'none')
  } else if (er_input$y_factor) {
    lambda_boxplot <- ggplot2::ggplot(data = bp_df,
                                      ggplot2::aes(x = method,
                                                   y = mean_auc,
                                                   fill = method_perm,
                                                   alpha = alpha)) +
      ggplot2::geom_violin() +
      ggplot2::labs(fill = "Method") +
      ggplot2::scale_alpha(guide = 'none')
  } else {
    lambda_boxplot <- ggplot2::ggplot(data = bp_df,
                                      ggplot2::aes(x = method,
                                                   y = mean_mse,
                                                   fill = method_perm,
                                                   alpha = alpha)) +
      ggplot2::geom_violint() +
      ggplot2::labs(fill = "Method") +
      ggplot2::scale_alpha(guide = 'none')
  }

  ggplot2::ggsave(pdf_file, lambda_boxplot, width = 20, height = 15, units = "in")

  ## Final plainER Run #########################################################
  er_output <- plainER(x = x,
                       y = y,
                       sigma = NULL,
                       delta = er_input$delta,
                       lambda = er_input$lambda,
                       rep_cv = er_input$rep_cv,
                       alpha_level = er_input$alpha_level,
                       thresh_fdr = er_input$thresh_fdr,
                       out_path = er_input$out_path)

  if (is.null(er_output)) {
    cat("plainER failed --- infeasible linear program \n")
    return ()
  }

  saveRDS(er_output, paste0(er_input$out_path, "final_delta_", er_input$delta, "_lambda_", er_input$lambda, ".rds"))
}

