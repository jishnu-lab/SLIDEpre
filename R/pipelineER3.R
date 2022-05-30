#' Essential Regression Pipeline - Step 5
#'
#' Run Essential Regression Pipeline with K-Fold Cross-Validation
#'     Step 5 - After selecting optimal \eqn{\delta} and optimal \eqn{\lambda}
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

  #####################################################################
  ##  Step 5: K-Fold Cross-Validation With Optimal Delta and Lambda  ##
  #####################################################################
  corr_bp_data <- NULL
  for (i in 1:length(er_input$lambda)) {
    lambda <- er_input$lambda[[i]]
    cat("LAMBDA = ", lambda, " . . . \n")
    foreach::foreach (j = 1:er_input$nreps, .combine = rbind) %dopar% {
      temp <- essregCV(k = er_input$k,
                       x = x,
                       y = y,
                       delta = er_input$best_delta,
                       perm_option = er_input$perm_option,
                       beta_est = er_input$beta_est,
                       sel_corr = er_input$sel_corr,
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
                       verbose = er_input$verbose,
                       thresh_fdr = er_input$thresh_fdr)
    } -> lambda_rep
    corr_bp_data[[length(corr_bp_data) + 1]] <- list("lambda" = lambda,
                                                     "result" = lambda_rep)
  }
  saveRDS(corr_bp_data, paste0(er_input$out_path, "corr_bp_data_s5.RDS"))

  ##################################################################
  ##           create boxplot of replicate correlations           ##
  ##################################################################
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
  pdf_file <- paste0(er_input$out_path, "/opt_delta_lambda_boxplot.pdf")
  dir.create(file.path(dirname(pdf_file)), showWarnings = F, recursive = T)
  lambda_boxplot <- ggplot2::ggplot(data = sel_corr_res,
                                    ggplot2::aes(x = lambda, y = spearman_corr, fill = method)) +
    ggplot2::geom_boxplot()
  ggplot2::ggsave(pdf_file, lambda_boxplot)

  ##################################################################
  ##                  run plainER to save output                  ##
  ##################################################################
  final_output <- vector(mode = "list", length = length(er_input$lambda))
  for (i in 1:length(er_input$lambda)) {
    lambda <- er_input$lambda[[i]]
    er_output <- plainER(x = x,
                         y = y,
                         sigma = NULL,
                         delta = er_input$best_delta,
                         beta_est = er_input$beta_est,
                         conf_int = er_input$conf_int,
                         pred = er_input$pred,
                         lambda = lambda,
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

    final_output[[i]] <- er_output
  }
  saveRDS(final_output, paste0(er_input$out_path, "final_er_output.RDS"))
}

