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

  ##  Step 5: K-Fold Cross-Validation With Optimal Delta and Lambda  ###########
  foreach::foreach (j = 1:er_input$nreps, .combine = rbind) %dopar% {
    temp <- essregCV(k = er_input$k,
                     x = x,
                     y = y,
                     delta = er_input$best_delta,
                     perm_option = er_input$perm_option,
                     beta_est = er_input$beta_est,
                     sel_corr = er_input$sel_corr,
                     y_factor = er_input$y_factor,
                     lambda = er_input$lambda,
                     rep_cv = er_input$rep_cv,
                     diagonal = er_input$diagonal,
                     merge = er_input$merge,
                     equal_var = er_input$equal_var,
                     alpha_level = er_input$alpha_level,
                     support = er_input$support,
                     change_all = er_input$change_all,
                     correction = er_input$correction,
                     thresh_fdr = er_input$thresh_fdr)
  } -> lambda_rep
  saveRDS(lambda_rep, paste0(er_input$out_path, "pipeline_step5.RDS"))

  ## create boxplot of replicate correlations ##################################
  bp_df <- lambda_rep %>%
    dplyr::filter(method == "plainER" | method == "plainER_y" | method == "priorER" | method == "priorER_y")
    dplyr::mutate(method = as.factor(method))
  pdf_file <- paste0(er_input$out_path, "/opt_delta_lambda_boxplot.pdf")
  dir.create(file.path(dirname(pdf_file)), showWarnings = F, recursive = T)
  lambda_boxplot <- ggplot2::ggplot(data = sel_corr_res,
                                    ggplot2::aes(x = method, y = spearman_corr, fill = method)) +
    ggplot2::geom_boxplot()
  ggplot2::ggsave(pdf_file, lambda_boxplot)

  ## Final plainER Run #########################################################
  er_output <- plainER(x = x,
                       y = y,
                       sigma = NULL,
                       delta = er_input$delta,
                       beta_est = er_input$beta_est,
                       conf_int = er_input$conf_int,
                       pred = er_input$pred,
                       lambda = er_input$lambda,
                       rep_cv = er_input$rep_cv,
                       diagonal = er_input$diagonal,
                       merge = er_input$merge,
                       equal_var = er_input$equal_var,
                       alpha_level = er_input$alpha_level,
                       support = er_input$support,
                       correction = er_input$correction,
                       thresh_fdr = er_input$thresh_fdr,
                       out_path = er_input$out_path)

  saveRDS(final_output, paste0(er_input$out_path, "final_delta", er_input$delta, "_lambda", er_input$lambda, ".rds"))
}

