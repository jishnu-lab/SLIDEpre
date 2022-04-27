#' Essential Regression Pipeline - Steps 1, 2
#'
#' Run Essential Regression Pipeline with K-Fold Cross-Validation
#'     Step 1 - Coarse Grid Search for \eqn{\delta}
#'     Step 2 - K-Fold Cross-Validation to find magnitude of \eqn{\delta}
#'
#' @importFrom magrittr '%>%'
#' @importFrom foreach '%dopar%'
#' @param yaml_path the path to a .yaml file containing all necessary parameters/arguments
#' for Essential Regression
#' @return nothing is returned, saves boxplot of cross-validation results for user to use
#' in selecting optimal \eqn{\delta}
#' @export

pipelineER1 <- function(yaml_path) {
  ## process arguments
  er_input <- yaml::yaml.load_file(yaml_path)
  x <- as.matrix(utils::read.csv(er_input$x_path, row.names = 1)) ## not standardized
  y <- as.matrix(utils::read.csv(er_input$y_path, row.names = 1)) ## not standardized

  deltas <- list(seq(0.0001, 0.001, 0.0001),
                 seq(0.001, 0.01, 0.001),
                 seq(0.01, 0.1, 0.01),
                 seq(0.1, 1, 0.1))

  ## Step 1: Coarse Delta Search
  foreach::foreach (i = 1:length(deltas)) %dopar% {
    plainER(y = y,
            x = x,
            sigma = NULL,
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
  corr_bp_data <- NULL
  for (i in 1:length(coarse_res)) {
    mag_delta <- coarse_res[[i]]$opt_delta
    cat("DELTA = ", mag_delta, " . . . \n")
    foreach::foreach (j = 1:er_input$nreps, .combine = rbind) %dopar% {
      temp <- essregCV(k = er_input$k,
                       x = x,
                       y = y,
                       y_factor = er_input$y_factor,
                       delta = mag_delta,
                       perm_option = er_input$perm_option,
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
    corr_bp_data[[length(corr_bp_data) + 1]] <- list("delta" = mag_delta,
                                                     "result" = delta_rep)
  }

  ## create boxplot of replicate correlations
  sel_corr_res <- NULL
  for (i in 1:length(corr_bp_data)) {
    bp_data <- corr_bp_data[[i]]
    bp_delta <- bp_data$delta
    bp_df <- bp_data$result %>%
      dplyr::filter(method == "plainER" | method == "plainER_y") %>%
      dplyr::mutate(delta = bp_delta)
    sel_corr_res <- rbind(sel_corr_res, bp_df)
  }
  sel_corr_res <- sel_corr_res %>%
    dplyr::mutate(delta = as.factor(delta),
                  method = as.factor(method))
  pdf_file <- paste0(er_input$out_path, "/delta_selection_boxplot.pdf")
  dir.create(file.path(dirname(pdf_file)), showWarnings = F)
  delta_boxplot <- ggplot2::ggplot(data = sel_corr_res,
                                   ggplot2::aes(x = delta, y = spearman_corr, fill = method)) +
    ggplot2::geom_boxplot()
  ggplot2::ggsave(pdf_file, delta_boxplot)
}

