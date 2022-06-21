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
#' @param steps an integer or string indicating which steps of the pipeline to perform: 1, 2, "all"
#' @return nothing is returned, saves boxplot of cross-validation results for user to use
#' in selecting optimal \eqn{\delta}
#' @export

pipelineER1 <- function(yaml_path, steps = "all") {
  ## process arguments
  er_input <- yaml::yaml.load_file(yaml_path)
  x <- as.matrix(utils::read.csv(er_input$x_path, row.names = 1)) ## not standardized
  y <- as.matrix(utils::read.csv(er_input$y_path, row.names = 1)) ## not standardized

  dir.create(file.path(er_input$out_path), showWarnings = F, recursive = T)

  if (er_input$y_factor) {
    y <- toCont(y, er_input$y_order)
    saveRDS(y, file = paste0(er_input$out_path, "pipeline1_y_mapping.rds"))
    orig_y <- y$cat_y
    y <- y$cont_y
  }

  deltas <- list(seq(0.0001, 0.001, 0.0001),
                 seq(0.001, 0.01, 0.001),
                 seq(0.01, 0.1, 0.01),
                 seq(0.1, 1, 0.1))

  if (steps == 1) {
    ## Step 1: Coarse Delta Search #############################################
    if (file.exists(paste0(er_input$out_path, "pipeline_step1.rds"))) {
      coarse_res <- readRDS(paste0(er_input$out_path, "pipeline_step1.rds"))
    } else {
      foreach::foreach (i = 1:length(deltas)) %dopar% {
        result <- plainER(y = y,
                          x = x,
                          sigma = NULL,
                          delta = deltas[[i]],
                          lambda = 0.5,
                          rep_cv = er_input$rep_cv,
                          alpha_level = er_input$alpha_level,
                          thresh_fdr = er_input$thresh_fdr,
                          out_path = er_input$out_path)
        if (is.null(result)) {
          cat("plainER failed --- infeasible linear program \n")
          return ()
        }
        result
      } -> coarse_res
      saveRDS(coarse_res, file = paste0(er_input$out_path, "pipeline_step1.rds"))
    }
  } else if (steps == 2) {
    ## Step 2: K-Fold Cross-Validation To Find Delta Magnitude #################
    mag_delta <- er_input$delta
    if (file.exists(paste0(er_input$out_path, "delta_", mag_delta, "pipeline_step2.rds"))) {
      delta_rep <- readRDS(paste0(er_input$out_path, "delta_", mag_delta, "pipeline_step2.rds"))
    } else {
      foreach::foreach (j = 1:er_input$nreps, .combine = rbind) %dopar% {
        if (file.exists(file = paste0(er_input$out_path, "replicate", j, "/output_table.rds"))) {
          readRDS(paste0(er_input$out_path, "replicate", j, "/output_table.rds"))
        } else {
          result <- NULL
          while (is.null(result)) {
            result <- essregCV(k = er_input$k,
                               x = x,
                               y = y,
                               delta = mag_delta,
                               y_factor = er_input$y_factor,
                               perm_option = er_input$perm_option,
                               sel_corr = er_input$sel_corr,
                               lambda = 0.5,
                               rep_cv = er_input$rep_cv,
                               alpha_level = er_input$alpha_level,
                               thresh_fdr = er_input$thresh_fdr,
                               out_path = paste0(er_input$out_path, "delta_", mag_delta, "/"),
                               lasso = er_input$lasso,
                               pcr = er_input$pcr,
                               plsr = er_input$plsr,
                               rep = j)
          }
          result
        }
      } -> delta_rep
      saveRDS(delta_rep, file = paste0(er_input$out_path, "delta_", mag_delta, "pipeline_step2.rds"))

      ## make CV plot
      if (!is.null(er_input$perm_option)) {
        final_res <- delta_rep %>%
          dplyr::mutate(perm = sub(".*_", "", method)) %>%
          dplyr::mutate(perm = ifelse(perm == er_input$perm_option, perm, "no_perm")) %>%
          dplyr::mutate(method = as.factor(method),
                        perm = as.factor(perm))
      } else {
        final_res <- delta_rep %>%
          dplyr::mutate(method = as.factor(method))
      }

      pdf_file <- paste0(er_input$out_path, "delta_", mag_delta, "_boxplot.pdf")
      dir.create(file.path(dirname(pdf_file)), showWarnings = F, recursive = T)
      if (er_input$sel_corr) {
        delta_boxplot <- ggplot2::ggplot(data = final_res,
                                         ggplot2::aes(x = method,
                                                      y = spear_corr,
                                                      fill = ifelse(!is.null(er_input$perm_option),
                                                                    perm, method))) +
          ggplot2::geom_boxplot()
      } else if (er_input$y_factor) {
        delta_boxplot <- ggplot2::ggplot(data = final_res,
                                         ggplot2::aes(x = method,
                                                      y = mean_auc,
                                                      fill = ifelse(!is.null(er_input$perm_option),
                                                                    perm, method))) +
          ggplot2::geom_boxplot()
      } else {
        delta_boxplot <- ggplot2::ggplot(data = final_res,
                                         ggplot2::aes(x = method,
                                                      y = mean_mse,
                                                      fill = ifelse(!is.null(er_input$perm_option),
                                                                    perm, method))) +
          ggplot2::geom_boxplot()
      }
      ggplot2::ggsave(pdf_file, delta_boxplot)
    }
  } else {
    ## Step 1: Coarse Delta Search #############################################
    if (file.exists(paste0(er_input$out_path, "pipeline_step1.rds"))) {
      coarse_res <- readRDS(file = paste0(er_input$out_path, "pipeline_step1.rds"))
    } else {
      foreach::foreach (i = 1:length(deltas)) %dopar% {
        result <- plainER(y = y,
                          x = x,
                          sigma = NULL,
                          delta = deltas[[i]],
                          lambda = 0.5,
                          rep_cv = er_input$rep_cv,
                          alpha_level = er_input$alpha_level,
                          thresh_fdr = er_input$thresh_fdr,
                          out_path = er_input$out_path)
        if (is.null(result)) {
          cat("plainER failed --- infeasible linear program \n")
          return ()
        }
        result
      } -> coarse_res
      saveRDS(coarse_res, file = paste0(er_input$out_path, "pipeline_step1.rds"))
    }

    ## Step 2: K-Fold Cross-Validation To Find Delta Magnitude #################
    corr_bp_data <- NULL
    for (i in 1:length(coarse_res)) {
      mag_delta <- coarse_res[[i]]$opt_delta
      magnitude <- deltas[[i]][1]
      cat("DELTA = ", mag_delta, " . . . \n")
      if (file.exists(paste0(er_input$out_path, "essregCV_delta_", mag_delta, ".rds"))) {
        delta_rep <- readRDS(paste0(er_input$out_path, "essregCV_delta_", mag_delta, ".rds"))
      } else {
        foreach::foreach (j = 1:er_input$nreps, .combine = rbind) %dopar% {
          if (file.exists(file = paste0(er_input$out_path, "delta_", mag_delta, "/replicate", j, "/output_table.rds"))) {
            readRDS(paste0(er_input$out_path, "delta", mag_delta, "/replicate", j, "/output_table.rds"))
          } else {
            result <- NULL
            while (is.null(result)) {
              result <- essregCV(k = er_input$k,
                                 x = x,
                                 y = y,
                                 delta = mag_delta,
                                 y_factor = er_input$y_factor,
                                 perm_option = er_input$perm_option,
                                 sel_corr = er_input$sel_corr,
                                 lambda = 0.5,
                                 rep_cv = er_input$rep_cv,
                                 alpha_level = er_input$alpha_level,
                                 thresh_fdr = er_input$thresh_fdr,
                                 out_path = paste0(er_input$out_path, "delta_", mag_delta, "/"),
                                 lasso = er_input$lasso,
                                 pcr = er_input$pcr,
                                 plsr = er_input$plsr,
                                 rep = j)
            }
            result
          }
        } -> delta_rep
        saveRDS(delta_rep, file = paste0(er_input$out_path, "essregCV_delta_", mag_delta, ".rds"))
      }
      ## make CV plot
      if (!is.null(er_input$perm_option)) {
        final_res <- delta_rep %>%
          dplyr::mutate(perm = sub(".*_", "", method)) %>%
          dplyr::mutate(perm = ifelse(perm == er_input$perm_option, perm, "no_perm")) %>%
          dplyr::mutate(method = as.factor(method),
                        perm = as.factor(perm))
        pdf_file <- paste0(er_input$out_path, "delta_", mag_delta, "_boxplot.pdf")
        dir.create(file.path(dirname(pdf_file)), showWarnings = F, recursive = T)
        if (er_input$sel_corr) {
          delta_boxplot <- ggplot2::ggplot(data = final_res,
                                           ggplot2::aes(x = method, y = spear_corr, fill = perm)) +
            ggplot2::geom_boxplot()
        } else if (er_input$y_factor) {
          delta_boxplot <- ggplot2::ggplot(data = final_res,
                                           ggplot2::aes(x = method, y = mean_auc, fill = perm)) +
            ggplot2::geom_boxplot()
        } else {
          delta_boxplot <- ggplot2::ggplot(data = final_res,
                                           ggplot2::aes(x = method, y = mean_mse, fill = perm)) +
            ggplot2::geom_boxplot()
        }

        ggplot2::ggsave(pdf_file, delta_boxplot, width = 20, height = 15, units = "in")
      } else {
        eval_res <- delta_rep %>%
          dplyr::mutate(method = as.factor(method))
        pdf_file <- paste0(er_input$out_path, "delta_", mag_delta, "_boxplot.pdf")
        dir.create(file.path(dirname(pdf_file)), showWarnings = F, recursive = T)
        if (er_input$sel_corr) {
          delta_boxplot <- ggplot2::ggplot(data = evall_res,
                                           ggplot2::aes(x = method, y = spear_corr, fill = perm)) +
            ggplot2::geom_boxplot()
        } else if (er_input$y_factor) {
          delta_boxplot <- ggplot2::ggplot(data = eval_res,
                                           ggplot2::aes(x = method, y = mean_auc, fill = method)) +
            ggplot2::geom_boxplot()
        } else {
          delta_boxplot <- ggplot2::ggplot(data = eval_res,
                                           ggplot2::aes(x = method, y = mean_mse, fill = method)) +
            ggplot2::geom_boxplot()
        }
        ggplot2::ggsave(pdf_file, delta_boxplot, width = 20, height = 15, units = "in")
      }

      corr_bp_data[[length(corr_bp_data) + 1]] <- list("delta" = mag_delta,
                                                       "result" = delta_rep)
    }
    saveRDS(corr_bp_data, file = paste0(er_input$out_path, "pipeline_step2.rds"))
  }

  ## create boxplot of replicate correlations ##################################
  final_res <- NULL
  for (i in 1:length(corr_bp_data)) {
    bp_data <- corr_bp_data[[i]]
    bp_delta <- bp_data$delta
    bp_df <- bp_data$result %>%
      dplyr::filter(method == "plainER" | method == "plainER_y") %>%
      dplyr::mutate(delta = bp_delta)
    final_res <- rbind(final_res, bp_df)
  }
  final_res <- final_res %>%
    dplyr::mutate(delta = as.factor(delta),
                  method = as.factor(method))
  pdf_file <- paste0(er_input$out_path, "/delta_selection_boxplot.pdf")
  dir.create(file.path(dirname(pdf_file)), showWarnings = F, recursive = T)
  if (er_input$sel_corr) {
    delta_boxplot <- ggplot2::ggplot(data = final_res,
                                     ggplot2::aes(x = delta, y = spear_corr, fill = method)) +
      ggplot2::geom_boxplot()
  } else if (er_input$y_factor) {
    delta_boxplot <- ggplot2::ggplot(data = final_res,
                                     ggplot2::aes(x = delta, y = mean_auc, fill = method)) +
      ggplot2::geom_boxplot()
  } else {
    delta_boxplot <- ggplot2::ggplot(data = final_res,
                                     ggplot2::aes(x = delta, y = mean_mse, fill = method)) +
      ggplot2::geom_boxplot()
  }

  ggplot2::ggsave(pdf_file, delta_boxplot, width = 20, height = 15, units = "in")
}
