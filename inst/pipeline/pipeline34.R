################################################################################
# Essential Regression Pipeline:                                               #
#     Step 3 - Fine Delta Grid Search                                          #
#     Step 4 - K-Fold Cross-Validation to find lambda                          #
################################################################################
## load libraries
library(EssReg)
library(dplyr)
library(doParallel)
library(foreach)

## FILL IN #####
yaml_path = "/Users/aerosengart/Documents/Das Lab/pipeline34.yaml"
priors = NULL
cores = 6
################

## set up parallelization
doParallel::registerDoParallel(cores = cores)

## process arguments
er_input <- yaml::yaml.load_file(yaml_path)
x <- as.matrix(utils::read.csv(er_input$x_path, header = F)) ## not standardized
y <- as.matrix(utils::read.csv(er_input$y_path, row.names = 1)) ## not standardized

## Step 3: Fine Delta Search ###################################################
if (is.null(delta_grid)) {
  d_lbd <- best_delta - best_delta / 2
  d_ubd <- best_delta + best_delta / 2
  delta_grid <- seq(d_lbd, d_ubd, best_delta / 100)
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
                         thresh_fdr = er_input$thresh_fdr)

best_delta <- fine_delta_er$opt_delta

## Step 4: Lambda Search #######################################################
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
                     perm_option = er_input$perm_option,
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
  corr_bp_data[[length(corr_bp_data) + 1]] <- list("lambda" = lambda,
                                                   "result" = lambda_rep)
}

## create boxplot of replicate correlations
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
  dplyr::mutate(delta = as.factor(lambda),
                method = as.factor(method))
pdf_file <- paste0(er_input$out_path, "/lambda_selection_boxplot.pdf")
dir.create(file.path(dirname(pdf_file)), showWarnings = F)
lambda_boxplot <- ggplot2::ggplot(data = sel_corr_res,
                                 ggplot2::aes(x = delta, y = spearman_corr, fill = method)) +
  ggplot2::geom_boxplot()
ggplot2::ggsave(pdf_file, lambda_boxplot)

