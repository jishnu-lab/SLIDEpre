#' Parse .yaml file for Essential Regression
#'
#' Read a .yaml file, parse arguments, and run Essential Regression.
#'
#' @param yaml_path a string path to the .yaml file
#' @param delta \eqn{\delta}, a numeric constant or vector of constants used in
#' thresholding during node purity testing
#' @param nreps number of replicates to do for cross-validation
#' @return none
#' @export

parseRun <- function(yaml_path, delta, nreps) {
  ## process arguments
  er_input <- yaml::yaml.load_file(yaml_path)
  x <- as.matrix(utils::read.csv(er_input$x_path, row.names = 1))
  y <- as.matrix(utils::read.csv(er_input$y_path, row.names = 1))

  ## cross-validation for lambda
  lambdas <- unlist(er_input$lambda)
  lambda_reps <- NULL
  lambda_df <- NULL
  for (i in 1:length(lambdas)) {
    lambda <- lambdas[i]
    foreach (j = 1:nreps) %dopar% {
      temp <- essregCV(k = er_input$k,
                       x = x,
                       y = y,
                       delta = delta,
                       lambda = lambda,
                       rep_cv = er_input$rep_cv,
                       alpha_level = er_input$alpha_level,
                       thresh_fdr = er_input$thresh_fdr)
    } -> lambda_rep
    lambda_reps[[length(lambda_reps) + 1]] <- lambda_rep
    l_df <- do.call(rbind.data.frame, lambda_rep)
    l_df$lambda <- lambda
    lambda_df <- rbind(lambda_df, l_df)
  }

  ## normalize after cross-validation because of training/validation sets
  if (er_input$process_input) {
    x <- scale(x, center = T, scale = T)
    y <- scale(y, center = T, scale = T)
  }
  ## final output
  er_output <- plainER(x = x,
                       y = y,
                       sigma = cor(x),
                       delta = delta,
                       lambda = er_input$lambda,
                       rep_cv = er_input$rep_cv,
                       alpha_level = er_input$alpha_level,
                       thresh_fdr = er_input$thresh_fdr,
                       out_path = er_input$out_path)

  saveRDS(er_output, paste0(er_input$out_path, 'er_output_', er_output$opt_delta, '.rds'))
}
