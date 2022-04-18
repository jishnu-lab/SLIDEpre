#' Parse .yaml file for Essential Regression
#'
#' Read a .yaml file, parse arguments, and run Essential Regression.
#'
#' @param yaml_path a string path to the .yaml file
#' @param delta \eqn{\delta}, a numeric constant or vector of constants used in
#' thresholding during node purity testing
#' @return none
#' @export

parseRun <- function(yaml_path, delta) {
  er_input <- yaml::yaml.load_file(yaml_path)
  x <- as.matrix(utils::read.csv(er_input$x_path, row.names = 1))
  y <- as.matrix(utils::read.csv(er_input$y_path, row.names = 1))

  if (er_input$process_input) {
    x <- scale(x, center = T, scale = T)
    y <- scale(y, center = T, scale = T)
  }

  er_output <- plainER(x = x,
                       y = y,
                       sigma = cor(x),
                       delta = delta,
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
                       verbose = er_input$verbose,
                       thresh_fdr = er_input$thresh_fdr,
                       out_path = er_input$out_path)

  saveRDS(er_output, paste0(er_input$out_path, 'er_output_', er_output$opt_delta, '.rds'))
}
