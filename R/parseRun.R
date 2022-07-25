#' Parse .yaml file for Essential Regression
#'
#' Read a .yaml file, parse arguments, and run Essential Regression.
#'
#' @param yaml_path a string path to the .yaml file
#' @return none
#' @export

parseRun <- function(yaml_path) {
  ## process arguments
  er_input <- yaml::yaml.load_file(yaml_path)
  x <- as.matrix(utils::read.csv(er_input$x_path, row.names = 1))
  y <- as.matrix(utils::read.csv(er_input$y_path, row.names = 1))

  dir.create(file.path(er_input$out_path), showWarnings = F, recursive = T)

  ## make continuous if categorical
  if (er_input$y_factor) {
    y <- toCont(y, er_input$y_order)
    saveRDS(y, file = paste0(er_input$out_path, "plainER_y_mapping.rds"))
    orig_y <- y$cat_y
    y <- y$cont_y
  }

  ## final output
  er_output <- plainER(x = x,
                       y = y,
                       sigma = cor(x),
                       delta = er_input$delta,
                       lambda = er_input$lambda,
                       rep_cv = er_input$rep_cv,
                       alpha_level = er_input$alpha_level,
                       thresh_fdr = er_input$thresh_fdr,
                       out_path = er_input$out_path)

  saveRDS(er_output, paste0(er_input$out_path, 'er_output_', er_output$opt_delta, '.rds'))
}
