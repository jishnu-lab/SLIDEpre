#' Run Shiny App
#'
#' Run Shiny App
#'
#' @export

runShiny <- function() {
  appDir <- system.file("shiny", package = "EssReg")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
