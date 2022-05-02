#' Categorical To Continuous Data
#'
#' Transform categorical response vector into continuous values
#'
#' @param y a response vector of dimension \eqn{n} of categorical values
#' @return a list containing the old (categorical) \eqn{y} and the new (continuous) \eqn{y}
#' @export

toCont <- function(y) {
  uniq_vals <- unique(y)
  num_uniq <- length(uniq_vals)
  cont_vals <- seq(1, num_uniq)
  corresp_vals <- cbind(uniq_vals, cont_vals)
  new_y <- plyr::mapvalues(y, from = uniq_vals, to = cont_vals)
  return (list("cat_y" = y,
               "cont_y" = new_y))
}
