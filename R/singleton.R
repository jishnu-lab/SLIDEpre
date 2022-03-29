#' Check if there exists an element of the given list has length equal to 1
#'
#' @param estPureIndices estimated indices of pure nodes
#' @return true or false depending upon existence of singleton element

singleton <- function(estPureIndices) {
  if (length(estPureIndices) == 0) {
    return(T)
  } else {
    ifelse(sum(sapply(estPureIndices, FUN = function(x) {length(x)}) == 1) > 0, T, F)
  }
}
