#' Calculate the error rates for an estimated matrix.
#'
#' @param estA an estimated matrix
#' @param trueA the true matrix of same dimensions as \code{estA}
#' @param flagPerm boolean for whether to use the optimal permuted version of \code{estA} for the calculations
#' @return the error rates for \code{estA}: specificity, sensitivity, false negative rate, and false positive rate
#' @export

getErrors <- function(estA, trueA, flagPerm = TRUE) {
  estGroup <- recoverGroup(estA)
  trueGroup <- recoverGroup(trueA)
  speAndSen <- calSPandSN(nrow(trueA), estGroup, trueGroup)
  if (flagPerm) {
    estAperm <- permA(estA, trueA)
  } else {
    estAperm <- estA
  }
  fr <- calFNRandFPR(estAperm, trueA)$error
  if (fr[1] == -1) {
    rates <- rep(-1, 3)
  } else {
    rates <- calArates(estAperm - trueA)
  }
  return(list(SPandSN = round(speAndSen,4), FNRandFPR = round(fr,4),rates = round(rates,4)))
}
