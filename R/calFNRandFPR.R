#' Calculate the FNR and FPR for the estimated matrix, \eqn{A}.
#'
#' @param estA an estimated matrix
#' @param trueA the true matrix of the same dimensions as \code{estA}
#' @return a list containing the FNRs and FPRs as well as the row indices of where the errors occur
#' @export

calFNRandFPR <- function(estA,trueA) {
  if (sum(dim(estA) == dim(trueA)) != 2) {
    # If estA and A don't have the same dimension, we return -1.
    return(list(error    = c(FPR = -1, FNR = -1, FPSR = -1, FNSR = -1),
                errorInd = c()))
  }
  errorInd <- c()
  p <- nrow(estA)
  K <- ncol(estA)
  TP <- FP <- TN <- FN <- FPS <- TPS <- TNS <- FNS <- 0
  for (i in 1:p) {
    for (j in 1:K) {
      if (trueA[i,j] != 0) {
        TN <- TN + 1
        if (estA[i,j] == 0) {
          errorInd <- c(errorInd, i) # record the row indices where errors occur
          FN <- FN + 1
        } else {  ### calculate the sign error rate
          if (trueA[i,j] > 0) {
            TPS <- TPS + 1
            if (estA[i,j] < 0)
              FNS <- FNS + 1
          } else {
            TNS <- TNS + 1
            if (estA[i,j] > 0)
              FPS <- FPS + 1
          }
        }
      } else {
        TP <- TP + 1
        if (estA[i,j] != 0) {
          FP <- FP + 1
          errorInd <- c(errorInd, i)
        }
      }
    }
  }
  if (TPS == 0)
    FNSR <- 0
  if (TNS == 0)
    FPSR <- 0
  return(list(error = c(FPR = FP / TP, FNR = FN / TN,FPSR = FPS / TNS, FNSR = FNS / TPS),
              errorInd = unique(errorInd)))
}
