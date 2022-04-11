#' Calculate specificity and sensitivity based on group partition
#'
#' @param p total number of variables
#' @param estGroup list of estimated variable indices
#' @param trueGroup list of true variable indices
#' @return the specificity and sensitivity
#' @export

calSPandSN <- function(p, estGroup, trueGroup) {
  TP <- TN <- FP <- FN <- TPS <- TNS <- FNS <- FPS <- 0
  estTable <- trueTable <- matrix(0, nrow=p, ncol=2)
  for (k in 1:p) {
    estTable[k,] <- checkElement(k,estGroup)
    trueTable[k,] <- checkElement(k,trueGroup)
  }
  for (i in 1:(p-1)) {
    flagiEst <- estTable[i,]
    flagiTrue <- trueTable[i,]
    for (j in (i+1):p) {
      flagjEst <- estTable[j,]
      flagjTrue <- trueTable[j,]
      if (flagiTrue[1] * flagjTrue[1] > 0 && flagjTrue[1] == flagiTrue[1]) {
        if (flagiEst[1] * flagjEst[1] > 0 && flagiEst[1] == flagjEst[1]) {
          TP <- TP + 1
          if (flagiTrue[2] == flagjTrue[2]) {
            if (flagiEst[2] == flagjEst[2]) {
              TPS <- TPS + 1
            } else {
              FNS <- FNS + 1
            }
          } else {
            if (flagiEst[2] == flagjEst[2]) {
              FPS <- FPS + 1
            } else {
              TNS <- TNS + 1
            }
          }
        } else {
          FN <- FN + 1
        }
      } else {
        if (flagiEst[1] * flagjEst[1] > 0 && flagiEst[1] == flagjEst[1]) {
          FP <- FP + 1
        } else {
          TN <- TN + 1
        }
      }
    }
  }
  ifelse (TN + FP == 0, SP <- 0, SP <- TN / (TN + FP))
  ifelse (TP + FN == 0, SN <- 0, SN <- TP / (TP + FN))
  ifelse (FNS + TPS == 0, SNS <- 0, SNS <- TPS / (FNS + TPS))
  ifelse (TNS + FPS == 0, SPS <- 0, SPS <- TNS / (FPS + TNS))
  return(c(SP = SP, SN = SN, SPS = SPS, SNS = SNS))
}
