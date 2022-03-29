#' Cross validation for choosing \eqn{\delta}. For each delta from the given grids,
#' first split the data into two data sets. Obtain \eqn{I}, \eqn{A_I} and \eqn{C} from data set 1.
#' Then calculate \eqn{A_I \cdot C \cdot A_I^\top} and choose \eqn{\delta} which minimizes the criterion
#' \eqn{A_I \cdot C \cdot A_I^\top - \Sigma(\text{data set 2})}
#'
#' @param x data matrix of dimensions \eqn{n \times p}
#' @param deltaGrids a vector of numerical constants over which to perform the search for the optimal \eqn{\delta}
#' @param diagonal a boolean indicating the diagonal structure of \eqn{C}
#' @param se_est estimated standard errors
#' @param merge a boolean indicating merge style
#' @return the selected optimal \eqn{\delta}

CV_Delta <- function(x, deltaGrids, diagonal, se_est, merge) {
  n <- nrow(x); p <- ncol(x)
  sampInd <- sample(n, floor(n / 2))
  X1 <- x[sampInd, ]
  X2 <- x[-sampInd, ]
  Sigma1 <- crossprod(X1) / nrow(X1);
  diag(Sigma1) <- 0
  Sigma2 <- crossprod(X2) / nrow(X2)

  result_Ms <- FindRowMax(abs(Sigma1))
  Ms <- result_Ms$M
  arg_Ms <- result_Ms$arg_M

  loss <- c()
  for (i in 1:length(deltaGrids)) {
    resultFitted <- CalFittedSigma(Sigma1, deltaGrids[i], Ms, arg_Ms, se_est, diagonal, merge)
    fittedValue <- resultFitted$fitted
    estPureVec <- resultFitted$pureVec
    if (is.null(dim(fittedValue)) && fittedValue == -1) {
      loss[i] <- Inf
    } else {
      denom <- length(estPureVec) * (length(estPureVec) - 1)
      loss[i] <- 2 * offSum(Sigma2[estPureVec, estPureVec], fittedValue, se_est[estPureVec]) / denom
    }
  }
  return(deltaGrids[which.min(loss)])
}
