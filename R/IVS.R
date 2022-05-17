#' Iterative Variable Selection
#'
#' A feature selection method suitable for variables with high correlation
#'
#' @importFrom magrittr '%>%'
#' @param y a response vector of dimension \eqn{n}
#' @param z a matrix of dimensions \eqn{p \times K}
#' @return a vector of selected variable indices
#' @export

IVS <- function(y, z, verbose = F){
  if (is.null(y) | is.null(z)){
    stop("y and z must be given")
  }
  pvalueVec  <- NULL

  for(i in 1:ncol(z)){
    pvalueVec <- rbind(pvalueVec,summary(lm(y~z[,i]))$coef[2,"Pr(>|t|)"])
  }

  ii <- which(pvalueVec<0.1)
  z <- z[,ii]

  ## Select the first variable
  minPvalI <- which.min(pvalueVec[ii])
  S <- minPvalI
  AdR_old <- 0
  Ad_new  <- 1

  while (Ad_new > AdR_old){
    Obj <- NULL
    for (i in c(1:ncol(z))[-S]){
      r_adjusted <- summary(lm(y ~ z[, cbind(S, i)]))$adj.r.squared
      colinear   <- summary(lm(z[, i] ~ z[, S]))$r.squared
      Obj <- rbind(Obj, cbind(i, r_adjusted - colinear))
    }

    S <- c(S,unname(Obj[which.max(Obj[,2]),1]))
    ifelse(length(S)==1,AdR_old <- summary(lm(y~z[,S]))$adj.r.squared,AdR_old <- summary(lm(y~z[,S[-length(S)]]))$adj.r.squared)
    Ad_new <- summary(lm(y~z[,S]))$adj.r.squared

    if (verbose == T) {
      print(paste("old Adj R2 is:",AdR_old))
      print(paste("new Adj R2 is:",Ad_new))
    }
  }

  a <- length(S)
  print(a)
  print(S)
  return(ii[S[-a]])
}
