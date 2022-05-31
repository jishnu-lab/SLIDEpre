#' Iterative Variable Selection
#'
#' A feature selection method suitable for variables with high correlation
#'
#' @importFrom magrittr '%>%'
#' @param y a response vector of dimension \eqn{n}
#' @param z a matrix of dimensions \eqn{p \times K}
#' @param priors a vector of indices indicating the important features
#' @param er_res an object returned by \code{plainER} or \code{priorER}
#' @return a vector of selected variable indices
#' @export

IVS <- function(y, z, priors = NULL, er_res = NULL, verbose = F) {
  if (is.null(y) | is.null(z)) {
    stop("y and z must be given")
  }
  pvalueVec <- NULL

  for(i in 1:ncol(z)) {
    pvalueVec <- rbind(pvalueVec, summary(lm(y ~ z[, i]))$coef[2, "Pr(>|t|)"])
  }

  if (!is.null(priors)) {
    #### find important features in each cluster
    er_read <- readER(er_res)
    clust_feats <- list()
    imp_clusts <- NULL
    for (i in 1:er_res$K) {
      cluster <- unlist(er_read$clusters[[i]])
      imp_feats_cluster <- intersect(cluster, priors)
      clust_feats[[length(clust_feats) + 1]] <- imp_feats_cluster
      if (length(imp_feats_cluster) > 0) {
        imp_clusts <- c(imp_clusts, i)
      }
    }
    ii <- c(which(pvalueVec < 0.1), imp_clusts) ## keep clusters with low enough p-value and imp_z
  } else {
    ii <- which(pvalueVec < 0.1)
  }

  if (length(ii) < 3) { ## if too few sig variables
    ii <- sort(pvalueVec, decreasing = FALSE, index.return=TRUE)
    ii <- ii$ix[1:3] ## use lowest 3
  }
  ii <- unique(ii)
  z <- as.matrix(z[, ii])

  ## Select the first variable
  minPvalI <- which.min(pvalueVec[ii])
  S <- minPvalI
  a <- length(S)
  AdR_old <- 0
  Ad_new  <- 1

  if (!is.null(priors)) {
    ## make weight for R^2
    loadings <- er_res$A
    z_imp_probs <- rep(0, ncol(z))
    for (i in 1:length(ii)) {
      column <- loadings[, ii[i]] ## get one column of loadings matrix A
      imp_col <- column[priors] ## get just rows of important features
      num_nonzero <- length(which(imp_col > 0)) ## get number of nonzero important features
      z_imp_probs[i] <- 1 + (num_nonzero / length(which(column != 0))) ## 1 + proportion of nonzero features are imp
    }

    while (abs((abs(Ad_new) - abs(AdR_old)) > 0.01) && (length(y) > (length(S) + 1)) && (length(S) < ncol(z))) {
      Obj <- NULL ## objective function = adjusted R^2 - collinearity R^2
      for (i in c(1:ncol(z))[-S]) {
        r_adjusted <- summary(lm(y ~ z[, cbind(S, i)]))$adj.r.squared
        colinear   <- summary(lm(z[, i] ~ z[, S]))$r.squared
        Obj <- rbind(Obj, cbind(i, (r_adjusted - colinear)))
      }

      S <- c(S, unname(Obj[which.max(Obj[, 2]), 1]))
      ifelse(length(S) == 1,
             AdR_old <- summary(lm(y ~ z[, S]))$adj.r.squared,
             AdR_old <- summary(lm(y ~ z[, S[-length(S)]]))$adj.r.squared)
      Ad_new <- summary(lm(y ~ z[,S]))$adj.r.squared
      a <- length(S)
      if (verbose == T) {
        print(paste("old Adj R2 is:", AdR_old))
        print(paste("new Adj R2 is:", Ad_new))
        print(ii[S[-a]])
      }
    }
  } else {
    while (abs((abs(Ad_new) - abs(AdR_old)) > 0.01) && (length(y) > (length(S) + 1)) && (length(S) < ncol(z))) {
      Obj <- NULL ## objective function = adjusted R^2 - colinearity R^2
      for (i in c(1:ncol(z))[-S]) {
        r_adjusted <- summary(lm(y ~ z[, cbind(S, i)]))$adj.r.squared
        colinear   <- summary(lm(z[, i] ~ z[, S]))$r.squared
        Obj <- rbind(Obj, cbind(i, (r_adjusted - colinear)))
      }

      S <- c(S, unname(Obj[which.max(Obj[, 2]), 1]))
      ifelse(length(S) == 1,
             AdR_old <- summary(lm(y ~ z[, S]))$adj.r.squared,
             AdR_old <- summary(lm(y ~ z[, S[-length(S)]]))$adj.r.squared)
      Ad_new <- summary(lm(y ~ z[,S]))$adj.r.squared
      a <- length(S)
      if (verbose == T) {
        print(paste("old Adj R2 is:", AdR_old))
        print(paste("new Adj R2 is:", Ad_new))
        print(ii[S[-a]])
      }
    }
  }

  return(ii[S[-a]])
}
