#' Iterative Variable Selection - Robustified
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

IVSRobust <- function(y, z, freq_tab) {
  pvalueVec <- NULL
  for(i in 1:ncol(z)) {
    pvalueVec <- rbind(pvalueVec, summary(lm(y ~ z[, i]))$coef[2, "Pr(>|t|)"])
  }

  p_val_df <- pvalueVec %>%
    as.data.frame() %>%
    dplyr::mutate(z = 1:length(pvalueVec))

  freq_tab <- freq_tab %>%
    as.data.frame() %>%
    dplyr::mutate(z = as.numeric(row.names(freq_tab)))

  df <- dplyr::full_join(p_val_df, freq_tab, by = "z")
  df[, 3] <- ifelse(is.na(df[, 3]), 0, df[, 3])
  colnames(df) <- c("p_val", "z", "freq")
  df <- df %>%
    dplyr::select(z, freq, p_val) %>%
    dplyr::arrange(-freq, p_val) %>%
    dplyr::filter(freq != 0 | p_val < 0.1)

  z <- as.matrix(z[, df[, "z"]])

  minPvalI <- which.max(df[, "freq"])
  S <- minPvalI
  a <- length(S)
  AdR_old <- 0
  Ad_new  <- 1

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
  }
  return(df[S[-a], "z"])
}
