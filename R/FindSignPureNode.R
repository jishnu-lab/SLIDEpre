#' Estimate the sign subpartition of pure node sets. If there is an element
#' of a list is empty, then a empty list will be put in that position.
#'
#' @param pureList a list of pure node indices
#' @param Sigma a sample correlation matrix of dimensions \eqn{p \times p}
#' @return a list of sign subpartitions of pure node indices

FindSignPureNode <- function(pureList, Sigma) {
  signPureList <- list()
  for (i in 1:length(pureList)) {
    purei <- sort(pureList[[i]])   ### For simulation purpose only.
    if (length(purei) != 1) {
      firstPure <- purei[1]
      pos <- firstPure
      neg <- c()
      for (j in 2:length(purei)) {
        purej <- purei[j]
        if (Sigma[firstPure, purej] < 0)
          neg <- c(neg, purej)
        else
          pos <- c(pos, purej)
      }
      if (length(neg) == 0)
        neg <- list()
      signPureList[[i]] <- list(pos = pos, neg = neg)
    } else
      signPureList[[i]] <- list(pos = purei, neg = list())
  }
  return(signPureList)
}
