#' Find the frequency with which each mixed variable appears in the results of Essential Regression and
#' return a table of the frequencies as well as the clusters in which each variable appears.
#'
#' @param er_res an object return by \code{\link{ER()}}
#' @return a matrix with the same number of rows as there are mixed variables in \code{er_res} and with
#' a column containing frequencies as well as a column for each cluster that contain a 1 if that variable
#' is found in that cluster and a 0 if not.
#' @export

getFreq <- function(er_res) {
  results <- readER(er_res)
  all_feats <- unlist(results$clusters)
  mix_vars <- results$mix_var
  freq_tab <- data.frame("var" = mix_vars)

  clusters <- results$clusters
  for (i in 1:length(clusters)) {
    clust <- clusters[[i]] %>% unlist()
    cl_col <- ifelse(mix_vars %in% clust, 1, 0) %>% as.data.frame()
    colnames(cl_col) <- paste0("cluster", i)
    freq_tab <- cbind(freq_tab, cl_col)
  }

  row_sums <- apply(freq_tab[, 2:ncol(freq_tab)], 1, sum)
  freq_tab$freq <- row_sums
  freq_tab <- freq_tab %>% as.data.frame() %>%
    dplyr::select(var, freq, everything())
  freq_tab <- freq_tab[order(freq_tab$freq), ]
  return (freq_tab)
}
