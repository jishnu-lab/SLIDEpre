#' Estimate list of pure node indices for given \eqn{\Sigma} and \eqn{\delta}.
#' This code is an implementation of Algorithm 1 from Bing et al. (2020).
#'
#' @param off_Sigma a sample correlation matrix of dimensions \eqn{p \times p}
#' @param delta \eqn{\delta}, a numerical constant
#' @param Ms the largest absolute values of each row of \code{off_Sigma}
#' @param arg_Ms vector of column indices at which the values in \code{Ms} are achieved in \code{off_Sigma}
#' @param se_est standard deviations of features (columns of \code{X})
#' @param merge boolean indicating merge style
#' @return a list including the list of estimated pure node indices and a vector
#' of the estimated pure node indices

siNetwork <- function(x, sigma, er_res, filename, equal_var = F, merge = F) {
  n <- nrow(x);  p <- ncol(x) #### feature matrix dimensions
  if (equal_var) {
    se_est <- rep(1, p)
  } else {
    se_est <- apply(x, 2, sd) #### get sd of columns for feature matrix
  }

  delta <- er_res$optDelta * sqrt(log(max(p, n)) / n)

  #### get node list
  nodes <- colnames(sigma)
  if (length(nodes) == 0) {
    nodes <- paste0("feat", seq(1, p))
    colnames(sigma) <- nodes
    rownames(sigma) <- nodes
  }

  #### get feature lists from er_res
  feats <- readER(er_res)
  clusters <- feats$clusters
  clust_unlist <- list()
  for (i in 1:length(clusters)) {
    clust <- clusters[[i]]
    clust_nodes <- unlist(clust)
    clust_nodes <- indName(clust_nodes, colnames(sigma), to_ind = F)
    clust_unlist[[length(clust_unlist) + 1]] <- clust_nodes
  }
  pures <- feats$pure_vars
  mixeds <- feats$mix_vars

  #### get absolute value of covariance matrix
  off_Sigma <- abs(sigma)

  #### set entries on main diagonal to 0
  diag(off_Sigma) <- 0

  #### calculate the maximal absolute value for each row of the given matrix
  result_Ms <- FindRowMax(off_Sigma)
  Ms <- result_Ms$M #### maximal abs values
  arg_Ms <- result_Ms$arg_M #### first index where max abs values are achieved

  #### get edges list
  edges <- c()
  for (i in 1:nrow(off_Sigma)) { #### loop through rows
    from_node <- i
    from_node_name <- nodes[i]
    row_i <- off_Sigma[i,]
    Si <- FindRowMaxInd(i, Ms[i], arg_Ms[i], row_i, delta, se_est)
    if (length(Si) > 0) {
      for (j in 1:length(Si)) {
        to_node <- Si[j]
        to_node_name <- nodes[to_node]
        weight <- round(sigma[from_node, to_node], 2)
        edges <- rbind(edges, c(from_node_name, to_node_name, weight))
      }
    }
  }
  colnames(edges) <- c("from", "to", "weight")
  edges <- as.data.frame(edges)

  #### get pure nodes
  node_nums <- seq(1, p)
  is_pure <- ifelse(node_nums %in% pures, "pure", ifelse(node_nums %in% mixeds, "mixed", "absent"))
  nodes <- data.frame("node" = nodes, "type" = is_pure)

  #### make network graph
  network <- igraph::graph.data.frame(edges, nodes, directed = T)
  V(network)$color <- ifelse(V(network)$type == "pure", "salmon", ifelse(V(network)$type == "mixed", "skyblue", "grey"))
  E(network)$color <- "black"
  pdf(file = filename)
  for (i in 1:length(clust_unlist)) {
    plot(network,
         edge.arrow.size = 0.2,
         vertex.label.cex = 0.05,
         vertex.label.font = 2,
         vertex.label.color = "black",
         vertex.size = p * 0.2,
         layout = igraph::layout.circle(network),
         mark.groups = clust_unlist[[i]],
         mark.border = NA,
         main = paste0("Cluster ", i))
  }
  dev.off()
  return("Finished!")
}
