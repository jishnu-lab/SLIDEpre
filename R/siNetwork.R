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

  delta <- er_res$opt_delta * sqrt(log(max(p, n)) / n)

  #### get node list
  nodes <- colnames(sigma)
  if (length(nodes) == 0) {
    nodes <- paste0("feat", seq(1, p))
    colnames(sigma) <- nodes
    rownames(sigma) <- nodes
  }
  node_names <- data.frame("name" = nodes, "position" = seq(1, p))

  #### get feature lists from er_res
  feats <- readER(er_res)
  clusters <- feats$clusters
  clust_unlist <- list()
  for (i in 1:length(clusters)) {
    clust <- clusters[[i]]
    clust_nodes <- unlist(clust)
    clust_nodes <- node_names[node_names$position %in% clust_nodes, 1]
    clust_unlist[[length(clust_unlist) + 1]] <- clust_nodes
  }
  pures <- feats$pure_vars
  pures <- node_names[node_names$position %in% pures, 1]
  mixeds <- feats$mix_vars
  mixeds <- node_names[node_names$position %in% mixeds, 1]

  #### get absolute value of covariance matrix
  abs_sigma <- abs(sigma)

  #### set entries on main diagonal to 0
  diag(abs_sigma) <- 0

  #### calculate the maximal absolute value for each row of the given matrix
  maxes <- findRowMax(abs_sigma)
  max_vals <- maxes$max_vals #### maximal abs values
  max_inds <- maxes$max_inds #### first index where max abs values are achieved

  #### get edges list
  edges <- c()
  for (i in 1:nrow(abs_sigma)) { #### loop through rows
    from_node <- i
    from_node_name <- node_names[node_names$position == i, 1]
    row_i <- abs_sigma[from_node_name, ]
    si <- findRowMaxInd(i = i,
                        max_val = max_vals[i],
                        max_ind = max_inds[i],
                        row_i = row_i,
                        delta = delta,
                        se_est = se_est)
    if (length(si) > 0) {
      for (j in 1:length(si)) {
        to_node <- si[j]
        to_node_name <- node_names[node_names$position == to_node, 1]
        weight <- round(sigma[from_node_name, to_node_name], 2)
        color <- ifelse(abs(sigma[from_node_name, to_node_name]) < max_vals[i], "#893FC9", "#3AAA30")
        edges <- rbind(edges, c(from_node_name, to_node_name, weight, color))
      }
    }
  }
  colnames(edges) <- c("from", "to", "weight", "color")
  edges <- as.data.frame(edges)

  #### redo node labels
  nodes <- data.frame("node" = node_names$name)
  is_pure <- ifelse(nodes$node %in% pures, "pure", ifelse(nodes$node %in% mixeds, "mixed", "absent"))
  nodes$type <- is_pure

  node_labs <- c()
  alt_sigma <- sigma
  diag(alt_sigma) <- 0
  for (i in 1:length(nodes$node)) {
    node_labs <- c(node_labs, paste0(nodes$node[i], "\n", round(max(abs(alt_sigma[nodes$node[i], ])), 2)))
  }
  nodes$labs <- node_labs

  #### make network graph
  network <- igraph::graph.data.frame(edges, nodes, directed = T)
  igraph::V(network)$color <- ifelse(igraph::V(network)$type == "pure", "salmon",
                                     ifelse(igraph::V(network)$type == "mixed", "skyblue", "grey"))
  igraph::E(network)$color <- "magenta"

  #### position edge labels
  layt <- igraph::layout.circle(network)
  ELx <- rep(0, igraph::ecount(network))
  ELy <- rep(0, igraph::ecount(network))
  for(i in 1:igraph::ecount(network)) {
    ends <- igraph::ends(network, i)
    end1 <- which(nodes$node == ends[1])
    end2 <- which(nodes$node == ends[2])
    from <- layt[end1, ]
    to <- layt[end2, ]
    dir_vec <- to - from
    ELx[i] = from[1] + dir_vec[1] * 0.1
    ELy[i] = from[2] + dir_vec[2] * 0.1
  }

  pdf(file = filename)
  for (i in 1:length(clust_unlist)) {
    plot(network,
         edge.arrow.size = 0.2,
         edge.label.font = 2,
         edge.label.color = "black",
         edge.label.cex = 0.1,
         edge.label = edges$weight,
         edge.label.x = ELx,
         edge.label.y = ELy,
         vertex.label = nodes$labs,
         vertex.label.cex = 0.05,
         vertex.label.font = 2,
         vertex.label.color = "black",
         vertex.size = p * 0.2,
         xlim = range(layt[, 1]),
         ylim = range(layt[, 2]),
         layout = layt,
         mark.groups = clust_unlist[[i]],
         mark.border = NA,
         main = paste0("Cluster ", i))
    legend(
      "bottomleft",
      legend = c("pure", "mixed", "absent"),
      pt.bg  = c("salmon", "skyblue", "grey"),
      pch    = c(21, 21, 21, 23, 23),
      cex    = 0.5,
      bty    = "n",
      title  = "Node/Edge Legend"
    )
  }
  dev.off()
  return("Finished!")
}
