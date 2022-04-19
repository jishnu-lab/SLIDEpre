#' Generate Node Network
#'
#' Create a pdf of the node network by cluster. EXPAND THIS.
#'
#' @param x a data matrix of dimensions \eqn{n \times p}
#' @param sigma a sample correlation matrix of dimensions \eqn{p \times p}
#' @param er_res the results from running Essential Regression with either
#' \link[priorER]{priorER()} or \link[plainER]{plainer()}
#' @param filename the output file name
#' @param equal_var a boolean flag indicating whether the columns of \code{x} have equal variance
#' @param merge a boolean flag indicating merge type

siNetwork <- function(x, sigma, er_res, filename, equal_var = F, merge = F) {
  n <- nrow(x);  p <- ncol(x) #### feature matrix dimensions
  if (equal_var) {
    se_est <- rep(1, p)
  } else {
    se_est <- apply(x, 2, stats::sd) #### get sd of columns for feature matrix
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

  grDevices::pdf(file = filename)
  for (i in 1:length(clust_unlist)) {
    graphics::plot(network,
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
    graphics::legend("bottomleft",
                     legend = c("pure", "mixed", "absent"),
                     pt.bg  = c("salmon", "skyblue", "grey"),
                     pch    = c(21, 21, 21, 23, 23),
                     cex    = 0.5,
                     bty    = "n",
                     title  = "Node/Edge Legend")
  }
  grDevices::dev.off()
  return("Finished!")
}
