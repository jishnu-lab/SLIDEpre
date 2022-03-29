#' Merge a group with other groups containing common nodes.
#'
#' @param groupList a list of groups of node indices
#' @param groupVec a vector of node indices
#' @return a list of the merged results

Merge <- function(groupList, groupVec) {
  # merge the new group with the previous ones which have common nodes
  if (length(groupList) != 0) {
    for (i in 1:length(groupList)) {
      common_nodes <- intersect(groupList[[i]], groupVec)
      if (length(common_nodes) != 0) {
        groupList[[i]] <- common_nodes
        return(groupList)
      }
    }
  }
  groupList <- append(groupList, list(groupVec))
  return(groupList)
}
