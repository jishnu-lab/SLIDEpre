#' Merge the provided vector of nodes into the list of vectors of nodes but use union rather than the intersection.
#'
#' @param groupList a list of groups of node indices
#' @param groupVec a vector of node indices
#' @return a list of the merged results

Merge_union <- function(groupList, groupVec) {
  # merge the new group with the previous ones which have common nodes
  if (length(groupList) != 0) {
    common_groups <- sapply(groupList, FUN = function(x, y) {
      length(intersect(x, y))
    }, y = groupVec)
    common_inds <- which(common_groups > 0)
    if (length(common_inds) > 0){
      new_group <- unlist(lapply(common_inds,
                                 FUN = function(x, y){y[[x]]}, y = groupList))
      remain_group <- lapply(which(common_groups == 0),
                             FUN = function(x, y){y[[x]]}, y = groupList)
      groupList <- append(remain_group, list(union(groupVec, new_group)))
      return(groupList)
    }
  }
  groupList <- append(groupList, list(groupVec))
  return(groupList)
}
