#' Calculate the inconsistency of a \code{hclust} object.
#'
#' @description Calculates the inconsistency of a 'hclust' object and splits the data into clusters using a inconsistency-based cut-off.
#' @param hclust_obj the input data - see details.
#' @param cut_point a value of inconsistency at which to split the clusters into groups.
#'
#' @details This function calcualtes the inconsistency of a  \code{hclust} object.
#' This provides the same functionality as the \code{inconsistent} function in Matlab and the 'scipy' library in Python.
#' The inconsistency can be calculated as \deqn{\frac{h_{i} - \bar{h}}{s}} where \eqn{h_{i}} is the height of the given subtree;
#' and \eqn{\bar{h}} and \eqn{s} are respectively the mean and standard deviation of the height of the subtree and any non-leaf subtrees.
#' Note that the mean and standard deviationare calculated over 1, 2, or 3 heights depending on whether each subtree has 0, 1, or 2 non-leaf subtrees (with the resulting degenerecies).
#' By definition, a subtree where both of its subtrees are leaves has an inconsistency of 0.
#'
#' @return The original object of class \code{hclust}, with additional components containing the calculated inconsistency at each merge point and,
#' if \code{cut_point} is specified, the resulting clusters.
#' \item{inconsistency}{ The calculated inconsistency at each merge point.}
#' \item{clusters}{ A data frame with 2 variables: the label and the assigned cluster. (Only present if a \code{cut_point} is specified.)}
#'
#' @seealso \code{\link{hclust}}
#' @export
#'
#' @examples
#' hc <- hclust(dist(USArrests), "ave")
#' hc_obj <- inconsistent(hc, 1)
#'
#' # the first 6 merges all involve leaves, so the inconsistency is 0.
#' head(with(hc_obj, cbind(merge, height, inconsistency)))
#'
#' # the last 6 merges are all non-leaves, so the inconsistency can be calculated
#' tail(with(hc_obj, cbind(merge, height, inconsistency)))
#'
#' # Each element of the tree is assigned to a cluster.
#' head(hc_obj$clusters)
#' table(hc_obj$clusters$cluster)


# function to mimic the "inconsistent" function in matlab
inconsistent <- function(hclust_obj,cut_point = NULL) {
  if(class(hclust_obj) != "hclust") stop("'hclust_obj' must be of type 'hclust'.")

  n <- nrow(hclust_obj$merge)

  hclust_data <- with(hclust_obj, data.frame(
    obj_left = merge[, 1],
    obj_right = merge[, 2],
    height = height
  ))

  # calcualte the inconsistency of the hclust_obj object
  hclust_data <- within(hclust_data,{
    cluster_left <- obj_left >  0
    cluster_right <- obj_right > 0

    height_left <- ifelse(cluster_left, height[abs(obj_left)], NA)
    height_right <- ifelse(cluster_right, height[abs(obj_right)], NA)
    height_mat <- cbind(cbind(height,height_left,height_right))
    n_heights <- rowSums(!is.na(height_mat))

    height_mean <- rowMeans(height_mat,na.rm = TRUE)
    height_sd <- sqrt(rowVar(height_mat,na.rm = TRUE))


    inconsistency <- ifelse(n_heights == 1 | height_sd < 1e-8 ,0, (height - height_mean) / height_sd)
    rm(height_mat)
  })
  hclust_obj$inconsistency <- hclust_data$inconsistency

  # now, if a cut_point is set, group the clusters up until the inconsistency threshold is reached.
  if (!is.null(cut_point)) {
    groups <- within(data.frame(label = hclust_obj$labels), {
      cluster <- 1:length(label)
    })
    splits <- rep(FALSE,n)
    for (i in 1:n) {
      if(hclust_data$inconsistency[i] > cut_point |  # check inconsistency of current point AND each subtree
          (hclust_data$obj_left[i] > 0 && splits[hclust_data$obj_left[i]]) |
          (hclust_data$obj_right[i] > 0 && splits[hclust_data$obj_right[i]])) {
        splits[i] <- TRUE
      } else {
        members <- PQPQ:::find_nodes(hclust_data$obj_left,hclust_data$obj_right, i)
        lowest_group <- min(groups$cluster[members])
        groups$cluster[members] <- min(groups$cluster[members])
      }
    }

    # clean up clusters
    groups$cluster <- as.numeric(factor(groups$cluster))

    hclust_obj$clusters <- groups

  }

  hclust_obj
}
