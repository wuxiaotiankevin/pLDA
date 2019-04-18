# # debug
# ddir <- '/home/xwu3/research/genetics/code/blood/clustering_comparison'
# load(paste0(ddir, '/output/plda_ntopics_10_lambda_100.rdata')) 
# library(gplots)
# plotGamma(fit)

# heatmap for gamma matrix

#' Rearrange columns of two matrices for match in a greedy way
#' 
#' Rearrange columns of two matrices for match in a greedy way
#' @export
#' @param mat1 Matrix to align to.
#' @param mat2 Matrix to be rearranged.
#' @return Vector of order to be applied to columns of mat2.
#' @examples 
#' \dontrun{
#' column_order <- greedyRearrangeColumn(mat1, mat2)
#' mat2 <- mat2[, column_order]
#' }
greedyRearrangeColumn <- function(mat1, mat2) {
  ncol <- ncol(mat1)
  not_picked_col <- 1:ncol
  res <- numeric()
  for (i in 1:ncol) {
    # calculate distance of current columns in mat1 with all not picked columns in mat2
    dist_col <- numeric()
    for (j in not_picked_col) {
      dist_col <- c(dist_col, sqrt(sum((mat1[, i]-mat2[, j])^2)))
    }
    idx_match <- which.min(dist_col)
    res <- c(res, not_picked_col[idx_match])
    not_picked_col <- not_picked_col[-idx_match]
  }
  return(res)
}



