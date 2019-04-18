# # debug
# ddir <- '/home/xwu3/research/genetics/code/blood/clustering_comparison'
# load(paste0(ddir, '/output/plda_ntopics_10_lambda_100.rdata')) 
# library(gplots)
# plotGamma(fit)

# heatmap for gamma matrix

#' Plot heatmap for gamma matrix
#' 
#' Plot heatmap for gamma matrix
#' @export
#' @param fit Model fit from plda()
#' @return Heatmap of gamma matrix.
#' @examples 
#' \dontrun{
#' plotGamma(fit=plda_model)
#' }
plotGamma <- function(fit) {
  mat.gamma <- fit$gamma
  colnames(mat.gamma) <- paste('topic_', 1:ncol(mat.gamma))
  heatmap.2(t(mat.gamma), density.info="none", 
            trace="none", 
            Colv=TRUE, 
            labCol = FALSE,
            Rowv=TRUE, 
            dendrogram="both", 
            main = 'Gamma Matrix')
}



