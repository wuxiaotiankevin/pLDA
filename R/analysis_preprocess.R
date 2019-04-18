# preprocess
# keep when vars/means>1 and row max > 10
#' Preprocess pLDA input data
#' 
#' Preprocess pLDA input data. Keep only genes with across samples variance/mean > var.mean.cutoff and across samples max > row.max.cutoff.
#' @export
#' @param data Input matrix. A matrix of positive integers where rows represents cells (documents) and column represents genes (words).
#' @param var.mean.cutoff Row variance/mean cutoff.
#' @param row.max.cutoff Row max cutoff.
#' @return Preprocessed data.
#' @examples
#' \dontrun{
#' preprocess(data=cell_by_gene_expr_matrix)
#' }
preprocess <- function(data, var.mean.cutoff=1, row.max.cutoff=10) {
  data <- data[which(apply(data, 1, sum)>0), ]
  vars <- apply(data, 1, var)
  means <- apply(data, 1, mean)
  keep1 <- which(vars/means > var.mean.cutoff)
  keep2 <- which(apply(data, 1, max)>row.max.cutoff)
  idx <- intersect(keep1, keep2)
  data[idx, ]
}

# quantile normalization
# norm <- normalizeQuantiles(dat)

# cpm
cpm <- function(dat, take.log=TRUE) {
  if (take.log) {
    log2(t(t(dat)/colSums(dat))*1e6 + 1)
  } else {
    t(t(dat)/colSums(dat))*1e6
  }
}