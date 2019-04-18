# clustering functions for plda, sc3, simlr

# plda
# input 
# fit: plda object
# nclust: number of clusters
# output
# save predicted group to ../output/pred_plda.rdata
#' @export
plda_kmeans <- function(fit, nclust, seed) {
  set.seed(seed)
  logit_gamma <- logit(fit$gamma)
  fit.kmeans <- kmeans(logit_gamma, nclust)
  y_plda <- factor(fit.kmeans$cluster)
  save(y_plda, file = '../output/pred_plda.rdata')
  return(y_plda)
}

# sc3
# input
# dat: count matrix
# y: true celltype labels
# output
# save predicted group to ../output/pred_sc3.rdata
#' @export
runSC3 <- function(dat, y, seed, nclust=NA, njobs=20, fpath=NA) {
  if (is.na(nclust)) {
    nclust <- length(unique(y))
  }
  if (is.na(fpath)) {
    fpath <- paste0('../output/pred_sc3_nclust_', nclust, '.rdata')
  }
  
  if (file.exists(fpath)) {
    load(fpath)
    return(y.sc3)
  } else {
    set.seed(seed)
    # create a SingleCellExperiment object
    sce <- SingleCellExperiment(
      assays = list(
        counts = as.matrix(dat),
        logcounts = log2(as.matrix(dat) + 1)
      ), 
      colData = y
    )
    
    # define feature names in feature_symbol column
    rowData(sce)$feature_symbol <- rownames(sce)
    # remove features with duplicated names
    sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
    
    # define spike-ins
    isSpike(sce, "ERCC") <- grepl("ERCC", rowData(sce)$feature_symbol)
    
    # run sc3
    sce <- sc3(sce, ks = nclust, biology = TRUE, n_cores=njobs)
    
    y.sc3 <- colData(sce)[, paste0('sc3_', nclust, '_clusters')]
    save(y.sc3, file = fpath)
    return(y.sc3)
  }
}

# SLIMR
#' @export
runSIMLR <- function(dat, y, seed, nclust=NA, fpath=NA) {
  if (is.na(nclust)) {
    nclust <- length(unique(y))
  }
  
  if (is.na(fpath)) {
    fpath <- paste0('../output/pred_simlr_nclust_', nclust, '.rdata')
  }
  
  if (file.exists(fpath)) {
    load(fpath)
    return(y.simlr)
  } else {
    set.seed(seed)
    fit.simlr <- SIMLR(X = log10(dat+1), c = nclust, cores.ratio = 0.1)
    y.simlr <- factor(fit.simlr$y$cluster)
    save(y.simlr, file = fpath)
    return(y.simlr)
  }
}
