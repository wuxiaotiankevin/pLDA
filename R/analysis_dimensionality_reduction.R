# input
# dat: preprocessed count matrix, gene by sample
# k: number of topics
# output
# sample by n_topics matrix, gamma
# saved lda fit object at ./output
#' @export
runLDA <- function(dat, k, fpath=NA) {
  if (is.na(fpath)) {
    fpath <- paste0('../output/lda_ntopics_', k, '.rdata')
  }
  # dat is gene by sample
  if (file.exists(fpath)) {
    load(fpath)
  } else {
    fit <- LDA(t(dat), k, method='VEM')
    save(fit, file = fpath)
  }
  fit@gamma # sample by n_topic
}

# input
# dat: preprocessed count matrix, gene by sample
# k: number of topics
# lambda: penalization parameter
# output
# sample by n_topics matrix, gamma
# saved plda fit object at ./output
# if (.Platform$OS.type=="windows") {
#   Sys.setenv(PATH="%PATH%") # solve couldn't find RBuildTools problem
#   source('D:/Dropbox (Brown)/WuLab/WuWuWu/code/functions_plda/plda_functions.R')
# }
# if (.Platform$OS.type=="unix") {
#   source('/home/xwu3/research/genetics/code/functions_plda/plda_functions.R')
# }
#' @export
runPLDA <- function(dat, k, lambda, fpath=NA) {
  if (is.na(fpath)) {
    fpath <- paste0('../output/plda_ntopics_', k, '_lambda_', lambda, '.rdata')
  }
  if (file.exists(fpath)) {
    load(fpath)
  } else {
    fit <- plda(t(dat), k, lambda)
    save(fit, file = fpath)
  }
  fit$gamma # sample by n_topic
}


# logit
# input matrix with entries in (0, 1)
# output logit of matrix
#' @export
logit <- function(x) {
  log(x/(1-x))
}

# input
# dat: input data for ICA, gene by sample
# k: number of components to be extracted
# output
# ICA result
# save ICA object at ./output
#' @export
runICA <- function(dat, k, ica_preprocess, fpath=NA) {
  if (is.na(fpath)) {
    fpath <- paste0('../output/ica_k_', k, '_', ica_preprocess, '.rdata')
  }
  if (file.exists(fpath)) {
    load(fpath)
  } else {
    fit <- fastICA(t(dat), k) # input data should be obs vs. variables
    save(fit, file = fpath)
  }
  fit$S
}

# input
# dat: input data for PCA, gene by sample
# k: number of PCs to extract
# output
# PCA result
# save PCA object at ./output
#' @export
runPCA <- function(dat, k, pca_preprocess, fpath=NA) {
  if (is.na(fpath)) {
    fpath <- paste0('../output/pca_', pca_preprocess, '.rdata')
  }
  if (file.exists(fpath)) {
    load(fpath)
  } else {
    fit <- prcomp(t(dat), center = TRUE, scale. = TRUE)
    save(fit, file = fpath)
  }
  fit$x[, 1:k]
}
