# input
# dat: count matrix, gene by sample
# lam_values: lambda values
# n_rep: repeated measure at each lambda value
# idx.hk: index of housekeeping genes

# output
# plot number of interesting genes for different lam values
# return table of lam_values and # interesting genes
get_n_interesting_genes <- function(dat, k, lam, idx_hk, cutoff, warm.start.log.beta, seed=2018, fpath=NA) {
  # if (is.na(fpath)) {
  #   fpath <- paste0('../output/plda_ntopics_', k, '_lambda_', lam, '.rdata')
  # }
  fpath <- file.path(fpath, paste0('plda_ntopics_', k, '_lambda_', lam, '.rdata'))
  set.seed(seed)
  # run plda and save
  if (file.exists(fpath)) {
    load(fpath)
  } else {
    fit <- plda(t(dat), k, lam, warm.start.log.beta=warm.start.log.beta)
    save(fit, file = fpath)
  }
  # get n interesting genes using cutoff
  bmat <- exp(fit$logProbW) # k by V
  nbmat <- sweep(bmat, 2, apply(bmat, 2, mean), FUN = "/") # normalize by dividing the mean expr
  vnb <- apply(nbmat, 2, var) # calculate variance
  n_interesting_fixed_cutoff <- sum(vnb>cutoff)
  idx_interesting_fixed_cutoff <- which(vnb>cutoff)
  # get n interesting genes using hk genes in each round
  cutoff_current <- mean(vnb[idx_hk])
  n_interesting_varied_cutoff <- sum(vnb>cutoff_current)
  idx_interesting_varied_cutoff <- which(vnb>cutoff_current)
  
  res.tmp <- data.frame(lam, 
                        n_interesting_fixed_cutoff, 
                        n_interesting_varied_cutoff)
  res.tmp$idx_interesting_fixed_cutoff <- list(idx_interesting_fixed_cutoff)
  res.tmp$idx_interesting_varied_cutoff <- list(idx_interesting_varied_cutoff)
  
  return(res.tmp)
}

#' pLDA Choose Shrinkage Parameter
#' 
#' Run pLDA under different lambda values. Plot and return output. Interesting genes are defined as those with variance in beta matrix greater than average variance of housekeeping genes.
#' @export
#' @param dat Input matrix. A matrix of positive integers where rows represent genes (words) and columns represent cells (documents).
#' @param k Number of topics.
#' @param lam_values A vector of shrinkage parameter values.
#' @param idx_hk Index (column number) of housekeeping gene.
#' @param njobs Number of cores used for the calculation. Default 1. 
#' @param n_rep Number of repetition for each lam_values. Default 1.
#' @param fdir Output file directory. Default ../output/
#' @param seed Set random seed. Default 2018.
#' @return Plot number of interesting genes for different lambda values. Return table of lam_values and corresponding number of interesting genes.
#' @examples
#' \dontrun{
#' plda_choose_lambda(dat=cell_by_gene_expr_matrix, k=10, lam_values=c(10, 10^2, 10^3),  idx_hk = 1:10)
#' }
plda_choose_lambda <- function(dat, k, lam_values, idx_hk, njobs=1, n_rep=1, fdir=NA, seed=2018) {
  require(plyr)
  set.seed(seed)
  if (is.na(fdir)) {
    fdir <- '../output/'
  }
  print(paste0('Calculating for ', k, ' topics.'))
  # create output fulder
  dir.create(fdir, showWarnings = FALSE)
  # run LDA and save
  if (file.exists(file.path(fdir, paste0('lda_ntopics_', k, '.rdata')))) {
    load(file.path(fdir, paste0('lda_ntopics_', k, '.rdata')))
  } else {
    fit <- LDA(t(dat), k, control = list(seed=seed))
    save(fit, file = file.path(fdir, paste0('lda_ntopics_', k, '.rdata')))
  }
  # warm start log beta
  warm.start.log.beta <- fit@beta
  
  # calculate cutoff
  bmat <- exp(fit@beta) # k by V
  nbmat <- sweep(bmat, 2, apply(bmat, 2, mean), FUN = "/") # normalize by dividing the mean expr
  vnb <- apply(nbmat, 2, var) # calculate variance
  cutoff <- mean(vnb[idx_hk])
  
  # run plda, get number of interesting genes
  res <- data.frame()
  if (njobs==1) {
    for (lam in lam_values) {
      set.seed(seed)
      print(paste0('Processing lambda=',lam, '...'))
      tmp <- get_n_interesting_genes(dat, k, lam, idx_hk, cutoff, warm.start.log.beta=warm.start.log.beta, fpath = fdir)
      res <- rbind(res, rbind.fill(tmp))
    }
  } else if (njobs>1) {
    set.seed(seed)
    print(paste0('Running ', njobs, ' jobs in parallel.'))
    tmp <- mclapply(1:length(lam_values), 
                    function(x) get_n_interesting_genes(dat, k, lam_values[x], idx_hk, cutoff, warm.start.log.beta=warm.start.log.beta, fpath = fdir), 
                    mc.cores = njobs)
    res <- rbind(res, do.call(rbind, tmp))
  }
  
  # save result
  save(res, file = file.path(fdir, paste0('choose_lambda_ntopics_', k, '_maxlam_', max(lam_values), '.rdata')))
  
  # make plot
  # fixed cutoff
  pdf(file.path(fdir, paste0('choose_lambda_ntopics_', k, '_maxlam_', max(lam_values), '_fixed_cutoff.pdf')))
  plot(res$lam, res$n_interesting_fixed_cutoff, 
       xlab = "lambda", ylab = "Number of genes selected",
       main = "Number of Genes Selected for Different lambda, fixed cutoff"
  )
  dev.off()
  # varied cutoff
  pdf(file.path(fdir, paste0('choose_lambda_ntopics_', k, '_maxlam_', max(lam_values), '_varied_cutoff.pdf')))
  plot(res$lam, res$n_interesting_varied_cutoff, 
       xlab = "lambda", ylab = "Number of genes selected",
       main = "Number of Genes Selected for Different lambda, varied cutoff"
  )
  dev.off()
  
  res
}
