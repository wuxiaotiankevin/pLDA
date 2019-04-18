#' Human housekeeping genes.
#'
#' A dataset containing 4061 human housekeeping genes.
#'
#' @format A data frame with 4061 rows and 2 variables:
#' \describe{
#'   \item{ensembl_gene_id}{Ensemble gene id.}
#'   \item{hgnc_symbol}{HGNC sympol.}
#' }
#' @source \url{https://doi.org/10.1016/j.tig.2013.05.010}
#' \url{https://www.tau.ac.il/~elieis/HKG/}
"hkGenes"

# load(paste0("/home/xwu3/research/genetics/data/housekeeping_genes/HK_genes.Rdata"))
# hkGenes <- g_list
# load(paste0("/home/xwu3/research/genetics/data/housekeeping_genes/HK_genes_shortlist.Rdata"))
# hkGenesShortlist <- g_list
# save(hkGenes, file = './data/hkGenes.rdata')
# save(hkGenesShortlist, file = './data/hkGenesShortlist.rdata')


#' Shortlist of human housekeeping genes.
#'
#' A dataset containing a shortlist of 12 human housekeeping genes.
#'
#' @format A data frame with 12 rows and 2 variables:
#' \describe{
#'   \item{ensembl_gene_id}{Ensemble gene id.}
#'   \item{hgnc_symbol}{HGNC sympol.}
#' }
#' @source \url{https://doi.org/10.1016/j.tig.2013.05.010}
#' \url{https://www.tau.ac.il/~elieis/HKG/}
"hkGenesShortlist"
