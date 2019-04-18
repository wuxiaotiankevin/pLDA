# input: trained beta, new data matrix with matched genes
# output: estimated beta

#' Estimate new gamma matrix
#'
#' Calculate topic by cell (gamma) matrix for new expression data
#' with estimated topic by gene matrix.
#'
#' @export
#' @param x Input matrix. A matrix of positive integers where rows represents cells (documents) and column represents genes (words).
#' @param log_beta Estimated log beta matrix from previous run of pLDA. Number of topics by vocabulary size.
#' @return Gamma matrix for the new expression data.
#' @examples
#' \dontrun{
#' plda_estimate_gamma_given_beta(x=new_expression_data, log_beta=plda_model$logProbW)
#' }
plda_estimate_gamma_given_beta <- function(x, log_beta, seed=2018, verbose=FALSE) {
  set.seed(seed)
  n_vocab <- ncol(x)
  n_docs <- nrow(x)
  genes <- colnames(x)
  k <- nrow(log_beta)

  # parse input
  docs.parse <- function(docs.line) {
    words <- 1:n_vocab-1
    counts <- docs.line
    idx <- which(counts!=0)
    words <- words[idx]
    counts <- counts[idx]
    docs.out <- list(words=words, counts=counts, dlength=length(words), total=sum(counts))
    return(docs.out)
  }

  docs <- apply(x, 1, docs.parse)
  corpus <- list(docs=docs, nterms=n_vocab, ndocs=n_docs)


  # other parameters
  estAlpha = TRUE
  MAX.ES.IT = 50
  ES.CONV.LIM = 1e-6
  EM.CONV = 1e-2 # 1e-4
  MAX.EM = 50
  alpha = 50/k;

  # init the model randomly
  cwinit = matrix(runif(k*corpus$nterms),k,corpus$nterms) + 1/corpus$nterms
  ldamodel = list(logProbW=matrix(rep(0,k*corpus$nterms),k,corpus$nterms), alpha = 1, ntopics=k, nterms=corpus$nterms)
  sstats = list(ndocs=0,classword=cwinit,k,corpus$nterms,classtotal=rowSums(cwinit), alpha.ss = 0)
  ldamodel = mstep.beta.start(ldamodel,sstats, 0)
  ldamodel$alpha = alpha

  # Initial beta to be the beta provided
  ldamodel$logProbW <- log_beta

  like.hist = c() # keep track of the likelihood
  likelihood.prev = 0
  numit = 0
  hasConverged = FALSE
  nTopics = ldamodel$ntopics
  nTerms = ldamodel$nterms
  phi = matrix(rep(0,nTopics*nTerms),nTopics, nTerms)

  # # # # Run variational expectation-maximization # # # #
  # while (!hasConverged){
    numit = numit + 1
    if (verbose) {
      print(sprintf("----- EM Iteration %i ----- ", numit))
    }

    # reset sufficient statistics and likelihood
    sstats$classword = matrix(rep(0,nTopics*nTerms), nrow=nTopics, nTerms)
    sstats$classtotal = rep(0,nTopics)
    sstats$ndocs = 0
    sstats$alpha.ss = 0
    likelihood = 0

    # # # do E-step # # #
    gammaOut <- matrix(NA, nrow = corpus$ndocs, ncol = nTopics)
    for (d in 1:corpus$ndocs){
      # # do posterior inference # #

      # initialize the document specific variables
      doc.oldlike = 0
      doc.length = corpus$docs[[d]]$dlength
      doc.totlen = corpus$docs[[d]]$total
      gammav = rep(ldamodel$alpha + doc.totlen/nTopics, nTopics)
      digamma.gam = rep(digamma(ldamodel$alpha + doc.totlen/nTopics), nTopics)
      phi = matrix(rep(1/nTopics, doc.length*nTopics), nrow=doc.length, ncol=nTopics)
      oldphi = phi[1,]

      # compute posterior dirichlet
      estep.converged = FALSE
      numits.es = 0;
      while (!estep.converged){
        numits.es = numits.es + 1
        # cpp implementation
        do_e_step_C(doc.length, oldphi, nTopics, phi,
                    digamma.gam, ldamodel$logProbW,
                    corpus$docs[[d]]$words,
                    gammav,
                    corpus$docs[[d]]$counts)

        # determine if the documents likelihood has converged
        doc.like <- compute_likelihood_C(ldamodel$ntopics,
                                         ldamodel$alpha,
                                         corpus$docs[[d]]$dlength,
                                         ldamodel$logProbW,
                                         corpus$docs[[d]]$words,
                                         corpus$docs[[d]]$counts,
                                         phi,
                                         gammav)
        convfrac = (doc.oldlike - doc.like) / doc.oldlike
        doc.oldlike = doc.like


        if (convfrac < ES.CONV.LIM || numits.es > MAX.ES.IT){
          estep.converged = TRUE
          # print(sprintf("leaving E-step after %i iterations and convfrac: %1.3e, doc-likelihood: %1.3e", numits.es, convfrac, doc.like))
          # plot(doc.histlike)
        }
      } # end while e-step has not converged

      # # update the sufficient statistics for the M-step # #
      gamma.sum = sum(gammav)
      sstats$alpha.ss = sstats$alpha.ss + sum(sapply(gammav,digamma))
      sstats$alpha.ss = sstats$alpha.ss - nTopics*digamma(gamma.sum)

      do_m_step_C(doc.length, nTopics,
                  corpus$docs[[d]]$counts,
                  phi,
                  corpus$docs[[d]]$words,
                  sstats$classword,
                  sstats$classtotal)

      sstats$ndocs = sstats$ndocs + 1
      likelihood = likelihood + doc.like
      gammaOut[d, ] <- gammav
    } # end for each document

  if (verbose) print("----- Finished -----")

  # calculate topic document matrix, gamma
  gammaOut/rowSums(gammaOut)
}
