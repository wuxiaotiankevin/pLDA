# Calculate left inverse of matrix X
left.inverse <- function(X) {
  solve(t(X) %*% X) %*% t(X)
}


# input: gamma matrix and corresponding dataset
# output: estimated beta

#' Estimate beta matrix given data and gamma matrix
#'
#' Calculate topic by gene (beta) matrix given data and gamma matrix.
#'
#' @export
#' @param x Input matrix. A matrix of positive integers where rows represents cells (documents) and column represents genes (words). Here the input matrix should be associated with gamma.
#' @param gamma Topic by sample matrix gamma associated with input matrix x.
#' @return log_beta matrix based on x and gamma.
#' @examples
#' \dontrun{
#' plda_estimate_beta(log_beta=plda_model$logProbW, x=new_expression_data)
#' }
plda_estimate_beta_given_gamma <- function(x, theta.true) {
  docs.word.freq <- sweep(x, 1, rowSums(x), FUN = '/')
  beta.from.true.gamma <- left.inverse(theta.true) %*% docs.word.freq
  beta.from.true.gamma
}

# The iteration does not converge
# plda_estimate_beta_given_gamma <- function(x, theta.true, converge_steps, seed=2018, verbose=FALSE) {
#   set.seed(seed)
#   n_vocab <- ncol(x)
#   n_docs <- nrow(x)
#   genes <- colnames(x)
#   k <- ncol(theta.true)
#   lambda <- 0
#   weight.beta <- rep(1, ncol(x))
#   gammaOut <- sweep(theta.true, 1, (apply(x, 1, sum) + k*50/k), FUN = '*')
#   
#   # parse input
#   docs.parse <- function(docs.line) {
#     words <- 1:n_vocab-1
#     counts <- docs.line
#     idx <- which(counts!=0)
#     words <- words[idx]
#     counts <- counts[idx]
#     docs.out <- list(words=words, counts=counts, dlength=length(words), total=sum(counts))
#     return(docs.out)
#   }
#   
#   docs <- apply(x, 1, docs.parse)
#   corpus <- list(docs=docs, nterms=n_vocab, ndocs=n_docs)
#   
#   # other parameters
#   estAlpha = TRUE
#   MAX.ES.IT = 50
#   ES.CONV.LIM = 1e-6
#   EM.CONV = 1e-2 # 1e-4
#   MAX.EM = 50
#   alpha = 50/k
#   
#   # init the model randomly
#   cwinit = matrix(runif(k*corpus$nterms),k,corpus$nterms) + 1/corpus$nterms
#   ldamodel = list(logProbW=matrix(rep(0,k*corpus$nterms),k,corpus$nterms), alpha = 1, ntopics=k, nterms=corpus$nterms)
#   sstats = list(ndocs=0,classword=cwinit,k,corpus$nterms,classtotal=rowSums(cwinit), alpha.ss = 0)
#   ldamodel = mstep.beta.start(ldamodel,sstats, 0)
#   ldamodel$alpha = alpha
#   
#   like.hist = c() # keep track of the likelihood
#   likelihood.prev = 0
#   numit = 0
#   hasConverged = FALSE
#   nTopics = ldamodel$ntopics
#   nTerms = ldamodel$nterms
#   phi = matrix(rep(0,nTopics*nTerms),nTopics, nTerms)
#   
#   # # # # Run variational expectation-maximization # # # #
#   # while (!hasConverged){
#   # converge_steps <- 100
#   for (counter in 1:(converge_steps+5)) {
#     numit = numit + 1
#     if (verbose) {
#       print(sprintf("----- EM Iteration %i ----- ", numit))
#     }
#     
#     # reset sufficient statistics and likelihood
#     sstats$classword = matrix(rep(0,nTopics*nTerms), nrow=nTopics, nTerms)
#     sstats$classtotal = rep(0,nTopics)
#     sstats$ndocs = 0
#     sstats$alpha.ss = 0
#     likelihood = 0
#     likelihood.beta = 0
#     
#     # # # do E-step # # #
#     for (d in 1:corpus$ndocs){
#       if (verbose) print(sprintf("Document: %d, numit: %d",d, numit))
#       # # do posterior inference # #
#       
#       # initialize the document specific variables
#       doc.oldlike = 0
#       doc.length = corpus$docs[[d]]$dlength
#       doc.totlen = corpus$docs[[d]]$total
#       
#       # if (numit==1) {
#       #   gammav = rep(ldamodel$alpha + doc.totlen/nTopics, nTopics)
#       #   digamma.gam = rep(digamma(ldamodel$alpha + doc.totlen/nTopics), nTopics)
#       # } else {
#       # set gamma to input
#       gammav = rep(ldamodel$alpha + doc.totlen/nTopics, nTopics) + 
#         (gammaOut[d, ] - rep(ldamodel$alpha + doc.totlen/nTopics, nTopics)) * min(numit/converge_steps, 1) 
#       digamma.gam = sapply(gammav, digamma)
#       # }
#       
#       print(gammav)
#       print(digamma.gam)
#       print(rep(ldamodel$alpha + doc.totlen/nTopics, nTopics))
#       print(rep(digamma(ldamodel$alpha + doc.totlen/nTopics), nTopics))
#       
#       # initialize phi
#       phi = matrix(rep(1/nTopics, doc.length*nTopics), nrow=doc.length, ncol=nTopics)
#       oldphi = phi[1,]
#       
#       # compute posterior dirichlet
#       estep.converged = FALSE
#       numits.es = 0;
#       while (!estep.converged){
#         numits.es = numits.es + 1
#         # cpp implementation
#         do_e_step_C(doc.length,
#                     oldphi,
#                     nTopics,
#                     phi,
#                     digamma.gam,
#                     ldamodel$logProbW,
#                     corpus$docs[[d]]$words,
#                     gammav,
#                     corpus$docs[[d]]$counts)
#         # if (numit>1) {
#         # reset gamma to input
#         # gammav = rep(ldamodel$alpha + doc.totlen/nTopics, nTopics) + 
#         #   (gammaOut[d, ] - rep(ldamodel$alpha + doc.totlen/nTopics, nTopics)) * min(numit/10, 1) 
#         # digamma.gam = sapply(gammav, digamma)
#         # }
#         
#         # determine if the documents likelihood has converged
#         doc.like <- compute_likelihood_C(ldamodel$ntopics,
#                                          ldamodel$alpha,
#                                          corpus$docs[[d]]$dlength,
#                                          ldamodel$logProbW,
#                                          corpus$docs[[d]]$words,
#                                          corpus$docs[[d]]$counts,
#                                          phi,
#                                          gammav)
#         convfrac = (doc.oldlike - doc.like) / doc.oldlike
#         doc.oldlike = doc.like
#         
#         if (convfrac < ES.CONV.LIM || numits.es > MAX.ES.IT){
#           estep.converged = TRUE
#           # print(sprintf("leaving E-step after %i iterations and convfrac: %1.3e, doc-likelihood: %1.3e", numits.es, convfrac, doc.like))
#           # plot(doc.histlike)
#         }
#       } # end while e-step has not converged
#       
#       # # update the sufficient statistics for the M-step # #
#       gamma.sum = sum(gammav)
#       sstats$alpha.ss = sstats$alpha.ss + sum(sapply(gammav,digamma))
#       sstats$alpha.ss = sstats$alpha.ss - nTopics*digamma(gamma.sum)
#       
#       print(phi)
#       
#       do_m_step_C(doc.length,
#                   nTopics,
#                   corpus$docs[[d]]$counts,
#                   phi,
#                   corpus$docs[[d]]$words,
#                   sstats$classword,
#                   sstats$classtotal)
#       
#       sstats$ndocs = sstats$ndocs + 1
#       likelihood = likelihood + doc.like
#       likelihood.beta = likelihood.beta + compute_beta_involved_likelihood_C(ldamodel$ntopics,
#                                                                              corpus$docs[[d]]$dlength,
#                                                                              ldamodel$logProbW,
#                                                                              corpus$docs[[d]]$words,
#                                                                              corpus$docs[[d]]$counts,
#                                                                              phi)
#       # gammaOut[d, ] <- gammav
#     } # end for each document
#     
#     # # # do M-step # # #
#     # estimate beta
#     ldamodel = mstep.beta.no.shrinkage(ldamodel, sstats, lambda, weight.beta)
#     
#     # add penalty term
#     beta <- exp(ldamodel$logProbW)
#     tmp <- sweep(beta, 2, apply(beta, 2, mean))
#     likelihood.penality <- lambda * sum(apply(tmp^2, 2, sum) * weight.beta)
#     likelihood.docs <- likelihood
#     likelihood <- likelihood.docs - likelihood.penality
#     
#     conv.em = (likelihood.prev - likelihood)/likelihood.prev
#     likelihood.prev = likelihood
#     like.hist[numit] = likelihood
#     
#     # browser()
#     # make sure we're iterating enough for the likelihood to converge'
#     if (conv.em < 0){
#       MAX.ES.IT = MAX.ES.IT*2
#     }
#     # if (((conv.em < EM.CONV && conv.em > 0)  || numit > MAX.EM) && numit > 2){
#     if (((abs(conv.em) < EM.CONV)  || numit > MAX.EM) && numit > 2){
#       # print(sprintf("Converged with conv = %0.3f and %i iterations",conv.em,numit))
#       hasConverged = TRUE
#     }
#     if (verbose) {
#       print(sprintf("likelihood: %1.4e, conv: %1.4e",likelihood, conv.em))
#       plot(like.hist)
#     }
#   }
#   
#   if (verbose) print("Finished!")
#   
#   ldamodel$logProbW
# }


# mstep.beta.no.shrinkage <- function(ldamodel,sstats, lam, weight){
#   # estimate beta (logProbW) according to equation (7) of C.Reed's tutorial
#   mstep_beta_C(ldamodel$ntopics, ldamodel$nterms,
#                sstats$classword, ldamodel$logProbW,
#                sstats$classtotal)
#   ldamodel
# }
