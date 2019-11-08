# This file contains all functions for penalized LDA
# main function is plda(x, k, lambda)
# x data matrix, row is docs, col is genes
# k number of topics
# lambda shrinkage parameter

# # make sure it can run in R markdown
# library(devtools)
# find_rtools()

## beta shrinkage
# calculate gradiant
d1 <- function(beta, lambda, a, weight) {
  k <- nrow(beta)
  V <- ncol(beta)
  b <- as.vector(beta)
  weight <- rep(weight, each=k)

  b[b==0] <- 1e-100

  cmeans <- colMeans(beta)
  cmeans.diff <- beta - matrix(rep(cmeans, k), nrow = k, ncol = V, byrow = T)
  # return
  as.vector(lambda * 2 * cmeans.diff) * weight - as.vector(a)/b # + 1/(t*b)

}

d2 <- function(beta, lambda, a, weight) {
  k <- nrow(beta)
  V <- ncol(beta)
  b <- as.vector(beta)
  block <- matrix(-2/k, nrow = k, ncol = k) + diag(2, k)
  res <- kronecker(Diagonal(V, lambda * weight), block)
  diff <- Diagonal(length(b), as.vector(a)/b^2)

  res + diff
}

# newton decrement
newton.decrement <- function(delta, d2f) {
  return(t(delta) %*% d2f %*% delta)
}

# solve Cx=d
oneRun <- function(beta, lambda, a, weight) {
  d1f <- d1(beta, lambda, a, weight)
  d2f <- d2(beta, lambda, a, weight)
  # get restirction matrix A
  k <- nrow(beta)
  V <- ncol(beta)
  A <- matrix(0, nrow = k, ncol = k*V)
  for (i in 1:k) {
    A[i, (0:(V-1))*k+i] <- 1
  }

  C <- cbind(rbind(d2f, A), rbind(t(A), matrix(0, nrow = k, ncol = k)))
  d <- c(-d1f, rep(0, k))
  x <- solve(C, d)
  delta <- x[1:(k*V)]
  while (min(as.vector(beta) + delta) < 0) {
    delta <- delta/2
    # print('rej')
  }
  beta.new <- matrix(as.vector(beta) + delta, nrow = k, ncol = V, byrow = F)
  n.d <- newton.decrement(delta, d2f)
  # print(paste0('Newton Decrement: ', n.d))
  return(list(beta.new, n.d))
}

est.b.shrinkage <- function(b.start, lambda, a, weight) {
  # beta shrinkage
  n.d <- 1
  itr <- 1
  while (n.d>1e-20 & itr<100) {
    b.next <- oneRun(b.start, lambda, a, weight)
    b.start <- b.next[[1]]
    itr <- itr + 1
    n.d <- b.next[[2]][1]
    if (is.nan(n.d)) {
      n.d <- 0
    }
  }
  return(b.next[[1]])
}

mstep.beta <- function(ldamodel,sstats, lam, weight){
  # estimate beta (logProbW) according to equation (7) of C.Reed's tutorial
  mstep_beta_C(ldamodel$ntopics, ldamodel$nterms,
               sstats$classword, ldamodel$logProbW,
               sstats$classtotal)
  b.start <- exp(ldamodel$logProbW)
  ldamodel$logProbW <- log(est.b.shrinkage(b.start, lam, a=sstats$classword, weight))
  ldamodel
}

mstep.beta.start <- function(ldamodel,sstats, lam){
  # estimate beta (logProbW) according to equation (7) of C.Reed's tutorial
  mstep_beta_C(ldamodel$ntopics, ldamodel$nterms,
               sstats$classword, ldamodel$logProbW,
               sstats$classtotal)
  ldamodel
}



################################################################
################################################################
#' Penalized Latent Dirichlet Allocation
#'
#' The main function that runs penalized LDA.
#' @export
#' @param x Input matrix. A matrix of positive integers where rows represents cells (documents) and column represents genes (words).
#' @param k Number of topics.
#' @param lambda Shrinkage parameter.
#' @param warm.start.seed Seed of LDA estimation used as pLDA warm start. pLDA uses parameters estimated from LDA as a warm start to boost speed.
#' @param verbose Default FALSE. Print progress.
#' @param warm.start.lda Default NULL. LDA object for warm start.
#' @param warm.start.log.beta Default NULL. log beta matrix for warm start.
#' @param pdata Default NULL. Table of phenotypic data.
#' @return plda() returns a list of the penalized LDA output. logProbW is the log of topic by gene matrix beta. gamma is the cell by topic matrix of topic frequencies for each cell. genes are gene names extracted from the column name of input count matrix.
#' @examples
#' \dontrun{
#' plda(x=cell_by_gene_expr_matrix, k=10, lambda=10^3)
#' }
plda <- function(x, 
                 k, 
                 lambda, 
                 warm.start.seed=2017, 
                 verbose = FALSE, 
                 warm.start.lda=NULL, 
                 warm.start.log.beta=NULL,
                 pdata = NULL) {
  n_vocab <- ncol(x)
  n_docs <- nrow(x)
  genes <- colnames(x)
  sample.names <- rownames(x)
  tot.lib.size.per.million <- sum(x) / 1e6
  lambda.base <- lambda
  lambda <- lambda * tot.lib.size.per.million

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
  alpha = 50/k

  # a warm start
  # run topicmodel LDA and use their output as our input
  if (!is.null(warm.start.lda)) {
    log.beta.warm.start <- warm.start.lda@beta
  } else if (!is.null(warm.start.log.beta)) {
    log.beta.warm.start <- warm.start.log.beta
  } else {
    lda_fit <- topicmodels::LDA(x, k, control = list(seed=warm.start.seed*2))
    log.beta.warm.start <- lda_fit@beta
  }

  # weights in beta penalty term
  # genes with large mean has smaller weight as they tend to have more variation
  weight.beta <- apply(x, 2, mean)
  weight.beta <- weight.beta / sum(weight.beta)
  weight.beta <- 1/weight.beta^2

  # init the model randomly
  cwinit = matrix(runif(k*corpus$nterms),k,corpus$nterms) + 1/corpus$nterms
  ldamodel = list(logProbW=matrix(rep(0,k*corpus$nterms),k,corpus$nterms), alpha = 1, ntopics=k, nterms=corpus$nterms)
  sstats = list(ndocs=0,classword=cwinit,k,corpus$nterms,classtotal=rowSums(cwinit), alpha.ss = 0)
  ldamodel = mstep.beta.start(ldamodel,sstats, 0)
  ldamodel$alpha = alpha
  # initialize the beta matrix with warm start
  ldamodel$logProbW <- log.beta.warm.start

  if (verbose) {
    print("Finished warm start.")
  }

  like.hist = c() # keep track of the likelihood
  likelihood.prev = 0
  numit = 0
  hasConverged = FALSE
  nTopics = ldamodel$ntopics
  nTerms = ldamodel$nterms
  phi = matrix(rep(0,nTopics*nTerms),nTopics, nTerms)

  # # # # Run variational expectation-maximization # # # #
  while (!hasConverged){
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
    likelihood.beta = 0

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
        do_e_step_C(doc.length,
                    oldphi,
                    nTopics,
                    phi,
                    digamma.gam,
                    ldamodel$logProbW,
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

      do_m_step_C(doc.length,
                  nTopics,
                  corpus$docs[[d]]$counts,
                  phi,
                  corpus$docs[[d]]$words,
                  sstats$classword,
                  sstats$classtotal)

      sstats$ndocs = sstats$ndocs + 1
      likelihood = likelihood + doc.like
      likelihood.beta = likelihood.beta + compute_beta_involved_likelihood_C(ldamodel$ntopics,
                                                                             corpus$docs[[d]]$dlength,
                                                                             ldamodel$logProbW,
                                                                             corpus$docs[[d]]$words,
                                                                             corpus$docs[[d]]$counts,
                                                                             phi)
      gammaOut[d, ] <- gammav
    } # end for each document

    # # # do M-step # # #
    # estimate beta
    ldamodel = mstep.beta(ldamodel, sstats, lambda, weight.beta)

    # estimate alpha
    if (estAlpha){
      D = sstats$ndocs
      alpha.iter = 0
      a.init = 100
      log.a = log(a.init)
      alpha.hasconv = FALSE
      while (!alpha.hasconv){
        alpha.iter = alpha.iter + 1
        a = exp(log.a)

        if (is.nan(a)){
          a.init = a.init*10
          print(sprintf("alpha became nan, initializing with alpha = %1.3e",a.init))
          a = a.init
          log.a = log(a)
        }

        f = D*(lgamma(nTopics*a) - nTopics*lgamma(a)) + (a-1)*sstats$alpha.ss
        df = D * (nTopics*digamma(nTopics*a) - nTopics*digamma(a)) + sstats$alpha.ss
        d2f = D * (nTopics*nTopics*trigamma(nTopics*a) - nTopics*trigamma(a))
        log.a = log.a - df/(d2f*a + df)
        # print(sprintf("alpha optimization: %1.3e  %1.3e   %1.3e", exp(log.a), f, df))
        if (abs(df) < 1e-5 || alpha.iter > 100){
          alpha.hasconv = TRUE
        }
      }
      ldamodel$alpha = exp(log.a)
    }

    # add penalty term
    beta <- exp(ldamodel$logProbW)
    tmp <- sweep(beta, 2, apply(beta, 2, mean))
    likelihood.penality <- lambda * sum(apply(tmp^2, 2, sum) * weight.beta)
    likelihood.docs <- likelihood
    likelihood <- likelihood.docs - likelihood.penality

    conv.em = (likelihood.prev - likelihood)/likelihood.prev
    likelihood.prev = likelihood
    like.hist[numit] = likelihood

    # make sure we're iterating enough for the likelihood to converge'
    if (conv.em < 0){
      MAX.ES.IT = MAX.ES.IT*2
    }
    # if (((conv.em < EM.CONV && conv.em > 0)  || numit > MAX.EM) && numit > 2){
    if (((abs(conv.em) < EM.CONV)  || numit > MAX.EM) && numit > 2){
      # print(sprintf("Converged with conv = %0.3f and %i iterations",conv.em,numit))
      hasConverged = TRUE
    }
    if (verbose) {
      print(sprintf("likelihood: %1.4e, conv: %1.4e",likelihood, conv.em))
      plot(like.hist)
    }

  }
  if (verbose) print("Finished!")

  # calculate topic document matrix, gamma
  ldamodel$gamma <- gammaOut/rowSums(gammaOut)

  # genes
  ldamodel$genes <- genes

  # decompose likelihoods
  ldamodel$likelihood.penality <- likelihood.penality
  ldamodel$likelihood.beta <- likelihood.beta
  ldamodel$likelihood.docs <- likelihood.docs
  
  # lambda value
  ldamodel$lambda.base <- lambda.base
  ldamodel$lambda <- lambda
  
  # add gene names to column names of beta
  colnames(ldamodel$logProbW) <- genes
  
  # add pdata to output
  ldamodel$pdata <- pdata
  
  # add column and row names to beta and gamma
  colnames(ldamodel$logProbW) <- genes
  rownames(ldamodel$gamma) <- sample.names
  rownames(ldamodel$logProbW) <- colnames(ldamodel$gamma) <- paste0('topic_', 1:nrow(ldamodel$logProbW))
  
  return(ldamodel)
}
