# simulate count table given theta, beta and document length
# theta is n_docs by n_topics
# beta is topic by word
# row sums of theta and beta are 1
# doc_length is document length of each documents
# return n_docs by n_vocab simulated data

#' Simulate expression count table
#'
#' Simulate expression count table given topic by gene, topic by cell matrices and library size.
#' @export
#' @param theta Theta matrix where rows represent cells and columns represent topics.
#' @param beta Beta matrix where rows represent topics and columns represent genes.
#' @param doc_length Library size of cell samples.
#' @return A cell by gene matrix of expression counts.
#' @examples
#' \dontrun{
#' counts <- simulateCounts(theta.true, beta.true, rep(100000, nrow(theta.true)))
#' }
simulateCounts <- function(theta, beta, doc_length) {
  n_docs <- nrow(theta)
  n_topics <- ncol(theta)
  if (n_topics != nrow(beta)) {
    print("Number of topics does not agree in theta and beta matrix.")
    return(0)
  }
  n_vocab <- ncol(beta)
  if (length(doc_length)==1) {
    doc_length <- rep(doc_length, n_docs)
  } else if (n_docs != length(doc_length)) {
    print("Number of documents does not agree with doc_length.")
    return(0)
  }

  # simulate data
  docs <- matrix(0, nrow = n_docs, ncol = n_vocab)
  for (i in 1:n_docs) {
    # draw topics for each word
    tops <- rmultinom(1, doc_length[i], theta[i, ])
    # draw words
    for (j in 1:n_topics) {
      docs[i, ] <- docs[i, ] + rmultinom(1, tops[j], beta[j, ])
    }
  }
  return(docs)
}
