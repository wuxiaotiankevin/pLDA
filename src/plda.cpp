#include <Rcpp.h>
#include <math.h>
#include <float.h>
using namespace Rcpp;

double trigamma(double x);
double digamma(double x);
double log_gamma(double x);
double log_sum_C(double a, double b);

// [[Rcpp::export]]
double compute_likelihood_C(int model_num_topics, 
                            double model_alpha, 
                            int doc_dlength, 
                            NumericMatrix model_logProbW, 
                            NumericVector doc_words, 
                            NumericVector doc_counts, 
                            NumericMatrix phi, 
                            NumericVector var_gamma)
{
  double likelihood = 0, digsum = 0, var_gamma_sum = 0, dig[model_num_topics];
  int k, n;
  
  for (k = 0; k < model_num_topics; k++)
  {
    dig[k] = digamma(var_gamma[k]);
    var_gamma_sum += var_gamma[k];
  }
  digsum = digamma(var_gamma_sum);
  
  likelihood =
    lgamma(model_alpha * model_num_topics)
    - model_num_topics * lgamma(model_alpha)
    - (lgamma(var_gamma_sum));
    for (k = 0; k < model_num_topics; k++)
    {
      likelihood +=
        (model_alpha - 1)*(dig[k] - digsum) + lgamma(var_gamma[k])
      - (var_gamma[k] - 1)*(dig[k] - digsum);
      
      for (n = 0; n < doc_dlength; n++)
      {
        if (phi(n, k) > 0)
        {
          likelihood += doc_counts[n] *
            (phi(n, k) * ((dig[k] - digsum) - log(phi(n, k))
                            + model_logProbW(k, doc_words[n]) ));
        }
      }
    }
    return(likelihood);
}

// [[Rcpp::export]]
double compute_beta_involved_likelihood_C(int model_num_topics, 
                                          int doc_dlength, 
                                          NumericMatrix model_logProbW, 
                                          NumericVector doc_words, 
                                          NumericVector doc_counts, 
                                          NumericMatrix phi)
{
  double likelihood = 0;
  int k, n;
  for (k = 0; k < model_num_topics; k++)
  {
    for (n = 0; n < doc_dlength; n++)
    {
      if (phi(n, k) > 0)
      {
        likelihood += doc_counts[n] * phi(n, k) * model_logProbW(k, doc_words[n]);
      }
    }
  }
  return(likelihood);
}


// [[Rcpp::export]]
void updateMat(NumericMatrix &x) {
  x(0, 0) = 1;
}

// [[Rcpp::export]]
double log_sum_C(double a, double b) {
  double small, large;
  if (a > b) {
    small = b;
    large = a;
  } else {
    small = a;
    large = b;
  }
  return(large + log(1 + exp(small - large)));
}

// [[Rcpp::export]]
void do_e_step_C(int doc_length, NumericVector &oldphi,
                 int nTopics, NumericMatrix &phi,
                 NumericVector &digamma_gam,
                 NumericMatrix &model_logProbW,
                 NumericVector &docs_words,
                 NumericVector &gammav,
                 NumericVector &docs_counts
) {
  int n, k;
  double phisum;
  
  for (n = 0; n < doc_length; n++) {
    phisum = 0;
    for (k = 0; k < nTopics; k++) {
      oldphi[k] = phi(n,k);
      phi(n,k) = digamma_gam[k] + model_logProbW(k, docs_words[n]); //corpus$docs[[d]]$words[[n]]+1
      
      if (k > 0) {
        phisum = log_sum_C(phisum, phi(n, k));
      } else {
        phisum = phi(n, k);
      }
    }
    
    for (k = 0; k < nTopics; k++) {
      phi(n,k) = exp(phi(n,k) - phisum);
      gammav[k] = gammav[k] + docs_counts[n] * (phi(n,k) - oldphi[k]);
      digamma_gam[k] = digamma(gammav[k]);
      /*
       if (is.na(gammav[k])){
       print('error with gammav, contains na')
       browser()
       }
       */
    }
  }
}

// [[Rcpp::export]]
void do_m_step_C(int &doc_length, int nTopics,
                 NumericVector &docs_counts,
                 NumericMatrix &phi,
                 NumericVector &docs_words,
                 NumericMatrix &sstats_classword,
                 NumericVector &sstats_classtotal) {
  int n, k;
  double phink;
  
  for (n = 0; n < doc_length; n++) {
    for (k = 0; k < nTopics; k++) {
      phink = docs_counts[n]*phi(n,k);
      sstats_classword(k,docs_words[n]) = sstats_classword(k, docs_words[n]) + phink;
      sstats_classtotal[k] = sstats_classtotal[k] + phink;
    }
  }
}

/*
 for (n in 1:doc.length ){
 for (k in 1:nTopics){
 phink = corpus$docs[[d]]$counts[n]*phi[n,k]
 sstats$classword[k,corpus$docs[[d]]$words[n] + 1] = sstats$classword[k,corpus$docs[[d]]$words[n] + 1] + phink
 sstats$classtotal[k] = sstats$classtotal[k] + phink
 }
 }
 */

// [[Rcpp::export]]
void mstep_beta_C(int &model_ntopics, int &model_nterms,
                  NumericMatrix &sstats_classword, 
                  NumericMatrix &model_logProbW,
                  NumericVector &sstats_classtotal) {
  int k, w;
  
  for (k = 0; k < model_ntopics; k++) {
    for(w =0; w < model_nterms; w++) {
      if (sstats_classword(k,w) > 0 ){
        model_logProbW(k,w) = log(sstats_classword(k,w)) - log(sstats_classtotal[k]);
      }
      else{
        model_logProbW(k,w) = -100;
      }
    }
  }
}

/*
# estimate beta (logProbW) according to equation (7) of C.Reed's tutorial
 for (k in 1:ldamodel$ntopics){
 for (w in 1:ldamodel$nterms){
 if (sstats$classword[k,w] > 0 ){
 ldamodel$logProbW[k,w] = log(sstats$classword[k,w]) - log(sstats$classtotal[k])
 }
 else{
 ldamodel$logProbW[k,w] = -100
 }
 }
 }
 */

double trigamma(double x)
{
  double p;
  int i;
  
  x=x+6;
  p=1/(x*x);
  p=(((((0.075757575757576*p-0.033333333333333)*p+0.0238095238095238)
         *p-0.033333333333333)*p+0.166666666666667)*p+1)/x+0.5*p;
  for (i=0; i<6 ;i++)
  {
    x=x-1;
    p=1/(x*x)+p;
  }
  return(p);
}


/*
 * taylor approximation of first derivative of the log gamma function
 *
 */
double digamma(double x)
{
  double p;
  x=x+6;
  p=1/(x*x);
  p=(((0.004166666666667*p-0.003968253986254)*p+
    0.008333333333333)*p-0.083333333333333)*p;
  p=p+log(x)-0.5/x-1/(x-1)-1/(x-2)-1/(x-3)-1/(x-4)-1/(x-5)-1/(x-6);
  return p;
}

double log_gamma(double x)
{
  double z=1/(x*x);
  
  x=x+6;
  z=(((-0.000595238095238*z+0.000793650793651)
        *z-0.002777777777778)*z+0.083333333333333)/x;
  z=(x-0.5)*log(x)-x+0.918938533204673+z-log(x-1)-
  log(x-2)-log(x-3)-log(x-4)-log(x-5)-log(x-6);
  return z;
}