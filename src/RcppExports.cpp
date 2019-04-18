// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// compute_likelihood_C
double compute_likelihood_C(int model_num_topics, double model_alpha, int doc_dlength, NumericMatrix model_logProbW, NumericVector doc_words, NumericVector doc_counts, NumericMatrix phi, NumericVector var_gamma);
RcppExport SEXP _pLDA_compute_likelihood_C(SEXP model_num_topicsSEXP, SEXP model_alphaSEXP, SEXP doc_dlengthSEXP, SEXP model_logProbWSEXP, SEXP doc_wordsSEXP, SEXP doc_countsSEXP, SEXP phiSEXP, SEXP var_gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type model_num_topics(model_num_topicsSEXP);
    Rcpp::traits::input_parameter< double >::type model_alpha(model_alphaSEXP);
    Rcpp::traits::input_parameter< int >::type doc_dlength(doc_dlengthSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type model_logProbW(model_logProbWSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type doc_words(doc_wordsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type doc_counts(doc_countsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type var_gamma(var_gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_likelihood_C(model_num_topics, model_alpha, doc_dlength, model_logProbW, doc_words, doc_counts, phi, var_gamma));
    return rcpp_result_gen;
END_RCPP
}
// compute_beta_involved_likelihood_C
double compute_beta_involved_likelihood_C(int model_num_topics, int doc_dlength, NumericMatrix model_logProbW, NumericVector doc_words, NumericVector doc_counts, NumericMatrix phi);
RcppExport SEXP _pLDA_compute_beta_involved_likelihood_C(SEXP model_num_topicsSEXP, SEXP doc_dlengthSEXP, SEXP model_logProbWSEXP, SEXP doc_wordsSEXP, SEXP doc_countsSEXP, SEXP phiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type model_num_topics(model_num_topicsSEXP);
    Rcpp::traits::input_parameter< int >::type doc_dlength(doc_dlengthSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type model_logProbW(model_logProbWSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type doc_words(doc_wordsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type doc_counts(doc_countsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type phi(phiSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_beta_involved_likelihood_C(model_num_topics, doc_dlength, model_logProbW, doc_words, doc_counts, phi));
    return rcpp_result_gen;
END_RCPP
}
// updateMat
void updateMat(NumericMatrix& x);
RcppExport SEXP _pLDA_updateMat(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type x(xSEXP);
    updateMat(x);
    return R_NilValue;
END_RCPP
}
// log_sum_C
double log_sum_C(double a, double b);
RcppExport SEXP _pLDA_log_sum_C(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(log_sum_C(a, b));
    return rcpp_result_gen;
END_RCPP
}
// do_e_step_C
void do_e_step_C(int doc_length, NumericVector& oldphi, int nTopics, NumericMatrix& phi, NumericVector& digamma_gam, NumericMatrix& model_logProbW, NumericVector& docs_words, NumericVector& gammav, NumericVector& docs_counts);
RcppExport SEXP _pLDA_do_e_step_C(SEXP doc_lengthSEXP, SEXP oldphiSEXP, SEXP nTopicsSEXP, SEXP phiSEXP, SEXP digamma_gamSEXP, SEXP model_logProbWSEXP, SEXP docs_wordsSEXP, SEXP gammavSEXP, SEXP docs_countsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type doc_length(doc_lengthSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type oldphi(oldphiSEXP);
    Rcpp::traits::input_parameter< int >::type nTopics(nTopicsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type digamma_gam(digamma_gamSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type model_logProbW(model_logProbWSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type docs_words(docs_wordsSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type gammav(gammavSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type docs_counts(docs_countsSEXP);
    do_e_step_C(doc_length, oldphi, nTopics, phi, digamma_gam, model_logProbW, docs_words, gammav, docs_counts);
    return R_NilValue;
END_RCPP
}
// do_m_step_C
void do_m_step_C(int& doc_length, int nTopics, NumericVector& docs_counts, NumericMatrix& phi, NumericVector& docs_words, NumericMatrix& sstats_classword, NumericVector& sstats_classtotal);
RcppExport SEXP _pLDA_do_m_step_C(SEXP doc_lengthSEXP, SEXP nTopicsSEXP, SEXP docs_countsSEXP, SEXP phiSEXP, SEXP docs_wordsSEXP, SEXP sstats_classwordSEXP, SEXP sstats_classtotalSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int& >::type doc_length(doc_lengthSEXP);
    Rcpp::traits::input_parameter< int >::type nTopics(nTopicsSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type docs_counts(docs_countsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type docs_words(docs_wordsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type sstats_classword(sstats_classwordSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type sstats_classtotal(sstats_classtotalSEXP);
    do_m_step_C(doc_length, nTopics, docs_counts, phi, docs_words, sstats_classword, sstats_classtotal);
    return R_NilValue;
END_RCPP
}
// mstep_beta_C
void mstep_beta_C(int& model_ntopics, int& model_nterms, NumericMatrix& sstats_classword, NumericMatrix& model_logProbW, NumericVector& sstats_classtotal);
RcppExport SEXP _pLDA_mstep_beta_C(SEXP model_ntopicsSEXP, SEXP model_ntermsSEXP, SEXP sstats_classwordSEXP, SEXP model_logProbWSEXP, SEXP sstats_classtotalSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int& >::type model_ntopics(model_ntopicsSEXP);
    Rcpp::traits::input_parameter< int& >::type model_nterms(model_ntermsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type sstats_classword(sstats_classwordSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type model_logProbW(model_logProbWSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type sstats_classtotal(sstats_classtotalSEXP);
    mstep_beta_C(model_ntopics, model_nterms, sstats_classword, model_logProbW, sstats_classtotal);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_pLDA_compute_likelihood_C", (DL_FUNC) &_pLDA_compute_likelihood_C, 8},
    {"_pLDA_compute_beta_involved_likelihood_C", (DL_FUNC) &_pLDA_compute_beta_involved_likelihood_C, 6},
    {"_pLDA_updateMat", (DL_FUNC) &_pLDA_updateMat, 1},
    {"_pLDA_log_sum_C", (DL_FUNC) &_pLDA_log_sum_C, 2},
    {"_pLDA_do_e_step_C", (DL_FUNC) &_pLDA_do_e_step_C, 9},
    {"_pLDA_do_m_step_C", (DL_FUNC) &_pLDA_do_m_step_C, 7},
    {"_pLDA_mstep_beta_C", (DL_FUNC) &_pLDA_mstep_beta_C, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_pLDA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}