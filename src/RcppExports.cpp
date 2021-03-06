// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// Rginv
arma::mat Rginv(const arma::mat& m);
RcppExport SEXP _autostsm_Rginv(SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(Rginv(m));
    return rcpp_result_gen;
END_RCPP
}
// gen_inv
arma::mat gen_inv(arma::mat& m);
RcppExport SEXP _autostsm_gen_inv(SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(gen_inv(m));
    return rcpp_result_gen;
END_RCPP
}
// kalman_filter
Rcpp::List kalman_filter(Rcpp::List& sp, const arma::mat& yt, const arma::mat& X, bool smooth);
RcppExport SEXP _autostsm_kalman_filter(SEXP spSEXP, SEXP ytSEXP, SEXP XSEXP, SEXP smoothSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type sp(spSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type yt(ytSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< bool >::type smooth(smoothSEXP);
    rcpp_result_gen = Rcpp::wrap(kalman_filter(sp, yt, X, smooth));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_autostsm_Rginv", (DL_FUNC) &_autostsm_Rginv, 1},
    {"_autostsm_gen_inv", (DL_FUNC) &_autostsm_gen_inv, 1},
    {"_autostsm_kalman_filter", (DL_FUNC) &_autostsm_kalman_filter, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_autostsm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
