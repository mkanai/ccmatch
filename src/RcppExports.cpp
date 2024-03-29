// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// ccmatch
DataFrame ccmatch(NumericMatrix x, int N, bool display_progress);
RcppExport SEXP ccmatch_ccmatch(SEXP xSEXP, SEXP NSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    __result = Rcpp::wrap(ccmatch(x, N, display_progress));
    return __result;
END_RCPP
}
