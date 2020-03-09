// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// cpp_IPOT
Rcpp::List cpp_IPOT(arma::vec mu, arma::vec nu, arma::mat C, double beta, int L, int maxiter, double abstol, bool printer);
RcppExport SEXP _T4transport_cpp_IPOT(SEXP muSEXP, SEXP nuSEXP, SEXP CSEXP, SEXP betaSEXP, SEXP LSEXP, SEXP maxiterSEXP, SEXP abstolSEXP, SEXP printerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type C(CSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type abstol(abstolSEXP);
    Rcpp::traits::input_parameter< bool >::type printer(printerSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_IPOT(mu, nu, C, beta, L, maxiter, abstol, printer));
    return rcpp_result_gen;
END_RCPP
}
// cpp_Sinkhorn_Bregman
Rcpp::List cpp_Sinkhorn_Bregman(arma::vec mu, arma::vec nu, arma::mat costm, double lambda, int maxiter, double abstol, bool printer);
RcppExport SEXP _T4transport_cpp_Sinkhorn_Bregman(SEXP muSEXP, SEXP nuSEXP, SEXP costmSEXP, SEXP lambdaSEXP, SEXP maxiterSEXP, SEXP abstolSEXP, SEXP printerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type costm(costmSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type abstol(abstolSEXP);
    Rcpp::traits::input_parameter< bool >::type printer(printerSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_Sinkhorn_Bregman(mu, nu, costm, lambda, maxiter, abstol, printer));
    return rcpp_result_gen;
END_RCPP
}
// cpp_Sinkhorn_Cuturi
Rcpp::List cpp_Sinkhorn_Cuturi(arma::vec a, arma::vec b, arma::mat costm, double lambda, int maxiter, double abstol, bool printer);
RcppExport SEXP _T4transport_cpp_Sinkhorn_Cuturi(SEXP aSEXP, SEXP bSEXP, SEXP costmSEXP, SEXP lambdaSEXP, SEXP maxiterSEXP, SEXP abstolSEXP, SEXP printerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type costm(costmSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type abstol(abstolSEXP);
    Rcpp::traits::input_parameter< bool >::type printer(printerSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_Sinkhorn_Cuturi(a, b, costm, lambda, maxiter, abstol, printer));
    return rcpp_result_gen;
END_RCPP
}
// extra_distmat
arma::mat extra_distmat(arma::mat& X);
RcppExport SEXP _T4transport_extra_distmat(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(extra_distmat(X));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_T4transport_cpp_IPOT", (DL_FUNC) &_T4transport_cpp_IPOT, 8},
    {"_T4transport_cpp_Sinkhorn_Bregman", (DL_FUNC) &_T4transport_cpp_Sinkhorn_Bregman, 7},
    {"_T4transport_cpp_Sinkhorn_Cuturi", (DL_FUNC) &_T4transport_cpp_Sinkhorn_Cuturi, 7},
    {"_T4transport_extra_distmat", (DL_FUNC) &_T4transport_extra_distmat, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_T4transport(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
