// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// log_Gaussian
double log_Gaussian(arma::vec y, arma::vec pars);
RcppExport SEXP _hmmcr_log_Gaussian(SEXP ySEXP, SEXP parsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pars(parsSEXP);
    rcpp_result_gen = Rcpp::wrap(log_Gaussian(y, pars));
    return rcpp_result_gen;
END_RCPP
}
// mcmc_routine
Rcpp::List mcmc_routine(const arma::vec par_index, const arma::vec par_init, const arma::vec par_prior, const arma::vec par_unique, const arma::vec y, const arma::mat X, const arma::mat Z, const arma::mat S, const arma::vec id, const std::string model, const int nr_states, const Rcpp::List groups, const int steps, const int burnin);
RcppExport SEXP _hmmcr_mcmc_routine(SEXP par_indexSEXP, SEXP par_initSEXP, SEXP par_priorSEXP, SEXP par_uniqueSEXP, SEXP ySEXP, SEXP XSEXP, SEXP ZSEXP, SEXP SSEXP, SEXP idSEXP, SEXP modelSEXP, SEXP nr_statesSEXP, SEXP groupsSEXP, SEXP stepsSEXP, SEXP burninSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type par_index(par_indexSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type par_init(par_initSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type par_prior(par_priorSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type par_unique(par_uniqueSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type id(idSEXP);
    Rcpp::traits::input_parameter< const std::string >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const int >::type nr_states(nr_statesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< const int >::type steps(stepsSEXP);
    Rcpp::traits::input_parameter< const int >::type burnin(burninSEXP);
    rcpp_result_gen = Rcpp::wrap(mcmc_routine(par_index, par_init, par_prior, par_unique, y, X, Z, S, id, model, nr_states, groups, steps, burnin));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_hmmcr_log_Gaussian", (DL_FUNC) &_hmmcr_log_Gaussian, 2},
    {"_hmmcr_mcmc_routine", (DL_FUNC) &_hmmcr_mcmc_routine, 14},
    {NULL, NULL, 0}
};

RcppExport void R_init_hmmcr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
