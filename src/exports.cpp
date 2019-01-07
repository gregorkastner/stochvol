#include "sampler.h"
#include "run-samplers.h"
#include <Rcpp.h>

using namespace Rcpp;

RcppExport SEXP _stochvol_sampler(const SEXP y_in, const SEXP draws_in,
  const SEXP burnin_in, const SEXP X_in,
  const SEXP bmu_in, const SEXP Bmu_in,
  const SEXP a0_in, const SEXP b0_in, const SEXP Bsigma_in,
  const SEXP thin_in, const SEXP timethin_in, const SEXP startpara_in,
  const SEXP startvol_in, const SEXP keeptau_in,
  const SEXP quiet_in, const SEXP para_in,
  const SEXP MHsteps_in, const SEXP B011_in, const SEXP B022_in,
  const SEXP mhcontrol_in, const SEXP gammaprior_in,
  const SEXP truncnormal_in, const SEXP offset_in,
  const SEXP dontupdatemu_in, const SEXP priordf_in,
  const SEXP priorbeta_in, const SEXP priorlatent0_in) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(sampler(y_in, draws_in,
          burnin_in, X_in,
          bmu_in, Bmu_in,
          a0_in, b0_in, Bsigma_in,
          thin_in, timethin_in, startpara_in,
          startvol_in, keeptau_in,
          quiet_in, para_in,
          MHsteps_in, B011_in, B022_in,
          mhcontrol_in, gammaprior_in,
          truncnormal_in, offset_in,
          dontupdatemu_in, priordf_in,
          priorbeta_in, priorlatent0_in));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _stochvol_svlsample_cpp(SEXP drawsSEXP, SEXP ySEXP, SEXP y_starSEXP, SEXP dSEXP, SEXP burninSEXP, SEXP thinparaSEXP, SEXP thinlatentSEXP, SEXP thintimeSEXP, SEXP phi_initSEXP, SEXP rho_initSEXP, SEXP sigma2_initSEXP, SEXP mu_initSEXP, SEXP h_initSEXP, SEXP prior_phi_aSEXP, SEXP prior_phi_bSEXP, SEXP prior_rho_aSEXP, SEXP prior_rho_bSEXP, SEXP prior_sigma2_shapeSEXP, SEXP prior_sigma2_rateSEXP, SEXP prior_mu_muSEXP, SEXP prior_mu_sigmaSEXP, SEXP verboseSEXP, SEXP stdevSEXP, SEXP gammapriorSEXP, SEXP strategySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type draws(drawsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type y_star(y_starSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type d(dSEXP);
    Rcpp::traits::input_parameter< const int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< const int >::type thinpara(thinparaSEXP);
    Rcpp::traits::input_parameter< const int >::type thinlatent(thinlatentSEXP);
    Rcpp::traits::input_parameter< const int >::type thintime(thintimeSEXP);
    Rcpp::traits::input_parameter< const double >::type phi_init(phi_initSEXP);
    Rcpp::traits::input_parameter< const double >::type rho_init(rho_initSEXP);
    Rcpp::traits::input_parameter< const double >::type sigma2_init(sigma2_initSEXP);
    Rcpp::traits::input_parameter< const double >::type mu_init(mu_initSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type h_init(h_initSEXP);
    Rcpp::traits::input_parameter< const double >::type prior_phi_a(prior_phi_aSEXP);
    Rcpp::traits::input_parameter< const double >::type prior_phi_b(prior_phi_bSEXP);
    Rcpp::traits::input_parameter< const double >::type prior_rho_a(prior_rho_aSEXP);
    Rcpp::traits::input_parameter< const double >::type prior_rho_b(prior_rho_bSEXP);
    Rcpp::traits::input_parameter< const double >::type prior_sigma2_shape(prior_sigma2_shapeSEXP);
    Rcpp::traits::input_parameter< const double >::type prior_sigma2_rate(prior_sigma2_rateSEXP);
    Rcpp::traits::input_parameter< const double >::type prior_mu_mu(prior_mu_muSEXP);
    Rcpp::traits::input_parameter< const double >::type prior_mu_sigma(prior_mu_sigmaSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const double >::type stdev(stdevSEXP);
    Rcpp::traits::input_parameter< const bool >::type gammaprior(gammapriorSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type strategy(strategySEXP);
    rcpp_result_gen = Rcpp::wrap(svlsample_cpp(draws, y, y_star, d, burnin, thinpara, thinlatent, thintime, phi_init, rho_init, sigma2_init, mu_init, h_init, prior_phi_a, prior_phi_b, prior_rho_a, prior_rho_b, prior_sigma2_shape, prior_sigma2_rate, prior_mu_mu, prior_mu_sigma, verbose, stdev, gammaprior, strategy));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_stochvol_sampler", (DL_FUNC) &sampler, 27},
    {"_stochvol_svlsample_cpp", (DL_FUNC) &_stochvol_svlsample_cpp, 25},
    {NULL, NULL, 0}
};

RcppExport void R_init_stochvol(DllInfo *dll) {
 R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
 R_useDynamicSymbols(dll, FALSE);

 // registering "update" to be available for other packages
 R_RegisterCCallable("stochvol", "update", (DL_FUNC) &update);
 R_RegisterCCallable("stochvol", "update_leverage", (DL_FUNC) &update_leverage);
}

