#include <RcppArmadillo.h>
#include "sampler.h"
#include "sampler-leverage.h"

using namespace Rcpp;

RcppExport SEXP _stochvol_svsample_cpp(SEXP y_inSEXP, SEXP drawsSEXP, SEXP burninSEXP, SEXP X_inSEXP, SEXP bmuSEXP, SEXP BmuSEXP, SEXP a0SEXP, SEXP b0SEXP, SEXP BsigmaSEXP, SEXP thinSEXP, SEXP timethinSEXP, SEXP startpara_inSEXP, SEXP startvol_inSEXP, SEXP keeptauSEXP, SEXP quietSEXP, SEXP paraSEXP, SEXP MHstepsSEXP, SEXP B011SEXP, SEXP B022SEXP, SEXP mhcontrolSEXP, SEXP gammapriorSEXP, SEXP truncnormalSEXP, SEXP offsetSEXP, SEXP dontupdatemuSEXP, SEXP priordf_inSEXP, SEXP priorbeta_inSEXP, SEXP priorlatent0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type y_in(y_inSEXP);
    Rcpp::traits::input_parameter< const int >::type draws(drawsSEXP);
    Rcpp::traits::input_parameter< const int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type X_in(X_inSEXP);
    Rcpp::traits::input_parameter< const double >::type bmu(bmuSEXP);
    Rcpp::traits::input_parameter< const double >::type Bmu(BmuSEXP);
    Rcpp::traits::input_parameter< const double >::type a0(a0SEXP);
    Rcpp::traits::input_parameter< const double >::type b0(b0SEXP);
    Rcpp::traits::input_parameter< const double >::type Bsigma(BsigmaSEXP);
    Rcpp::traits::input_parameter< const int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< const int >::type timethin(timethinSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type startpara_in(startpara_inSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type startvol_in(startvol_inSEXP);
    Rcpp::traits::input_parameter< const bool >::type keeptau(keeptauSEXP);
    Rcpp::traits::input_parameter< const bool >::type quiet(quietSEXP);
    Rcpp::traits::input_parameter< const int >::type para(paraSEXP);
    Rcpp::traits::input_parameter< const int >::type MHsteps(MHstepsSEXP);
    Rcpp::traits::input_parameter< const double >::type B011(B011SEXP);
    Rcpp::traits::input_parameter< const double >::type B022(B022SEXP);
    Rcpp::traits::input_parameter< const double >::type mhcontrol(mhcontrolSEXP);
    Rcpp::traits::input_parameter< const bool >::type gammaprior(gammapriorSEXP);
    Rcpp::traits::input_parameter< const bool >::type truncnormal(truncnormalSEXP);
    Rcpp::traits::input_parameter< const double >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< const bool >::type dontupdatemu(dontupdatemuSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type priordf_in(priordf_inSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type priorbeta_in(priorbeta_inSEXP);
    Rcpp::traits::input_parameter< const double >::type priorlatent0(priorlatent0SEXP);
    rcpp_result_gen = Rcpp::wrap(svsample_cpp(y_in, draws, burnin, X_in, bmu, Bmu, a0, b0, Bsigma, thin, timethin, startpara_in, startvol_in, keeptau, quiet, para, MHsteps, B011, B022, mhcontrol, gammaprior, truncnormal, offset, dontupdatemu, priordf_in, priorbeta_in, priorlatent0));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _stochvol_svlsample_cpp(SEXP ySEXP, SEXP drawsSEXP, SEXP burninSEXP, SEXP designmatrixSEXP, SEXP thinparaSEXP, SEXP thinlatentSEXP, SEXP thintimeSEXP, SEXP startparaSEXP, SEXP startlatentSEXP, SEXP prior_phi_aSEXP, SEXP prior_phi_bSEXP, SEXP prior_rho_aSEXP, SEXP prior_rho_bSEXP, SEXP prior_sigma2_shapeSEXP, SEXP prior_sigma2_rateSEXP, SEXP prior_mu_muSEXP, SEXP prior_mu_sigmaSEXP, SEXP prior_beta_muSEXP, SEXP prior_beta_sigmaSEXP, SEXP verboseSEXP, SEXP offsetSEXP, SEXP stdevSEXP, SEXP gammapriorSEXP, SEXP strategySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const int >::type draws(drawsSEXP);
    Rcpp::traits::input_parameter< const int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type designmatrix(designmatrixSEXP);
    Rcpp::traits::input_parameter< const int >::type thinpara(thinparaSEXP);
    Rcpp::traits::input_parameter< const int >::type thinlatent(thinlatentSEXP);
    Rcpp::traits::input_parameter< const int >::type thintime(thintimeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type startpara(startparaSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type startlatent(startlatentSEXP);
    Rcpp::traits::input_parameter< const double >::type prior_phi_a(prior_phi_aSEXP);
    Rcpp::traits::input_parameter< const double >::type prior_phi_b(prior_phi_bSEXP);
    Rcpp::traits::input_parameter< const double >::type prior_rho_a(prior_rho_aSEXP);
    Rcpp::traits::input_parameter< const double >::type prior_rho_b(prior_rho_bSEXP);
    Rcpp::traits::input_parameter< const double >::type prior_sigma2_shape(prior_sigma2_shapeSEXP);
    Rcpp::traits::input_parameter< const double >::type prior_sigma2_rate(prior_sigma2_rateSEXP);
    Rcpp::traits::input_parameter< const double >::type prior_mu_mu(prior_mu_muSEXP);
    Rcpp::traits::input_parameter< const double >::type prior_mu_sigma(prior_mu_sigmaSEXP);
    Rcpp::traits::input_parameter< const double >::type prior_beta_mu(prior_beta_muSEXP);
    Rcpp::traits::input_parameter< const double >::type prior_beta_sigma(prior_beta_sigmaSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const double >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< const double >::type stdev(stdevSEXP);
    Rcpp::traits::input_parameter< const bool >::type gammaprior(gammapriorSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type strategy(strategySEXP);
    rcpp_result_gen = Rcpp::wrap(svlsample_cpp(y, draws, burnin, designmatrix, thinpara, thinlatent, thintime, startpara, startlatent, prior_phi_a, prior_phi_b, prior_rho_a, prior_rho_b, prior_sigma2_shape, prior_sigma2_rate, prior_mu_mu, prior_mu_sigma, prior_beta_mu, prior_beta_sigma, verbose, offset, stdev, gammaprior, strategy));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_stochvol_svsample_cpp", (DL_FUNC) &_stochvol_svsample_cpp, 27},
    {"_stochvol_svlsample_cpp", (DL_FUNC) &_stochvol_svlsample_cpp, 24},
    {NULL, NULL, 0}
};

RcppExport void R_init_stochvol(DllInfo *dll) {
 R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
 R_useDynamicSymbols(dll, FALSE);

 // registering "update" to be available for other packages
 R_RegisterCCallable("stochvol", "update_sv", (DL_FUNC) &update_sv);
 R_RegisterCCallable("stochvol", "update_svl", (DL_FUNC) &update_svl);
}

