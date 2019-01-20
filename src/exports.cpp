#include "update_functions.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

// svsample_cpp
Rcpp::List svsample_cpp(const arma::vec& y_in, const int draws, const int burnin, const arma::mat& X_in, const double bmu, const double Bmu, const double a0, const double b0, const double Bsigma, const int thin, const int timethin, const Rcpp::List& startpara_in, const arma::vec& startvol_in, const bool keeptau, const bool quiet, const int para, const int MHsteps, const double B011, const double B022, const double mhcontrol, const bool gammaprior, const bool truncnormal, const double offset, const bool dontupdatemu, const arma::vec& priordf_in, const arma::vec& priorbeta_in, const double priorlatent0);
RcppExport SEXP _stochvol_svsample_cpp(SEXP y_inSEXP, SEXP drawsSEXP, SEXP burninSEXP, SEXP X_inSEXP, SEXP bmuSEXP, SEXP BmuSEXP, SEXP a0SEXP, SEXP b0SEXP, SEXP BsigmaSEXP, SEXP thinSEXP, SEXP timethinSEXP, SEXP startpara_inSEXP, SEXP startvol_inSEXP, SEXP keeptauSEXP, SEXP quietSEXP, SEXP paraSEXP, SEXP MHstepsSEXP, SEXP B011SEXP, SEXP B022SEXP, SEXP mhcontrolSEXP, SEXP gammapriorSEXP, SEXP truncnormalSEXP, SEXP offsetSEXP, SEXP dontupdatemuSEXP, SEXP priordf_inSEXP, SEXP priorbeta_inSEXP, SEXP priorlatent0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y_in(y_inSEXP);
    Rcpp::traits::input_parameter< const int >::type draws(drawsSEXP);
    Rcpp::traits::input_parameter< const int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X_in(X_inSEXP);
    Rcpp::traits::input_parameter< const double >::type bmu(bmuSEXP);
    Rcpp::traits::input_parameter< const double >::type Bmu(BmuSEXP);
    Rcpp::traits::input_parameter< const double >::type a0(a0SEXP);
    Rcpp::traits::input_parameter< const double >::type b0(b0SEXP);
    Rcpp::traits::input_parameter< const double >::type Bsigma(BsigmaSEXP);
    Rcpp::traits::input_parameter< const int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< const int >::type timethin(timethinSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type startpara_in(startpara_inSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type startvol_in(startvol_inSEXP);
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
    Rcpp::traits::input_parameter< const arma::vec& >::type priordf_in(priordf_inSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type priorbeta_in(priorbeta_inSEXP);
    Rcpp::traits::input_parameter< const double >::type priorlatent0(priorlatent0SEXP);
    rcpp_result_gen = Rcpp::wrap(svsample_cpp(y_in, draws, burnin, X_in, bmu, Bmu, a0, b0, Bsigma, thin, timethin, startpara_in, startvol_in, keeptau, quiet, para, MHsteps, B011, B022, mhcontrol, gammaprior, truncnormal, offset, dontupdatemu, priordf_in, priorbeta_in, priorlatent0));
    return rcpp_result_gen;
END_RCPP
}
// svlsample_cpp
Rcpp::List svlsample_cpp(const arma::vec& y, const int draws, const int burnin, const arma::mat& X, const int thinpara, const int thinlatent, const int thintime, const Rcpp::List& theta_init, const arma::vec& h_init, const double prior_phi_a, const double prior_phi_b, const double prior_rho_a, const double prior_rho_b, const double prior_sigma2_shape, const double prior_sigma2_rate, const double prior_mu_mu, const double prior_mu_sigma, const double prior_beta_mu, const double prior_beta_sigma, const bool verbose, const double offset, const double stdev, const bool gammaprior, const Rcpp::CharacterVector& strategy);
RcppExport SEXP _stochvol_svlsample_cpp(SEXP ySEXP, SEXP drawsSEXP, SEXP burninSEXP, SEXP XSEXP, SEXP thinparaSEXP, SEXP thinlatentSEXP, SEXP thintimeSEXP, SEXP theta_initSEXP, SEXP h_initSEXP, SEXP prior_phi_aSEXP, SEXP prior_phi_bSEXP, SEXP prior_rho_aSEXP, SEXP prior_rho_bSEXP, SEXP prior_sigma2_shapeSEXP, SEXP prior_sigma2_rateSEXP, SEXP prior_mu_muSEXP, SEXP prior_mu_sigmaSEXP, SEXP prior_beta_muSEXP, SEXP prior_beta_sigmaSEXP, SEXP verboseSEXP, SEXP offsetSEXP, SEXP stdevSEXP, SEXP gammapriorSEXP, SEXP strategySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const int >::type draws(drawsSEXP);
    Rcpp::traits::input_parameter< const int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const int >::type thinpara(thinparaSEXP);
    Rcpp::traits::input_parameter< const int >::type thinlatent(thinlatentSEXP);
    Rcpp::traits::input_parameter< const int >::type thintime(thintimeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type theta_init(theta_initSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type h_init(h_initSEXP);
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
    rcpp_result_gen = Rcpp::wrap(svlsample_cpp(y, draws, burnin, X, thinpara, thinlatent, thintime, theta_init, h_init, prior_phi_a, prior_phi_b, prior_rho_a, prior_rho_b, prior_sigma2_shape, prior_sigma2_rate, prior_mu_mu, prior_mu_sigma, prior_beta_mu, prior_beta_sigma, verbose, offset, stdev, gammaprior, strategy));
    return rcpp_result_gen;
END_RCPP
}

// validate (ensure exported C++ functions exist before calling them)
static int _stochvol_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("void(*update_sv)(const arma::vec&,arma::vec&,arma::vec&,double&,arma::vec&,arma::ivec&,const bool,const double,const double,const double,const double,const double,const double,const double,const double,const double,const bool,const bool,const double,const int,const int,const bool,const double)");
        signatures.insert("void(*update_svl)(const arma::vec&,const arma::vec&,const arma::vec&,double&,double&,double&,double&,arma::vec&,arma::vec&,const arma::vec&,const arma::vec&,const arma::vec&,const arma::vec&,const double,const bool,const arma::ivec&)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _stochvol_RcppExport_registerCCallable() { 
    R_RegisterCCallable("stochvol", "update_sv", (DL_FUNC)update_sv);
    R_RegisterCCallable("stochvol", "update_svl", (DL_FUNC)update_svl);
    R_RegisterCCallable("stochvol", "_stochvol_RcppExport_validate", (DL_FUNC)_stochvol_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_stochvol_svsample_cpp", (DL_FUNC) &_stochvol_svsample_cpp, 27},
    {"_stochvol_svlsample_cpp", (DL_FUNC) &_stochvol_svlsample_cpp, 24},
    {"_stochvol_RcppExport_registerCCallable", (DL_FUNC) &_stochvol_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_stochvol(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
