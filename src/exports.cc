#include "single_update.h"
#include "sampling_main.h"
#include "utils_latent_states.h"
#include <RcppArmadillo.h>
#include <string>
#include <set>

using namespace Rcpp;

RcppExport SEXP _stochvol_svsample_cpp(
    SEXP y_inSEXP,
    SEXP drawsSEXP,
    SEXP burninSEXP,
    SEXP XSEXP,
    SEXP priorspec_inSEXP,
    SEXP thinparaSEXP,
    SEXP thinlatentSEXP,
    SEXP keeptime_inSEXP,
    SEXP startparaSEXP,
    SEXP startlatentSEXP,
    SEXP keeptauSEXP,
    SEXP quietSEXP,
    SEXP correct_model_specificationSEXP,
    SEXP interweaveSEXP,
    SEXP offsetSEXP,
    SEXP expert_inSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y_in(y_inSEXP);
    Rcpp::traits::input_parameter< const int >::type draws(drawsSEXP);
    Rcpp::traits::input_parameter< const int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type priorspec_in(priorspec_inSEXP);
    Rcpp::traits::input_parameter< const int >::type thinpara(thinparaSEXP);
    Rcpp::traits::input_parameter< const int >::type thinlatent(thinlatentSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type keeptime_in(keeptime_inSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type startpara(startparaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type startlatent(startlatentSEXP);
    Rcpp::traits::input_parameter< const bool >::type keeptau(keeptauSEXP);
    Rcpp::traits::input_parameter< const bool >::type quiet(quietSEXP);
    Rcpp::traits::input_parameter< const bool >::type correct_model_specification(correct_model_specificationSEXP);
    Rcpp::traits::input_parameter< const bool >::type interweave(interweaveSEXP);
    Rcpp::traits::input_parameter< const double >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type expert_in(expert_inSEXP);
    rcpp_result_gen = Rcpp::wrap(stochvol::svsample_cpp(
          y_in,
          draws,
          burnin,
          X,
          priorspec_in,
          thinpara,
          thinlatent,
          keeptime_in,
          startpara,
          startlatent,
          keeptau,
          quiet,
          correct_model_specification,
          interweave,
          offset,
          expert_in));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _stochvol_svlsample_cpp(
    SEXP y_inSEXP,
    SEXP drawsSEXP,
    SEXP burninSEXP,
    SEXP XSEXP,
    SEXP priorspec_inSEXP,
    SEXP thinparaSEXP,
    SEXP thinlatentSEXP,
    SEXP keeptime_inSEXP,
    SEXP startparaSEXP,
    SEXP startlatentSEXP,
    SEXP keeptauSEXP,
    SEXP quietSEXP,
    SEXP correct_model_specificationSEXP,
    SEXP interweaveSEXP,
    SEXP offsetSEXP,
    SEXP expert_inSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y_in(y_inSEXP);
    Rcpp::traits::input_parameter< const int >::type draws(drawsSEXP);
    Rcpp::traits::input_parameter< const int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type priorspec_in(priorspec_inSEXP);
    Rcpp::traits::input_parameter< const int >::type thinpara(thinparaSEXP);
    Rcpp::traits::input_parameter< const int >::type thinlatent(thinlatentSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type keeptime_in(keeptime_inSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type startpara(startparaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type startlatent(startlatentSEXP);
    Rcpp::traits::input_parameter< const bool >::type keeptau(keeptauSEXP);
    Rcpp::traits::input_parameter< const bool >::type quiet(quietSEXP);
    Rcpp::traits::input_parameter< const bool >::type correct_model_specification(correct_model_specificationSEXP);
    Rcpp::traits::input_parameter< const bool >::type interweave(interweaveSEXP);
    Rcpp::traits::input_parameter< const double >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type expert_in(expert_inSEXP);
    rcpp_result_gen = Rcpp::wrap(stochvol::svlsample_cpp(
          y_in,
          draws,
          burnin,
          X,
          priorspec_in,
          thinpara,
          thinlatent,
          keeptime_in,
          startpara,
          startlatent,
          keeptau,
          quiet,
          correct_model_specification,
          interweave,
          offset,
          expert_in));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _stochvol_get_omori_constants () {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(stochvol::get_omori_constants());
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP run_testthat_tests ();

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _stochvol_Export_registerCCallable() { 
    R_RegisterCCallable("stochvol", "update_vanilla_sv", (DL_FUNC)stochvol::update_vanilla_sv);
    R_RegisterCCallable("stochvol", "update_general_sv", (DL_FUNC)stochvol::update_general_sv);
    R_RegisterCCallable("stochvol", "update_sv", (DL_FUNC)stochvol::update_sv);
    R_RegisterCCallable("stochvol", "update_svl", (DL_FUNC)stochvol::update_svl);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_stochvol_svsample_cpp", (DL_FUNC) &_stochvol_svsample_cpp, 16},
    {"_stochvol_svlsample_cpp", (DL_FUNC) &_stochvol_svlsample_cpp, 16},
    {"_stochvol_get_omori_constants", (DL_FUNC) &_stochvol_get_omori_constants, 0},
    {"run_testthat_tests", (DL_FUNC) &run_testthat_tests, 0},
    {"_stochvol_Export_registerCCallable", (DL_FUNC) &_stochvol_Export_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_stochvol(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

