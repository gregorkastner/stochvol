/*
 * R package stochvol by
 *     Gregor Kastner Copyright (C) 2013-2021
 *     Darjus Hosszejni Copyright (C) 2019-2021
 *
 *  This file is part of the R package stochvol: Efficient Bayesian
 *  Inference for Stochastic Volatility Models.
 *
 *  The R package stochvol is free software: you can redistribute it
 *  and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation, either version 2 or
 *  any later version of the License.
 *
 *  The R package stochvol is distributed in the hope that it will be
 *  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with the R package stochvol. If that is not the case, please
 *  refer to <http://www.gnu.org/licenses/>.
 */

/*
 * exports.cc
 *
 * Transform Rcpp implementations to .Call()-conform ones.
 */

#include <RcppArmadillo.h>
#include "single_update.h"
#include "sampling_main.h"
#include "utils_latent_states.h"
#include "utils_main.h"

using namespace Rcpp;

// svsample_fast_cpp
RcppExport SEXP _stochvol_svsample_fast_cpp(
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
    SEXP print_settingsSEXP,
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
    Rcpp::traits::input_parameter< const Rcpp::List >::type print_settings(print_settingsSEXP);
    Rcpp::traits::input_parameter< const bool >::type correct_model_specification(correct_model_specificationSEXP);
    Rcpp::traits::input_parameter< const bool >::type interweave(interweaveSEXP);
    Rcpp::traits::input_parameter< const double >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type expert_in(expert_inSEXP);
    rcpp_result_gen = Rcpp::wrap(stochvol::svsample_fast_cpp(
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
          print_settings,
          correct_model_specification,
          interweave,
          offset,
          expert_in));
    return rcpp_result_gen;
END_RCPP
}

// svsample_general_cpp
RcppExport SEXP _stochvol_svsample_general_cpp(
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
    SEXP print_settingsSEXP,
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
    Rcpp::traits::input_parameter< const Rcpp::List >::type print_settings(print_settingsSEXP);
    Rcpp::traits::input_parameter< const bool >::type correct_model_specification(correct_model_specificationSEXP);
    Rcpp::traits::input_parameter< const bool >::type interweave(interweaveSEXP);
    Rcpp::traits::input_parameter< const double >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type expert_in(expert_inSEXP);
    rcpp_result_gen = Rcpp::wrap(stochvol::svsample_general_cpp(
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
          print_settings,
          correct_model_specification,
          interweave,
          offset,
          expert_in));
    return rcpp_result_gen;
END_RCPP
}

// get_omori_constants
RcppExport SEXP _stochvol_get_omori_constants () {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(stochvol::get_omori_constants());
    return rcpp_result_gen;
END_RCPP
}

// Declare testthat function for the C++ unit tests
RcppExport SEXP run_testthat_tests ();

// Register entry points for exported C++ functions
RcppExport SEXP _stochvol_Export_registerCCallable() {
    R_RegisterCCallable("stochvol", "update_fast_sv", (DL_FUNC)stochvol::update_fast_sv);
    R_RegisterCCallable("stochvol", "update_t_error", (DL_FUNC)stochvol::update_t_error);
    R_RegisterCCallable("stochvol", "update_general_sv", (DL_FUNC)stochvol::update_general_sv);
    R_RegisterCCallable("stochvol", "update_regressors", (DL_FUNC)stochvol::update_regressors);
    R_RegisterCCallable("stochvol", "list_to_priorspec", (DL_FUNC)stochvol::list_to_priorspec);
    R_RegisterCCallable("stochvol", "list_to_general_sv", (DL_FUNC)stochvol::list_to_general_sv);
    R_RegisterCCallable("stochvol", "list_to_fast_sv", (DL_FUNC)stochvol::list_to_fast_sv);
    R_RegisterCCallable("stochvol", "update_sv", (DL_FUNC)stochvol::update_sv);
    return R_NilValue;
}

// List of exposed functions in the shared library
static const R_CallMethodDef CallEntries[] = {
    {"_stochvol_svsample_fast_cpp", (DL_FUNC) &_stochvol_svsample_fast_cpp, 16},
    {"_stochvol_svsample_general_cpp", (DL_FUNC) &_stochvol_svsample_general_cpp, 16},
    {"_stochvol_get_omori_constants", (DL_FUNC) &_stochvol_get_omori_constants, 0},
    {"run_testthat_tests", (DL_FUNC) &run_testthat_tests, 0},
    {"_stochvol_Export_registerCCallable", (DL_FUNC) &_stochvol_Export_registerCCallable, 0},
    {NULL, NULL, 0}
};

// Gets executed when loading the shared library
RcppExport void R_init_stochvol(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

