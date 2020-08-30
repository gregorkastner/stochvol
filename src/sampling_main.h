#ifndef _SAMPLING_MAIN_H_
#define _SAMPLING_MAIN_H_

#include <RcppArmadillo.h>

namespace stochvol {

Rcpp::List svsample_cpp(
    const arma::vec& y_in,
    const int draws,
    const int burnin,
    const arma::mat& X_in,
    const Rcpp::List& priorspec_in,
    const int thinpara,
    const int thinlatent,
    const Rcpp::CharacterVector& keeptime_in,
    const Rcpp::List& startpara,
    const arma::vec& startlatent,
    const bool keeptau,
    const bool quiet,
    const bool correct_model_specification,
    const bool interweave,
    const double offset,
    const Rcpp::List& expert);

Rcpp::List svlsample_cpp (
    const arma::vec& y_in,
    const int draws,
    const int burnin,
    const arma::mat& X_in,
    const Rcpp::List& priorspec_in,
    const int thinpara,
    const int thinlatent,
    const Rcpp::CharacterVector& keeptime_in,
    const Rcpp::List& startpara,
    const arma::vec& startlatent,
    const bool keeptau,
    const bool quiet,
    const bool correct_model_specification,
    const bool interweave,
    const double offset,
    const Rcpp::List& expert);

}

#endif  // _SAMPLING_MAIN_H_

