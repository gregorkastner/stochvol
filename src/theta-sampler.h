#ifndef THETA_SAMPLER_H
#define THETA_SAMPLER_H

#include <Rcpp.h>
#include "parameterization.hpp"

void draw_theta_rwMH(
    double& phi,
    double& rho,
    double& sigma2,
    double& mu,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericVector& h,
    const Rcpp::NumericVector& ht,
    const Rcpp::NumericVector& prior_phi,
    const Rcpp::NumericVector& prior_rho,
    const Rcpp::NumericVector& prior_sigma2,
    const Rcpp::NumericVector& prior_mu,
    const Parameterization centering,
    const double stdev,
    const bool gammaprior);

#endif  // THETA_SAMPLER_H
