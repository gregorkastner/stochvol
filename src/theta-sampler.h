#ifndef THETA_SAMPLER_H
#define THETA_SAMPLER_H

#include <Rcpp.h>
#include "parameterization.hpp"

Rcpp::NumericVector draw_theta_rwMH(const double phi, const double rho,
                                    const double sigma2, const double mu,
                                    const Rcpp::NumericVector y,
                                    const Rcpp::NumericVector h,
                                    const Rcpp::NumericVector prior_phi,
                                    const Rcpp::NumericVector prior_rho,
                                    const Rcpp::NumericVector prior_sigma2,
                                    const Rcpp::NumericVector prior_mu,
                                    const Parameterization centering,
                                    const double stdev,
                                    const bool gammaprior);

#endif  // THETA_SAMPLER_H
