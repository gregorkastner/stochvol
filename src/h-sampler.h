#ifndef H_SAMPLER_H
#define H_SAMPLER_H

#include <Rcpp.h>
#include "parameterization.hpp"

Rcpp::NumericVector draw_h_auxiliary(const Rcpp::NumericVector y_star,
                                     const Rcpp::NumericVector d,
                                     const Rcpp::NumericVector s,
                                     const double phi,
                                     const double rho,
                                     const double sigma2,
                                     const double mu,
                                     const double priormu_mu,
                                     const double priormu_sigma,
                                     const Parameterization centering);

Rcpp::NumericVector draw_latent_auxiliaryMH(const Rcpp::NumericVector y,
                                            const Rcpp::NumericVector y_star,
                                            const Rcpp::NumericVector d,
                                            const Rcpp::NumericVector h,
                                            const double phi,
                                            const double rho,
                                            const double sigma2,
                                            const double mu,
                                            const double priormu_mu,
                                            const double priormu_sigma);

#endif  // H_SAMPLER_H
