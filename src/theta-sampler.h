#ifndef THETA_SAMPLER_H
#define THETA_SAMPLER_H

#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericVector draw_theta_rwMH(const double phi, const double rho,
                                    const double sigma2, const double mu,
                                    const Rcpp::NumericVector y,
                                    const Rcpp::NumericVector h,
                                    const Rcpp::NumericVector prior_phi,
                                    const Rcpp::NumericVector prior_rho,
                                    const Rcpp::NumericVector prior_sigma2,
                                    const Rcpp::NumericVector prior_mu,
                                    const Rcpp::CharacterVector centering,
                                    const double stdev = .1);

// [[Rcpp::export]]
Rcpp::NumericVector draw_theta_auxiliary(const double phi, const double rho,
                                         const double sigma2, const double mu,
                                         const Rcpp::NumericVector y_star,
                                         const Rcpp::NumericVector d,
                                         const Rcpp::NumericVector s,
                                         const Rcpp::NumericVector prior_phi,
                                         const Rcpp::NumericVector prior_rho,
                                         const Rcpp::NumericVector prior_sigma2,
                                         const Rcpp::NumericVector prior_mu,
                                         const Rcpp::DataFrame mixing_constants);

#endif  // THETA_SAMPLER_H
