#ifndef H_SAMPLER_H
#define H_SAMPLER_H

#include <Rcpp.h>

Rcpp::NumericVector draw_h_auxiliary(const Rcpp::NumericVector y_star,
                                     const Rcpp::NumericVector d,
                                     const Rcpp::NumericVector s,
                                     const double phi,
                                     const double rho,
                                     const double sigma2,
                                     const double mu,
                                     const Rcpp::CharacterVector centering,
                                     const Rcpp::DataFrame mixing_constants);

Rcpp::NumericVector draw_latent_auxiliaryMH(const Rcpp::NumericVector y,
                                            const Rcpp::NumericVector y_star,
                                            const Rcpp::NumericVector d,
                                            const Rcpp::NumericVector h,
                                            const double phi,
                                            const double rho,
                                            const double sigma2,
                                            const double mu,
                                            const Rcpp::DataFrame mixing_constants);

#endif  // H_SAMPLER_H
