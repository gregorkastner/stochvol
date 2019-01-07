#ifndef MIXTURE_STATE_SAMPLER_H
#define MIXTURE_STATE_SAMPLER_H

#include <Rcpp.h>

Rcpp::NumericVector draw_s_auxiliary(const Rcpp::NumericVector y_star,
                                     const Rcpp::NumericVector d,
                                     const Rcpp::NumericVector h,
                                     const double phi,
                                     const double rho,
                                     const double sigma2,
                                     const double mu,
                                     const Rcpp::CharacterVector centering);

Rcpp::NumericMatrix mixture_state_post_dist(const Rcpp::NumericVector eps_star,
                                            const Rcpp::NumericVector eta,
                                            const Rcpp::NumericVector d,
                                            const double mu,
                                            const double sigma2,
                                            const double rho,
                                            const Rcpp::CharacterVector centering);

#endif  // MIXTURE_STATE_SAMPLER_H

