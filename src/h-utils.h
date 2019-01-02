#ifndef H_UTILS_H
#define H_UTILS_H

#include <Rcpp.h>

double h_log_posterior(const Rcpp::NumericVector h,
                       const Rcpp::NumericVector y,
                       const double phi,
                       const double rho,
                       const double sigma2,
                       const double mu);

double h_aux_log_posterior(const Rcpp::NumericVector h,
                           const Rcpp::NumericVector y_star,
                           const Rcpp::NumericVector d,
                           const Rcpp::NumericVector mixing_a,
                           const Rcpp::NumericVector mixing_b,
                           const Rcpp::NumericVector mixing_m,
                           const Rcpp::NumericVector mixing_p,
                           const Rcpp::NumericVector mixing_v,
                           const double phi,
                           const double rho,
                           const double sigma2,
                           const double mu);

#endif  // H_UTILS_H
