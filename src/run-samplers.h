#ifndef RUN_SAMPLERS_H
#define RUN_SAMPLERS_H

#include <Rcpp.h>

Rcpp::List svlsample_cpp (
    const int draws,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericVector& y_star,
    const Rcpp::NumericVector& d,
    const int burnin,
    const int thinpara,
    const int thinlatent,
    const int thintime,
    const double phi_init,
    const double rho_init,
    const double sigma2_init,
    const double mu_init,
    const Rcpp::NumericVector& h_init,
    const double prior_phi_a,
    const double prior_phi_b,
    const double prior_rho_a,
    const double prior_rho_b,
    const double prior_sigma2_shape,
    const double prior_sigma2_rate,
    const double prior_mu_mu,
    const double prior_mu_sigma,
    const bool verbose,
    const double stdev,
    const bool gammaprior,
    const Rcpp::CharacterVector& strategy);

#endif  // RUN_SAMPLERS_H

