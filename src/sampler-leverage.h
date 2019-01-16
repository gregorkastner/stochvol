#ifndef RUN_SAMPLERS_H
#define RUN_SAMPLERS_H

#include <RcppArmadillo.h>

Rcpp::List svlsample_cpp (
    const Rcpp::NumericVector& y,
    const int draws,
    const int burnin,
    const Rcpp::NumericMatrix& X,
    const int thinpara,
    const int thinlatent,
    const int thintime,
    const Rcpp::List& theta_init,
    const Rcpp::NumericVector& h_init,
    const double prior_phi_a,
    const double prior_phi_b,
    const double prior_rho_a,
    const double prior_rho_b,
    const double prior_sigma2_shape,
    const double prior_sigma2_rate,
    const double prior_mu_mu,
    const double prior_mu_sigma,
    const double prior_beta_mu,
    const double prior_beta_sigma,
    const bool verbose,
    const double offset,
    const double stdev,
    const bool gammaprior,
    const Rcpp::CharacterVector& strategy);

void update_svl (
    const Rcpp::NumericVector& y,
    const Rcpp::NumericVector& y_star,
    const Rcpp::NumericVector& d,
    Rcpp::NumericVector& theta,
    Rcpp::NumericVector& h,
    Rcpp::NumericVector& ht,
    const double prior_phi_a,
    const double prior_phi_b,
    const double prior_rho_a,
    const double prior_rho_b,
    const double prior_sigma2_shape,
    const double prior_sigma2_rate,
    const double prior_mu_mu,
    const double prior_mu_sigma,
    const double stdev,
    const bool gammaprior,
    const Rcpp::IntegerVector& strategy);

#endif  // RUN_SAMPLERS_H

