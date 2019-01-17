#ifndef _SAMPLER_H_
#define _SAMPLER_H_

// [[Rcpp::interfaces(r)]]

#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::List svsample_cpp(
    const Rcpp::NumericVector& y_in,
    const int draws,
    const int burnin,
    const Rcpp::NumericMatrix& X_in,
    const double bmu,
    const double Bmu,
    const double a0,
    const double b0,
    const double Bsigma,
    const int thin,
    const int timethin,
    const Rcpp::List& startpara_in,
    const Rcpp::NumericVector& startvol_in,
    const bool keeptau,
    const bool quiet,
    const int para,
    const int MHsteps,
    const double B011,
    const double B022,
    const double mhcontrol,
    const bool gammaprior,
    const bool truncnormal,
    const double offset,
    const bool dontupdatemu,
    const Rcpp::NumericVector& priordf_in,
    const Rcpp::NumericVector& priorbeta_in,
    const double priorlatent0);

// [[Rcpp::export]]
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

#endif  // _SAMPLER_H_

