#ifndef UPDATE_FUNCTIONS_H
#define UPDATE_FUNCTIONS_H

// [[Rcpp::interfaces(cpp)]]

#include <Rcpp.h>

// a single MCMC update (normal errors):
// [[Rcpp::export]]
void update_sv(
    const Rcpp::NumericVector& data,
    Rcpp::NumericVector& curpara,
    Rcpp::NumericVector& h,
    double& h0,
    Rcpp::NumericVector& mixprob,
    Rcpp::IntegerVector& r,
    const bool centered_baseline,
    const double C0,
    const double cT,
    const double Bsigma,
    const double a0,
    const double b0,
    const double bmu,
    const double Bmu,
    const double B011inv,
    const double B022inv,
    const bool Gammaprior,
    const bool truncnormal,
    const double MHcontrol,
    const int MHsteps,
    const int parameterization,
    const bool dontupdatemu,
    const double priorlatent0);

// a single MCMC update (t errors):
void update_terr(
    const Rcpp::NumericVector&,
    Rcpp::NumericVector&,
    double &,
    const double,
    const double);

// [[Rcpp::export]]
void update_svl (
    const Rcpp::NumericVector& y,
    const Rcpp::NumericVector& y_star,
    const Rcpp::NumericVector& d,
    double& phi,
    double& rho,
    double& sigma2,
    double& mu,
    Rcpp::NumericVector& h,
    Rcpp::NumericVector& ht,
    const Rcpp::NumericVector& prior_phi,
    const Rcpp::NumericVector& prior_rho,
    const Rcpp::NumericVector& prior_sigma2,
    const Rcpp::NumericVector& prior_mu,
    const double stdev,
    const bool gammaprior,
    const Rcpp::IntegerVector& strategy);

#endif  // UPDATE_FUNCTIONS_H
