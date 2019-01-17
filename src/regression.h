#ifndef REGRESSION_H
#define REGRESSION_H

#include <Rcpp.h>

// Step (b): sample mu, phi, sigma - __CENTERED__ version:
Rcpp::NumericVector regressionCentered(
       double h0, const Rcpp::NumericVector &h,
       double mu, double phi, double sigma,
       double C0, double cT, double Bsigma,
       double a0, double b0,
       double bmu, double Bmu,
       double B011inv, double B022inv,
       bool gammaprior, bool truncnormal,
       double MHcontrol, int MHsteps,
       const bool dontupdatemu, const double priorlatent0);

// Step (b): sample mu, phi, sigma - __NONCENTERED__ version:
Rcpp::NumericVector regressionNoncentered(
       const Rcpp::NumericVector &data,
       double h0, const Rcpp::NumericVector &h,
       const Rcpp::IntegerVector& r,
       double mu, double phi, double sigma,
       double Bsigma, double a0, double b0,
       double bmu, double Bmu,
       bool truncnormal, int MHsteps,
       const bool dontupdatemu, const double priorlatent0);

#endif  // REGRESSION_H

