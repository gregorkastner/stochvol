#ifndef _SAMPLER_H_
#define _SAMPLER_H_

// Main sampling steps and helper functions

#include <Rcpp.h>

// Main sampler (as called from R):
SEXP sampler(const SEXP y_in, const SEXP draws_in,
  const SEXP burnin_in, const SEXP X_in,
  const SEXP bmu_in, const SEXP Bmu_in,
  const SEXP a0_in, const SEXP b0_in, const SEXP Bsigma_in,
  const SEXP thin_in, const SEXP timethin_in, const SEXP startpara_in,
  const SEXP startvol_in, const SEXP keeptau_in,
  const SEXP quiet_in, const SEXP para_in,
  const SEXP MHsteps_in, const SEXP B011_in, const SEXP B022_in,
  const SEXP mhcontrol_in, const SEXP gammaprior_in,
  const SEXP truncnormal_in, const SEXP offset_in,
  const SEXP dontupdatemu_in, const SEXP priordf_in,
  const SEXP priorbeta_in, const SEXP priorlatent0_in);

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
       const int * const r,
       double mu, double phi, double sigma,
       double Bsigma, double a0, double b0,
       double bmu, double Bmu,
       bool truncnormal, int MHsteps,
       const bool dontupdatemu, const double priorlatent0);

// a single MCMC update (normal errors):
void update(const Rcpp::NumericVector &, double *, double *, double &,
            double *, int *, const bool, const double,
	    const double, const double, const double, const double,
	    const double, const double, const double, const double,
	    const bool, const bool, const double, const int, const int,
	    const bool, const double);

// a single MCMC update (t errors):
void update_terr(const Rcpp::NumericVector &,
                 double *, double &, const double, const double);

#endif
