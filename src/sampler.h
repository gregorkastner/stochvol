#ifndef _SAMPLER_H_
#define _SAMPLER_H_

// Main sampling steps and helper functions

#include <Rcpp.h>

// Main sampler (as called from R):
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
