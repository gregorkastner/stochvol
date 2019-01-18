#ifndef REGRESSION_H
#define REGRESSION_H

#include <RcppArmadillo.h>

// Step (b): sample mu, phi, sigma - __CENTERED__ version:
arma::vec regressionCentered(
    double h0,
    const arma::vec &h,
    double mu,
    double phi,
    double sigma,
    double C0,
    double cT,
    double Bsigma,
    double a0,
    double b0,
    double bmu,
    double Bmu,
    double B011inv,
    double B022inv,
    bool gammaprior,
    bool truncnormal,
    double MHcontrol,
    int MHsteps,
    const bool dontupdatemu,
    const double priorlatent0);

// Step (b): sample mu, phi, sigma - __NONCENTERED__ version:
arma::vec regressionNoncentered(
    const arma::vec& data,
    double h0,
    const arma::vec& h,
    const arma::ivec& r,
    double mu,
    double phi,
    double sigma,
    double Bsigma,
    double a0,
    double b0,
    double bmu,
    double Bmu,
    bool truncnormal,
    int MHsteps,
    const bool dontupdatemu,
    const double priorlatent0);

#endif  // REGRESSION_H

