#ifndef _SINGLE_UPDATES_H_
#define _SINGLE_UPDATES_H_

#include <RcppArmadillo.h>
#include <adaptation.hpp>

// a single MCMC update (normal errors):
void update_sv(
    const arma::vec& data,
    arma::vec& curpara,
    arma::vec& h,
    double& h0,
    arma::vec& mixprob,
    arma::ivec& r,
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
    const arma::vec& data,
    arma::vec& tau,
    double& nu,
    const double lambda);

void update_svl(
    const arma::vec& y,
    const arma::vec& y_star,
    const arma::ivec& d,
    double& phi,
    double& rho,
    double& sigma2,
    double& mu,
    double& h0,
    arma::vec& h,
    arma::vec& ht,
    stochvol::AdaptationCollection& adaptation,
    const arma::vec& prior_phi,
    const arma::vec& prior_rho,
    const arma::vec& prior_sigma2,
    const arma::vec& prior_mu,
    const bool use_mala,
    const bool gammaprior,
    const bool correct,
    const arma::ivec& strategy,
    const bool dontupdatemu);

#endif  // _SINGLE_UPDATES_H_
