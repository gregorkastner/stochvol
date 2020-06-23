#ifndef _SAMPLING_PARAMETERS_H_
#define _SAMPLING_PARAMETERS_H_

// Sampling the time independent model parameters

#include <RcppArmadillo.h>
#include <adaptation.hpp>
#include "type_definitions.h"

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

bool draw_theta(
    double& phi,
    double& rho,
    double& sigma2,
    double& mu,
    const arma::vec& y,
    const arma::vec& h,
    const arma::vec& ht,
    const arma::vec& exp_h_half,
    arma::vec& exp_h_half_proposal_nc,
    const arma::vec& prior_phi,
    const arma::vec& prior_rho,
    const arma::vec& prior_sigma2,
    const arma::vec& prior_mu,
    const Parameterization centering,
    const stochvol::Adaptation::Result& adaptation_proposal,
    const bool gammaprior,
    const Proposal sampler);

bool draw_thetamu_rwMH(
    double& phi,
    double& rho,
    double& sigma2,
    const double mu,
    const arma::vec& y,
    const arma::vec& h,
    const arma::vec& ht,
    const arma::vec& exp_h_half,
    arma::vec& exp_h_half_proposal_nc,
    const arma::vec& prior_phi,
    const arma::vec& prior_rho,
    const arma::vec& prior_sigma2,
    const Parameterization centering,
    const stochvol::Adaptation::Result& adaptation_proposal,
    const bool gammaprior);

#endif  // THETA_SAMPLER_H

