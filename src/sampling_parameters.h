#ifndef _SAMPLING_PARAMETERS_H_
#define _SAMPLING_PARAMETERS_H_

// Sampling the time independent model parameters

#include <RcppArmadillo.h>
#include <adaptation.hpp>
#include "type_definitions.h"

namespace stochvol {

struct ReturnRegression {
  double mu;
  double phi;
  double sigma;
};

// Step (b): sample mu, phi, sigma - __CENTERED__ version:
ReturnRegression regression_centered(
    const double h0,
    const arma::vec &h,
    double mu,
    double phi,
    double sigma,
    const double C0,
    const double cT,
    const double Bsigma,
    const double a0,
    const double b0,
    const double bmu,
    const double Bmu,
    const double B011inv,
    const double B022inv,
    const bool gammaprior,
    const bool truncnormal,
    const double MHcontrol,
    const int MHsteps,
    const bool dontupdatemu,
    const double priorlatent0);

// Step (b): sample mu, phi, sigma - __NONCENTERED__ version:
ReturnRegression regression_noncentered(
    const arma::vec& data,
    const double h0,
    const arma::vec& h,
    const arma::ivec& r,
    double mu,
    double phi,
    double sigma,
    const double Bsigma,
    const double a0,
    const double b0,
    const double bmu,
    const double Bmu,
    const bool truncnormal,
    const int MHsteps,
    const bool dontupdatemu,
    const double priorlatent0);

bool draw_theta(
    double& phi,
    double& rho,
    double& sigma2,
    double& mu,
    const arma::vec& y,
    const double h0,
    const arma::vec& h,
    const arma::vec& ht,
    const arma::vec& exp_h_half,
    arma::vec& exp_h_half_proposal_nc,
    const arma::vec& prior_phi,
    const arma::vec& prior_rho,
    const arma::vec& prior_sigma2,
    const arma::vec& prior_mu,
    const Parameterization centering,
    const ProposalDiffusionKen& adaptation_proposal,
    const bool gammaprior,
    const Proposal sampler);

bool draw_thetamu_rwMH(
    double& phi,
    double& rho,
    double& sigma2,
    const double mu,
    const arma::vec& y,
    const double h0,
    const arma::vec& h,
    const arma::vec& ht,
    const arma::vec& exp_h_half,
    arma::vec& exp_h_half_proposal_nc,
    const arma::vec& prior_phi,
    const arma::vec& prior_rho,
    const arma::vec& prior_sigma2,
    const Parameterization centering,
    const ProposalDiffusionKen& adaptation_proposal,
    const bool gammaprior);

}

#endif  // THETA_SAMPLER_H

