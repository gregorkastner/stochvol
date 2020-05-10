#ifndef THETA_SAMPLER_H
#define THETA_SAMPLER_H

#include <RcppArmadillo.h>
#include <adaptation.hpp>
#include "parameterization.h"

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
    const stochvol::Adaptation<4>::Result& adaptation_proposal,
    const bool gammaprior,
    const Proposal sampler);

//void draw_thetamu_rwMH(
//    double& phi,
//    double& rho,
//    double& sigma2,
//    const double mu,
//    const arma::vec& y,
//    const arma::vec& h,
//    const arma::vec& ht,
//    const arma::vec& prior_phi,
//    const arma::vec& prior_rho,
//    const arma::vec& prior_sigma2,
//    const Parameterization centering,
//    const arma::mat& proposal_chol,
//    const arma::mat& proposal_chol_inv,
//    const bool gammaprior);

#endif  // THETA_SAMPLER_H
