#ifndef THETA_SAMPLER_H
#define THETA_SAMPLER_H

#include <RcppArmadillo.h>
#include "parameterization.hpp"

void draw_theta_rwMH(
    double& phi,
    double& rho,
    double& sigma2,
    double& mu,
    const arma::vec& y,
    const arma::vec& h,
    const arma::vec& ht,
    const arma::vec& prior_phi,
    const arma::vec& prior_rho,
    const arma::vec& prior_sigma2,
    const arma::vec& prior_mu,
    const Parameterization centering,
    const double stdev,
    const bool gammaprior);

#endif  // THETA_SAMPLER_H
