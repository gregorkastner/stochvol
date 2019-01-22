#ifndef H_SAMPLER_H
#define H_SAMPLER_H

#include <RcppArmadillo.h>
#include "parameterization.hpp"

arma::vec draw_h_auxiliary(
    const arma::vec& y_star,
    const arma::ivec& d,
    const arma::vec& s,
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const double priormu_mu,
    const double priormu_sigma,
    const Parameterization centering);

arma::vec draw_latent_auxiliaryMH(
    const arma::vec& y,
    const arma::vec& y_star,
    const arma::ivec& d,
    const arma::vec& h,
    const arma::vec& ht,
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const double priormu_mu,
    const double priormu_sigma);

#endif  // H_SAMPLER_H
