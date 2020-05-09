#ifndef H_SAMPLER_H
#define H_SAMPLER_H

#include <RcppArmadillo.h>
#include "parameterization.h"

arma::vec draw_latent(
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
    const double priormu_sigma,
    const bool correct);

arma::vec draw_h_auxiliary(
    const arma::vec& y_star,
    const arma::ivec& d,
    const arma::uvec& z,
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const double priormu_mu,
    const double priormu_sigma,
    const Parameterization centering);

arma::vec correct_latent_auxiliaryMH(
    const arma::vec& y,
    const arma::vec& y_star,
    const arma::ivec& d,
    const arma::vec& h,
    const arma::vec& ht,
    const arma::vec& proposed,
    const double phi,
    const double rho,
    const double sigma2,
    const double mu);

#endif  // H_SAMPLER_H
