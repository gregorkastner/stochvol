#ifndef _UTILS_LATENT_STATES_H_
#define _UTILS_LATENT_STATES_H_

#include <RcppArmadillo.h>

double h_log_posterior(
    const arma::vec& h,
    const arma::vec& y,
    const double phi,
    const double rho,
    const double sigma2,
    const double mu);

double h_aux_log_posterior(
    const arma::vec& h,
    const arma::vec& y_star,
    const arma::ivec& d,
    const double phi,
    const double rho,
    const double sigma2,
    const double mu);

#endif  // H_UTILS_H
