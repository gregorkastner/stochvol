#ifndef H_UTILS_H
#define H_UTILS_H

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
