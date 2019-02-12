#ifndef THETA_UTILS_H
#define THETA_UTILS_H

#include <RcppArmadillo.h>
#include "parameterization.h"

double theta_log_likelihood(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const arma::vec& y,
    const arma::vec& h,
    const arma::vec& ht,
    const Parameterization centering);

double theta_log_prior(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const arma::vec& prior_phi,
    const arma::vec& prior_rho,
    const arma::vec& prior_sigma2,
    const arma::vec& prior_mu,
    const bool gammaprior);

arma::vec theta_transform(
    const double f,
    const double r,
    const double s,
    const double m);

arma::vec theta_transform_inv(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu);

double theta_transform_log_det_jac(
    const double f,
    const double r,
    const double s,
    const double m);

double theta_transform_inv_log_det_jac(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu);

arma::vec theta_propose(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const arma::vec& y,
    const arma::vec& h,
    const arma::vec& ht,
    const arma::mat& proposal_chol,
    const arma::mat& proposal_chol_inv);

double theta_log_likelihood_c(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const arma::vec& y,
    const arma::vec& h);

double theta_log_likelihood_nc(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const arma::vec& y,
    const arma::vec& h);

#endif  // THETA_UTILS_H
