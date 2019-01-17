#ifndef AUG_KALMAN_FILTER_H
#define AUG_KALMAN_FILTER_H

#include <Rcpp.h>
#include "parameterization.hpp"

// Augmented Kalman filter
// 
// Implementation of the augmented Kalman filter for the special
// case of the current model. The augmented Kalman filter is a
// modified version of the standard Kalman filter that
// considers the the error terms' correlation and the constant
// \eqn{\mu} in the state equation.
// Used to calculate the likelihood of \eqn{\theta=(\phi,\sigma,\rho)},
// and to calculate the posterior of \eqn{\mu}.
// For a specification see Appendix B of Nakajima, 2009.
// @author hdarjus \email{hdarjus@gmail.com}
// @param dPhi parameter \eqn{\phi}
// @param dSigma2 square of parameter \eqn{\sigma}
// @param dRho parameter \eqn{\rho}
// @param va numeric vector of mixture states constants
// @param vb numeric vector of mixture states constants
// @param vm numeric vector of mixture states constants
// @param vv numeric vector of mixture states constants
// @param vD vector of signs
// @param vYStar vector of y*
// @param dMuMu prior mean of \eqn{\mu}
// @param dSigma2Mu prior variance of \eqn{\mu}
// @return List with elements \code{D, J1, L, f, F, hts, v, Q, q, jt22, h1var}, from which
//   \code{D, J1, L, f, F, h1var, v} are numeric vectors of length T,
//   and \code{Q, q, jt22, h1var} are numbers. All of them are just partial results helping
//   later calculations.
// @references Nakajima, Jouchi, and Yasuhiro Omori.
//   "Leverage, heavy-tails and correlated jumps in stochastic volatility models."
//   Computational Statistics & Data Analysis 53.6 (2009): 2335-2353.
Rcpp::List aug_kalman_filter(
    const double phi,
    const double rho,
    const double sigma2,
    const Rcpp::NumericVector& a,
    const Rcpp::NumericVector& b,
    const Rcpp::NumericVector& m,
    const Rcpp::NumericVector& v,
    const Rcpp::NumericVector& d,
    const Rcpp::NumericVector& y_star,
    const double mu_mu,
    const double sigma2_mu,
    const Parameterization centering);

Rcpp::List aug_kalman_filter_c(
    const double phi,
    const double rho,
    const double sigma2,
    const Rcpp::NumericVector& a,
    const Rcpp::NumericVector& b,
    const Rcpp::NumericVector& m,
    const Rcpp::NumericVector& v,
    const Rcpp::NumericVector& d,
    const Rcpp::NumericVector& y_star,
    const double mu_mu,
    const double sigma2_mu);

Rcpp::List aug_kalman_filter_nc(
    const double phi,
    const double rho,
    const double sigma2,
    const Rcpp::NumericVector& a,
    const Rcpp::NumericVector& b,
    const Rcpp::NumericVector& m,
    const Rcpp::NumericVector& v,
    const Rcpp::NumericVector& d,
    const Rcpp::NumericVector& y_star,
    const double mu_mu,
    const double sigma2_mu);

#endif  // AUG_KALMAN_FILTER_H
