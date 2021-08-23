/*
 * R package stochvol by
 *     Gregor Kastner Copyright (C) 2013-2018
 *     Gregor Kastner and Darjus Hosszejni Copyright (C) 2019-
 *  
 *  This file is part of the R package stochvol: Efficient Bayesian
 *  Inference for Stochastic Volatility Models.
 *  
 *  The R package stochvol is free software: you can redistribute it
 *  and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation, either version 2 or
 *  any later version of the License.
 *  
 *  The R package stochvol is distributed in the hope that it will be
 *  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 *  General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with the R package stochvol. If that is not the case, please
 *  refer to <http://www.gnu.org/licenses/>.
 */

/*
 * densities.h
 * 
 * Inline definition of some non-normalized densities and transformations
 * thereof such as acceptance rates or other proportions.
 * These implementations are more efficient than their Rmath.h counterparts
 * because there is no validation, no normalization, and some expressions can
 * be pre-computed for these functions.
 */

#ifndef _DENSITIES_H_
#define _DENSITIES_H_

#include <Rmath.h>
#include <limits>

namespace stochvol {

// non-normalized log-density for N(mu, sigma^2)
// with log_sd as input
inline
double logdnorm2(
    const double x,
    const double mu,
    const double sigma,
    const double log_sigma = 0) {
  const double z = (x - mu) / sigma;
  return -.5 * z * z - log_sigma;
}

// non-normalized log-density for N(mu, sigma^2)
inline
double logdnorm(
    double x,
    double mu = 0,
    double sigma = 1) {
  return logdnorm2(x, mu, sigma, std::log(sigma));
}

// non-normalized log-density of the gamma distribution
inline
double logdgamma(
    const double x,
    const double alpha,
    const double beta) {
  return (alpha - 1.) * std::log(x) - beta * x;
}

// non-normalized log-density of the inverse gamma distribution
inline
double logdinvgamma(
    const double x,
    const double alpha,
    const double beta) {
  return (-alpha - 1.) * std::log(x) - beta / x;
}

// non-normalized log-density for Beta(a, b)
inline
double logdbeta(
    const double x,
    const double a,
    const double b) {
  return (a - 1) * std::log(x) + (b - 1) * std::log(1 - x);
}

// acceptance ratio for prior matching when sampling sigma
// Proposal is InvGamma(-0.5, 0)
// Target is Gamma(.5, 1/(2*Bsigma))
inline
double logacceptrateGamma(
    const double xnew,
    const double xold,
    const double Bsigma) {
  return (xold - xnew) / (2 * Bsigma);
}

// acceptance ratio for log normal random walk
inline
double logacceptrateRW(
    const double xnew,
    const double xold,
    const double Bsigma,
    const int T,
    const double z) {
  return .5 * (T * std::log(xold / xnew) + (xold - xnew) / Bsigma + (1 / xold - 1 / xnew) * z);
}

// proportion of two beta-distributions with same parameters
// (evaluated at two different points)
inline
double propBeta(
    const double x,
    const double y,
    const double a,
    const double b) {
  return std::pow(x / y, a - 1) * std::pow((1 - x) / (1 - y), b - 1);
}

// full conditional non-normalized posterior log-density of the
// degrees of freedom parameter nu
inline
double logdnu(
    const double nu,
    const double sum_tau,
    const double lambda,
    const int n) {
  using nl = std::numeric_limits<double>;
  static constexpr double negative_infinity = nl::has_infinity ? -nl::infinity() : nl::lowest();
  return nu > 2 ?
    .5 * nu * (-sum_tau + n * std::log(.5*(nu-2))) - n*std::lgamma(.5*nu) - (nu-2)*lambda:
    negative_infinity;
}

// first derivative of logdnu
inline
double dlogdnu(
    const double nu,
    const double sum_tau,
    const double lambda,
    const int n) {
  return .5 * (n*(nu/(nu-2) + std::log(.5*(nu-2)) - ::Rf_digamma(.5*nu)) - sum_tau) - lambda;
}

// second derivative of logdnu
inline
double ddlogdnu(
    const double nu,
    const int n) {
  return .25 * n * (2*(nu-4) * std::pow(nu-2, -2) - ::Rf_trigamma(.5*nu));
}

}

#endif

