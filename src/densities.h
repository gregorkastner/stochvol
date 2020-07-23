#ifndef _DENSITIES_H_
#define _DENSITIES_H_

#include <Rmath.h>
#include <limits>

// Some (non-normalized) densities and transformations thereof
// (such as acceptance rates or other proportions)

namespace stochvol {

// non-normalized log-density for N(mu, sigma^2)
inline
double logdnorm(
    double x,
    double mu = 0,
    double sigma = 1) {
  return -std::log(sigma)-((x-mu)*(x-mu)/(2*sigma*sigma));
}

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
  return -(alpha - 1.) * std::log(x) - beta / x;
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
    const double sumtau,
    const double lambda,
    const int n) {
  using nl = std::numeric_limits<double>;
  static constexpr double negative_infinity = nl::has_infinity ? -nl::infinity() : nl::lowest();
  return nu > 2 ? .5 * nu * (n * std::log(.5 * nu) - sumtau) - n * std::lgamma(.5 * nu) - (nu - 2) * lambda : negative_infinity;
}

// first derivative of logdnu
inline
double dlogdnu(
    const double nu,
    const double sumtau,
    const double lambda,
    const int n) {
  return .5 * (n * (1 + std::log(.5 * nu) - ::Rf_digamma(.5 * nu)) - sumtau) - lambda;
}

// second derivative of logdnu
inline
double ddlogdnu(
    const double nu,
    const int n) {
  return .25 * n * (2 / nu - ::Rf_trigamma(.5 * nu));
}

}

#endif

