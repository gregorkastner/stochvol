#ifndef _DENSITIES_H_
#define _DENSITIES_H_

#include <Rmath.h>

// Some (non-normalized) densities and transformations thereof
// (such as acceptance rates or other proportions)

// non-normalized log-density for N(mu, sigma^2)
inline double logdnorm(
    double x,
    double mu = 0,
    double sigma = 1) {
 return -log(sigma)-((x-mu)*(x-mu)/(2*sigma*sigma));
}

// non-normalized log-density for Beta(a, b)
inline double logdbeta(
    double x,
    double a,
    double b) {
 return (a-1)*log(x)+(b-1)*log(1-x);
}

// acceptance ratio for prior matching when sampling sigma
// Proposal is InvGamma(-0.5, 0)
// Target is Gamma(.5, 1/(2*Bsigma))
inline double logacceptrateGamma(
    double xnew,
    double xold,
    double Bsigma) {
 return (xold-xnew)/(2*Bsigma);
}

// acceptance ratio for log normal random walk
inline double logacceptrateRW(
    double xnew,
    double xold,
    double Bsigma,
    int T,
    double z) {
 return .5*(T*log(xold/xnew)+(xold-xnew)/Bsigma+(1/xold-1/xnew)*z);
}

// proportion of two beta-distributions with same parameters
// (evaluated at two different points)
inline double propBeta(
    double x,
    double y,
    double a,
    double b) {
 return pow(x/y, a-1)*pow((1-x)/(1-y), b-1);
}

// full conditional non-normalized posterior log-density of the
// degrees of freedom parameter nu
inline double logdnu(
    double nu,
    double sumtau,
    int n) {
 return .5 * nu * (n*log(.5*nu) - sumtau) - n*lgamma(.5*nu);
}

// first derivative of logdnu
inline double dlogdnu(
    double nu,
    double sumtau,
    int n) {
 return .5 * (n * ( 1 + log(.5*nu) - Rf_digamma(.5*nu)) - sumtau);
}

// second derivative of logdnu
inline double ddlogdnu(
    double nu,
    int n) {
 return .25 * n * (2/nu - Rf_trigamma(.5*nu));
}


#endif
