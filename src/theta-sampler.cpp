#include <Rcpp.h>
#include <cmath>
#include "theta-sampler.h"
#include "theta-utils.h"
#include "parameterization.hpp"

using namespace Rcpp;

void draw_theta_rwMH(
    double& phi,
    double& rho,
    double& sigma2,
    double& mu,
    const NumericVector& y,
    const NumericVector& h,
    const NumericVector& ht,
    const NumericVector& prior_phi,
    const NumericVector& prior_rho,
    const NumericVector& prior_sigma2,
    const NumericVector& prior_mu,
    const Parameterization centering,
    const double stdev,
    const bool gammaprior) {
  const NumericVector proposed = theta_propose(phi, rho, sigma2, mu, y, h, ht, stdev);
  const double phi_prop = proposed[0], rho_prop = proposed[1], sigma2_prop = proposed[2],
    mu_prop = proposed[3], prop_old_logdens = proposed[4], prop_new_logdens = proposed[5];
  const double log_acceptance = (theta_log_prior(phi_prop, rho_prop, sigma2_prop, mu_prop, prior_phi, prior_rho, prior_sigma2, prior_mu, gammaprior) +
    theta_log_likelihood(phi_prop, rho_prop, sigma2_prop, mu_prop, y, h, ht, centering)) -
    (theta_log_prior(phi, rho, sigma2, mu, prior_phi, prior_rho, prior_sigma2, prior_mu, gammaprior) +
    theta_log_likelihood(phi, rho, sigma2, mu, y, h, ht, centering)) -
    (prop_new_logdens - prop_old_logdens);
  if (log_acceptance > 0 || exp(log_acceptance) > R::runif(0, 1)) {
    phi = phi_prop;
    rho = rho_prop;
    sigma2 = sigma2_prop;
    mu = mu_prop;
  }
}

