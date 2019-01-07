#include <Rcpp.h>
#include <string>
#include <cmath>
#include "theta-sampler.h"
#include "theta-utils.h"
using namespace Rcpp;

NumericVector draw_theta_rwMH(const double phi, const double rho,
                              const double sigma2, const double mu,
                              const NumericVector y, const NumericVector h,
                              const NumericVector prior_phi,
                              const NumericVector prior_rho,
                              const NumericVector prior_sigma2,
                              const NumericVector prior_mu,
                              const CharacterVector centering,
                              const double stdev,
                              const bool gammaprior) {
  NumericVector result;
  std::string scentering = as<std::string>(centering);
  const NumericVector proposed = theta_propose(phi, rho, sigma2, mu, y, h, stdev);
  const double phi_prop = proposed(0), rho_prop = proposed(1), sigma2_prop = proposed(2),
    mu_prop = proposed(3), prop_old_logdens = proposed(4), prop_new_logdens = proposed(5);
  const double log_acceptance = (theta_log_prior(phi_prop, rho_prop, sigma2_prop, mu_prop, prior_phi, prior_rho, prior_sigma2, prior_mu, gammaprior) +
    theta_log_likelihood(phi_prop, rho_prop, sigma2_prop, mu_prop, y, h, centering)) -
    (theta_log_prior(phi, rho, sigma2, mu, prior_phi, prior_rho, prior_sigma2, prior_mu, gammaprior) +
    theta_log_likelihood(phi, rho, sigma2, mu, y, h, centering)) -
    (prop_new_logdens - prop_old_logdens);
  if (log_acceptance > 0 || exp(log_acceptance) > R::runif(0, 1)) {
    result = NumericVector::create(phi_prop, rho_prop, sigma2_prop, mu_prop);
  } else {
    result = NumericVector::create(phi, rho, sigma2, mu);
  }
  return result;
}

