#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <iterator>
#include <numeric>
#include <string>
#include "aug-kalman-filter.h"  // includes RcppArmadillo
#include "h-sampler.h"
#include "h-utils.h"
#include "mixture-state-sampler.h"
#include "simulation-smoother.h"

using namespace Rcpp;
 
NumericVector draw_h_auxiliary(const NumericVector y_star,
                               const NumericVector d,
                               const NumericVector s,
                               const double phi,
                               const double rho,
                               const double sigma2,
                               const double mu,
                               const CharacterVector centering,
                               const DataFrame mixing_constants) {
  //const double mu_mu = prior_mu(0);  // TODO
  //const double mu_sigma2 = pow(prior_mu(1), 2);
  const double mu_mu = -10;
  const double mu_sigma2 = 100;
  
  const NumericVector mixing_a = as<NumericVector>(mixing_constants["a"])[s];
  const NumericVector mixing_b = as<NumericVector>(mixing_constants["b"])[s];
  const NumericVector mixing_m = as<NumericVector>(mixing_constants["m"])[s];
  const NumericVector mixing_v = as<NumericVector>(mixing_constants["v"])[s];
  
  const List filter_result = aug_kalman_filter(phi, rho, sigma2, mixing_a, mixing_b, mixing_m, mixing_v, d, y_star, mu_mu, mu_sigma2, centering);
  
  const List smoothing_result = simulation_smoother(mu, filter_result, centering);
  const NumericVector eta = smoothing_result["eta"];
  const double eta0 = as<NumericVector>(smoothing_result["eta0"])(0);

  const int n = as<NumericVector>(filter_result["D"]).size();
  NumericVector h = rep(0.0, n);
  NumericVector dt;
  std::string scentering = as<std::string>(centering);
  if (scentering == "centered") {
    h(0) = mu + eta0;
    dt = mu*(1-phi) + rho*sqrt(sigma2)*d*mixing_a*exp(mixing_m/2);
  } else if (scentering == "non-centered") {
    h(0) = eta0;
    dt = rho*d*mixing_a*exp(mixing_m/2);
  }

  for (int i = 0; i < n-1; i++) {
    h(i+1) = dt(i) + phi*h(i) + eta(i);
  }

  return h;
}

NumericVector draw_latent_auxiliaryMH(const NumericVector y,
                                  const NumericVector y_star,
                                  const NumericVector d,
                                  const NumericVector h,
                                  const double phi,
                                  const double rho,
                                  const double sigma2,
                                  const double mu,
                                  //const CharacterVector centering,
                                  const DataFrame mixing_constants) {
  
  // Draw h from AUX
  const NumericVector s = draw_s_auxiliary(y_star, d, h, phi, rho, sigma2, mu, "centered", mixing_constants);
  const NumericVector proposed = draw_h_auxiliary(y_star, d, s, phi, rho, sigma2, mu, "centered", mixing_constants);

  // Calculate MH acceptance ratio
  const NumericVector mixing_a = mixing_constants["a"];
  const NumericVector mixing_b = mixing_constants["b"];
  const NumericVector mixing_m = mixing_constants["m"];
  const NumericVector mixing_p = mixing_constants["p"];
  const NumericVector mixing_v = mixing_constants["v"];
  const double log_acceptance = h_log_posterior(proposed, y, phi, rho, sigma2, mu) - h_log_posterior(h, y, phi, rho, sigma2, mu) -
    (h_aux_log_posterior(proposed, y_star, d, mixing_a, mixing_b, mixing_m, mixing_p, mixing_v, phi, rho, sigma2, mu) -
     h_aux_log_posterior(h, y_star, d, mixing_a, mixing_b, mixing_m, mixing_p, mixing_v, phi, rho, sigma2, mu));
  NumericVector result;
  if (log_acceptance > 0 || exp(log_acceptance) > R::runif(0, 1)) {
    result = proposed;
  } else {
    result = h;
  }

  return result;
}

