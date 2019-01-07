#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <string>
#include "aug-kalman-filter.h"
#include "h-sampler.h"
#include "h-utils.h"
#include "mixture-state-sampler.h"
#include "simulation-smoother.h"
#include "auxmix.h"

using namespace Rcpp;
 
NumericVector draw_h_auxiliary(const NumericVector y_star,
                               const NumericVector d,
                               const NumericVector s,
                               const double phi,
                               const double rho,
                               const double sigma2,
                               const double mu,
                               const double priormu_mu,
                               const double priormu_sigma,
                               const CharacterVector centering) {
  NumericVector mixing_a(s.length()); std::transform(s.cbegin(), s.cend(), mixing_a.begin(), [](const int selem) -> double {return mix_a[selem];});
  NumericVector mixing_b(s.length()); std::transform(s.cbegin(), s.cend(), mixing_b.begin(), [](const int selem) -> double {return mix_b[selem];});
  NumericVector mixing_m(s.length()); std::transform(s.cbegin(), s.cend(), mixing_m.begin(), [](const int selem) -> double {return mix_mean[selem];});
  NumericVector mixing_v(s.length()); std::transform(s.cbegin(), s.cend(), mixing_v.begin(), [](const int selem) -> double {return mix_var[selem];});
  
  const List filter_result = aug_kalman_filter(phi, rho, sigma2, mixing_a, mixing_b, mixing_m, mixing_v, d, y_star, priormu_mu, pow(priormu_sigma, 2), centering);
  
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
                                  const double priormu_mu,
                                  const double priormu_sigma) {
                                  //const CharacterVector centering,
  
  // Draw h from AUX
  const NumericVector s = draw_s_auxiliary(y_star, d, h, phi, rho, sigma2, mu, "centered");
  const NumericVector proposed = draw_h_auxiliary(y_star, d, s, phi, rho, sigma2, mu, priormu_mu, priormu_sigma, "centered");

  // Calculate MH acceptance ratio
  const double log_acceptance = h_log_posterior(proposed, y, phi, rho, sigma2, mu) - h_log_posterior(h, y, phi, rho, sigma2, mu) -
    (h_aux_log_posterior(proposed, y_star, d, phi, rho, sigma2, mu) -
     h_aux_log_posterior(h, y_star, d, phi, rho, sigma2, mu));
  NumericVector result;
  if (log_acceptance > 0 || exp(log_acceptance) > R::runif(0, 1)) {
    result = proposed;
  } else {
    result = h;
  }

  return result;
}

