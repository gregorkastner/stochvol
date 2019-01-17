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
#include "parameterization.hpp"

using namespace Rcpp;
 
void draw_h_auxiliary(
    NumericVector& h,
    List& cache,
    const NumericVector& y_star,
    const NumericVector& d,
    const NumericVector& s,
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const double priormu_mu,
    const double priormu_sigma,
    const Parameterization centering) {
  
  NumericVector mixing_a = cache["mixing_a"];
  NumericVector mixing_b = cache["mixing_b"];
  NumericVector mixing_m = cache["mixing_m"];
  NumericVector mixing_v = cache["mixing_v"];
  std::transform(s.cbegin(), s.cend(), mixing_a.begin(), [](const int selem) -> double {return mix_a[selem];});
  std::transform(s.cbegin(), s.cend(), mixing_b.begin(), [](const int selem) -> double {return mix_b[selem];});
  std::transform(s.cbegin(), s.cend(), mixing_m.begin(), [](const int selem) -> double {return mix_mean[selem];});
  std::transform(s.cbegin(), s.cend(), mixing_v.begin(), [](const int selem) -> double {return mix_var[selem];});
  
  aug_kalman_filter(cache, phi, rho, sigma2, mixing_a, mixing_b, mixing_m, mixing_v, d, y_star, priormu_mu, pow(priormu_sigma, 2), centering);
  
  double eta0;
  simulation_smoother(eta0, cache, mu, centering, sigma2);
  const NumericVector eta = cache["eta"];

  const int n = y_star.length();
  NumericVector dt;
  switch (centering) {
    case Parameterization::CENTERED:
    h[0] = mu + eta0;
    dt = mu*(1-phi) + rho*sqrt(sigma2)*d*mixing_a*exp(mixing_m/2);
    break;
    case Parameterization::NONCENTERED:
    h[0] = eta0;
    dt = rho*d*mixing_a*exp(mixing_m/2);
    break;
  }

  for (int i = 0; i < n-1; i++) {
    h[i+1] = dt[i] + phi*h[i] + eta[i];
  }

  return;
}

void draw_latent_auxiliaryMH(
    NumericVector& h,
    List& cache,
    const NumericVector& y,
    const NumericVector& y_star,
    const NumericVector& d,
    const NumericVector& ht,
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const double priormu_mu,
    const double priormu_sigma) {
    //const CharacterVector centering,
  
  // Draw h from AUX
  NumericVector s = cache["s"];
  NumericVector proposed = cache["proposed"];
  draw_s_auxiliary(s, cache, y_star, d, h, ht, phi, rho, sigma2, mu, Parameterization::CENTERED);
  draw_h_auxiliary(proposed, cache, y_star, d, s, phi, rho, sigma2, mu, priormu_mu, priormu_sigma, Parameterization::CENTERED);

  // Calculate MH acceptance ratio
  const double log_acceptance = h_log_posterior(proposed, y, phi, rho, sigma2, mu) - h_log_posterior(h, y, phi, rho, sigma2, mu) -
    (h_aux_log_posterior(proposed, y_star, d, phi, rho, sigma2, mu) -
     h_aux_log_posterior(h, y_star, d, phi, rho, sigma2, mu));
  if (log_acceptance > 0 || exp(log_acceptance) > R::runif(0, 1)) {
    h = proposed;
  }

  return;
}

