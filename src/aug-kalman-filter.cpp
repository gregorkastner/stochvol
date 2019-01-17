#include <Rcpp.h>
#include "aug-kalman-filter.h"
#include "parameterization.hpp"
#include "auxmix.h"

using namespace Rcpp;

void aug_kalman_filter(
    List& cache,
    const double phi,
    const double rho,
    const double sigma2,
    const NumericVector& mixing_a,
    const NumericVector& mixing_b,
    const NumericVector& mixing_m,
    const NumericVector& mixing_v,
    const NumericVector& d,  // TODO integervector
    const NumericVector& y_star,
    const double mu_mu,
    const double sigma2_mu,
    const Parameterization centering) {

  switch (centering) {
    case Parameterization::CENTERED:
    aug_kalman_filter_c(cache, phi, rho, sigma2, mixing_a, mixing_b, mixing_m, mixing_v, d, y_star, mu_mu, sigma2_mu);
    break;
    case Parameterization::NONCENTERED:
    aug_kalman_filter_nc(cache, phi, rho, sigma2, mixing_a, mixing_b, mixing_m, mixing_v, d, y_star, mu_mu, sigma2_mu);
    break;
  }
  return;
}

void aug_kalman_filter_c(
    List& cache,
    const double phi,
    const double rho,
    const double sigma2,
    const NumericVector& a,
    const NumericVector& b,
    const NumericVector& m,
    const NumericVector& v,
    const NumericVector& d,
    const NumericVector& y_star,
    const double mu_mu,
    const double sigma2_mu) {
  const int n = y_star.size();
  
  // Init returned values
  NumericVector D = cache["D"];
  NumericVector J = cache["J1"];
  NumericVector L = cache["L"];
  NumericVector f = cache["f"];
  NumericVector F = cache["F"];
  NumericVector h_ts = cache["hts"];
  double Q = 1/sigma2_mu;
  double q = mu_mu * Q;

  const double sigma = sqrt(sigma2);
  const NumericVector gamma_ts = d * rho * sigma * exp(m/2);
  h_ts = b * v * gamma_ts;
  const double j_t22 = sigma2 * (1-rho*rho);
  const double h1_var = sigma2/(phi*phi > 1-1e-8 ? 1e-8 : (1-phi*phi));
  
  double a_star = 0;
  double A_star = -1;
  double P = h1_var;
  
  
  // Init for the loop
  double K;
  
  // Main loop
  for (int i = 0; i < n; i++) {
    D[i] = P + v[i]*v[i];
    K = (phi*P + h_ts[i]*v[i])/D[i];
    L[i] = phi - K;
    J[i] = h_ts[i] - K*v[i];
    P = phi*P*L[i] + h_ts[i]*J[i] + j_t22;
    f[i] = y_star[i] - m[i] - a_star;
    F[i] = -A_star;
    a_star = a[i]*gamma_ts[i] + phi*a_star + K*f[i];
    A_star = phi*(A_star+1) - 1 + K*F[i];
    
    q += F[i]*f[i]/D[i];
    Q += F[i]*F[i]/D[i];
  }

  cache["Q"] = Q;
  cache["q"] = q;
  cache["jt22"] = j_t22;
  cache["h1var"] = h1_var;
  
  return;
}

void aug_kalman_filter_nc(
    List& cache,
    const double phi,
    const double rho,
    const double sigma2,
    const NumericVector& a,
    const NumericVector& b,
    const NumericVector& m,
    const NumericVector& v,
    const NumericVector& d,
    const NumericVector& y_star,
    const double mu_mu,
    const double sigma2_mu) {
  
  // Init returned values
  NumericVector D = cache["D"];
  NumericVector J = cache["J1"];
  NumericVector L = cache["L"];
  NumericVector f = cache["f"];
  NumericVector F = cache["F"];
  NumericVector h_ts = cache["hts"];
  double Q = 1/sigma2_mu;
  double q = mu_mu * Q;

  const int n = y_star.size();
  
  const double sigma = sqrt(sigma2);
  const NumericVector gamma_ts = d * rho * exp(m/2);
  h_ts = b * v * gamma_ts;
  const double j_t22 = (1-rho*rho);
  const double h1_var = 1/(phi*phi > 1-1e-8 ? 1e-8 : (1-phi*phi));
  
  double a_star = 0;
  double A_star = 0;
  double P = h1_var;
  
  // Init for the loop
  double K;
  
  // Main loop
  for (int i = 0; i < n; i++) {
    D[i] = sigma2*P + v[i]*v[i];
    K = (sigma*phi*P + h_ts[i]*v[i])/D[i];
    L[i] = phi - sigma*K;
    J[i] = h_ts[i] - K*v[i];
    P = phi*P*L[i] + h_ts[i]*J[i] + j_t22;
    f[i] = y_star[i] - m[i] - sigma*a_star;
    F[i] = 1.0 - sigma*A_star;
    a_star = a[i]*gamma_ts[i] + phi*a_star + K*f[i];
    A_star = phi*A_star + K*F[i];
    
    q += F[i]*f[i]/D[i];
    Q += F[i]*F[i]/D[i];
  }

  cache["Q"] = Q;
  cache["q"] = q;
  cache["jt22"] = j_t22;
  cache["h1var"] = h1_var;
  
  return;
}
