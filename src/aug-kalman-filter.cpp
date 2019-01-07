#include <Rcpp.h>
#include "aug-kalman-filter.h"
#include "parameterization.hpp"

using namespace Rcpp;

List aug_kalman_filter(const double phi, const double rho, const double sigma2,
                       const NumericVector a, const NumericVector b, const NumericVector m, const NumericVector v,
                       const NumericVector d, const NumericVector y_star,
                       const double mu_mu, const double sigma2_mu, const Parameterization centering) {
  List result;
  switch (centering) {
    case Parameterization::CENTERED:
    result = aug_kalman_filter_c(phi, rho, sigma2, a, b, m, v, d, y_star, mu_mu, sigma2_mu);
    break;
    case Parameterization::NONCENTERED:
    result = aug_kalman_filter_nc(phi, rho, sigma2, a, b, m, v, d, y_star, mu_mu, sigma2_mu);
    break;
  }
  return result;
}

List aug_kalman_filter_c(const double phi, const double rho, const double sigma2,
                         const NumericVector a, const NumericVector b, const NumericVector m, const NumericVector v,
                         const NumericVector d, const NumericVector y_star,
                         const double mu_mu, const double sigma2_mu) {
  const R_xlen_t n = y_star.size();
  
  const double sigma = sqrt(sigma2);
  const NumericVector gamma_ts = d * rho * sigma * exp(m/2);
  const NumericVector h_ts = b * v * gamma_ts;
  const double j_t22 = sigma2 * (1-rho*rho);
  const double h1_var = sigma2/(phi*phi > 1-1e-8 ? 1e-8 : (1-phi*phi));
  
  double a_star = 0;
  double A_star = -1;
  double P = h1_var;
  
  // Init returned values
  NumericVector D(n);
  NumericVector J(n);
  NumericVector L(n);
  NumericVector f(n);
  NumericVector F(n);
  double Q = 1/sigma2_mu;
  double q = mu_mu * Q;
  
  // Init for the loop
  double K;
  
  // Main loop
  for (int i = 0; i < int(n); i++) {
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
  
  return List::create(
    _["sigma2"] = sigma2,
    _["D"] = D,
    _["J1"] = J,
    _["L"] = L,
    _["f"] = f,
    _["F"] = F,
    _["hts"] = h_ts,
    _["v"] = v,
    _["Q"] = Q,
    _["q"] = q,
    _["jt22"] = j_t22,
    _["h1var"] = h1_var
  );
}

List aug_kalman_filter_nc(const double phi, const double rho, const double sigma2,
                          const NumericVector a, const NumericVector b, const NumericVector m, const NumericVector v,
                          const NumericVector d, const NumericVector y_star,
                          const double mu_mu, const double sigma2_mu) {
  const R_xlen_t n = y_star.size();
  
  const double sigma = sqrt(sigma2);
  const NumericVector gamma_ts = d * rho * exp(m/2);
  const NumericVector h_ts = b * v * gamma_ts;
  const double j_t22 = (1-rho*rho);
  const double h1_var = 1/(phi*phi > 1-1e-8 ? 1e-8 : (1-phi*phi));
  
  double a_star = 0;
  double A_star = 0;
  double P = h1_var;
  
  // Init returned values
  NumericVector D(n);
  NumericVector J(n);
  NumericVector L(n);
  NumericVector f(n);
  NumericVector F(n);
  double Q = 1/sigma2_mu;
  double q = mu_mu * Q;
  
  // Init for the loop
  double K;
  
  // Main loop
  for (int i = 0; i < int(n); i++) {
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
  
  return List::create(
    _["sigma2"] = sigma2,
    _["D"] = D,
    _["J1"] = J,
    _["L"] = L,
    _["f"] = f,
    _["F"] = F,
    _["hts"] = h_ts,
    _["v"] = v,
    _["Q"] = Q,
    _["q"] = q,
    _["jt22"] = j_t22,
    _["h1var"] = h1_var
  );
}
