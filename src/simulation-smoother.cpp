#include <Rcpp.h>
#include "simulation-smoother.h"
#include "parameterization.hpp"

using namespace Rcpp;

void simulation_smoother(
    double& eta0,
    Rcpp::List& cache,
    const double mu,
    const Parameterization centering,
    const double sigma2) {

  Rcpp::NumericVector eta = cache["eta"];
  switch (centering) {
    case Parameterization::CENTERED:
    simulation_smoother_c(eta, eta0, mu, cache);
    break;
    case Parameterization::NONCENTERED:
    simulation_smoother_nc(eta, eta0, mu, cache, sigma2);
    break;
  }
  return;
}

void simulation_smoother_c(
    Rcpp::NumericVector& eta,
    double& eta0,
    const double mu,
    const List& filter_results) {
  const NumericVector D_ret = filter_results["D"];
  const NumericVector J = filter_results["J1"];
  const NumericVector L = filter_results["L"];
  const NumericVector H_ts = filter_results["hts"];
  const NumericVector v = filter_results["mixing_v"];
  const double j_t22 = filter_results["jt22"];
  const double H_1_var = filter_results["h1var"];
  const NumericVector F = filter_results["F"];
  const NumericVector f = filter_results["f"];
  const NumericVector E = f - mu * F;
  
  double r = 0;
  double U = 0;
  const NumericVector D_inv = 1/D_ret;
  double hjpj, C, Eps, V;
  for (int i = D_ret.size()-1; i >= 0; i--) {
    hjpj = H_ts[i]*J[i] + j_t22;
    C = (pow(H_ts[i], 2) + j_t22) - pow(H_ts[i]*v[i], 2)*D_inv[i] - U*pow(H_ts[i] * J[i] + j_t22, 2);
    Eps = R::rnorm(0, sqrt(C));
    V = H_ts[i]*v[i]*D_inv[i] + U*L[i]*hjpj;
    eta[i] = H_ts[i]*v[i]*E[i]*D_inv[i] + hjpj*r + Eps;
    
    r = E[i]*D_inv[i] + L[i]*r - V*Eps/C;
    U = D_inv[i] + U*pow(L[i], 2) + pow(V, 2)/C;
  }
  
  // Case i = 0
  C = H_1_var*(1 - U*H_1_var);
  Eps = R::rnorm(0, sqrt(C));
  eta0 = H_1_var*r + Eps;
}

void simulation_smoother_nc(
    Rcpp::NumericVector& eta,
    double& eta0,
    const double mu,
    const List& filter_results,
    const double sigma2) {
  const NumericVector D_ret = filter_results["D"];
  const NumericVector J = filter_results["J1"];
  const NumericVector L = filter_results["L"];
  const NumericVector H_ts = filter_results["hts"];
  const NumericVector v = filter_results["mixing_v"];
  const double j_t22 = filter_results["jt22"];
  const double H_1_var = filter_results["h1var"];
  const NumericVector F = filter_results["F"];
  const NumericVector f = filter_results["f"];
  const NumericVector E = f - mu * F;
  
  const double sigma = sqrt(sigma2);
  
  double r = 0;
  double U = 0;
  const NumericVector D_inv = 1/D_ret;
  double hjpj, C, Eps, V;
  for (int i = D_ret.size()-1; i >= 0; i--) {
    hjpj = H_ts[i]*J[i] + j_t22;
    C = (pow(H_ts[i], 2) + j_t22) - pow(H_ts[i]*v[i], 2)*D_inv[i] - U*pow(H_ts[i] * J[i] + j_t22, 2);
    Eps = R::rnorm(0, sqrt(C));
    V = sigma*H_ts[i]*v[i]*D_inv[i] + U*L[i]*hjpj;
    eta[i] = H_ts[i]*v[i]*E[i]*D_inv[i] + hjpj*r + Eps;
    
    r = sigma*E[i]*D_inv[i] + L[i]*r - V*Eps/C;
    U = sigma2*D_inv[i] + U*pow(L[i], 2) + pow(V, 2)/C;
  }
  
  // Case i = 0
  C = H_1_var*(1 - U*H_1_var);
  Eps = R::rnorm(0, sqrt(C));
  eta0 = H_1_var*r + Eps;
  
  return;
}
