// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <string>
#include <cmath>
#include "theta-sampler.h"
#include "theta-utils.h"
#include "aug-kalman-filter.h"
using namespace Rcpp;

NumericVector draw_theta_rwMH(const double phi, const double rho,
                              const double sigma2, const double mu,
                              const NumericVector y, const NumericVector h,
                              const NumericVector prior_phi,
                              const NumericVector prior_rho,
                              const NumericVector prior_sigma2,
                              const NumericVector prior_mu,
                              const CharacterVector centering,
                              const double stdev) {
  NumericVector result;
  std::string scentering = as<std::string>(centering);
  const NumericVector proposed = theta_propose(phi, rho, sigma2, mu, y, h, stdev);
  const double phi_prop = proposed(0), rho_prop = proposed(1), sigma2_prop = proposed(2),
    mu_prop = proposed(3), prop_old_logdens = proposed(4), prop_new_logdens = proposed(5);
  const double log_acceptance = (theta_log_prior(phi_prop, rho_prop, sigma2_prop, mu_prop, prior_phi, prior_rho, prior_sigma2, prior_mu) +
    theta_log_likelihood(phi_prop, rho_prop, sigma2_prop, mu_prop, y, h, centering)) -
    (theta_log_prior(phi, rho, sigma2, mu, prior_phi, prior_rho, prior_sigma2, prior_mu) +
    theta_log_likelihood(phi, rho, sigma2, mu, y, h, centering)) -
    (prop_new_logdens - prop_old_logdens);
  if (log_acceptance > 0 || exp(log_acceptance) > R::runif(0, 1)) {
    result = NumericVector::create(phi_prop, rho_prop, sigma2_prop, mu_prop);
  } else {
    result = NumericVector::create(phi, rho, sigma2, mu);
  }
  return result;
}

NumericVector draw_theta_auxiliary(const double phi, const double rho,
                                   const double sigma2, const double mu,
                                   const NumericVector y_star,
                                   const NumericVector d,
                                   const NumericVector s,
                                   const NumericVector prior_phi,
                                   const NumericVector prior_rho,
                                   const NumericVector prior_sigma2,
                                   const NumericVector prior_mu,
                                   const DataFrame mixing_constants) {
  //::Rf_error("nyi");
  
  Environment stats = Environment::namespace_env("stats");
  Function optim = stats["optim"];
  Environment numDeriv = Environment::namespace_env("numDeriv");
  Function genD = numDeriv["genD"];
  Environment Matrix = Environment::namespace_env("Matrix");
  Function nearPD = Matrix["nearPD"];
  const NumericVector mixing_a = as<NumericVector>(mixing_constants["a"])[s];
  const NumericVector mixing_b = as<NumericVector>(mixing_constants["b"])[s];
  const NumericVector mixing_m = as<NumericVector>(mixing_constants["m"])[s];
  const NumericVector mixing_v = as<NumericVector>(mixing_constants["v"])[s];
  NumericVector _auxtheta_old = NumericVector::create(phi, rho, sigma2);
  const arma::vec auxtheta_old(_auxtheta_old.begin(), 3);
  const double eps = 1.0e-5;
  const List optim_result = optim(_["par"] = _auxtheta_old,
                                       _["fn"] = InternalFunction(&auxtheta_log_posterior),
                                       _["gr"] = R_NilValue,
                                       _["a"] = mixing_a,
                                       _["b"] = mixing_b,
                                       _["mm"] = mixing_m,
                                       _["v"] = mixing_v,
                                       _["d"] = d,
                                       _["y_star"] = y_star,
                                       _["prior_phi"] = prior_phi,
                                       _["prior_rho"] = prior_rho,
                                       _["prior_sigma2"] = prior_sigma2,
                                       _["prior_mu"] = prior_mu,
                                       _["method"] = "L-BFGS-B",
                                       _["control"] = List::create(_["fnscale"] = -1,
                                                                   _["trace"] = 0),
                                       _["lower"] = NumericVector::create(-1+eps, -1+eps, 0+eps),
                                       _["upper"] = NumericVector::create(1-eps, 1-eps, R_PosInf));
  //print(optim_result);
  const double optim_value = optim_result["value"];
  NumericVector _auxtheta_hat = optim_result["par"];
  //print(_auxtheta_hat);
  const double phi_hat = _auxtheta_hat(0), rho_hat = _auxtheta_hat(1), sigma2_hat = _auxtheta_hat(2);
  NumericMatrix _deriv_result = as<List>(genD(_["func"] = InternalFunction(&auxtheta_log_posterior),
                                              _["x"] = _auxtheta_hat,
                                              _["method.args"] = List::create(_["d"] = eps/2, _["eps"] = eps/2),
                                              _["a"] = mixing_a,
                                              _["b"] = mixing_b,
                                              _["mm"] = mixing_m,
                                              _["v"] = mixing_v,
                                              _["d"] = d,
                                              _["y_star"] = y_star,
                                              _["prior_phi"] = prior_phi,
                                              _["prior_rho"] = prior_rho,
                                              _["prior_sigma2"] = prior_sigma2,
                                              _["prior_mu"] = prior_mu))["D"];
  const arma::vec deriv_result(_deriv_result.begin(), 9);
  const arma::vec auxtheta_hat(_auxtheta_hat.begin(), 3);
  arma::vec grad = deriv_result.subvec(0, 2);
  arma::mat hessian(3, 3);
  hessian(0, 0) = deriv_result(3 + 0);
  hessian(1, 1) = deriv_result(3 + 2);
  hessian(2, 2) = deriv_result(3 + 5);
  hessian(0, 1) = hessian(1, 0) = deriv_result(3 + 1);
  hessian(0, 2) = hessian(2, 0) = deriv_result(3 + 3);
  hessian(1, 2) = hessian(2, 1) = deriv_result(3 + 4);
  //print(wrap(grad));
  try {
    hessian = -arma::mat(NumericVector(as<S4>(as<List>(nearPD(_["x"] = wrap(-hessian), _["conv.tol"] = wrap(1e-8), _["eig.tol"] = wrap(1e-5), _["maxit"] = wrap(150)))["mat"]).slot("x")).begin(), 3, 3);
  } catch (...) {
    hessian = -arma::eye(3, 3)/100;
    grad = grad*0;
    //Rcout << "+";
  }
  arma::mat sigma_star = arma::inv_sympd(-hessian);
  arma::vec auxtheta_star = auxtheta_hat + sigma_star*grad;
  //print(wrap(sigma_star));
  //print(wrap(auxtheta_star));
  arma::mat sigma_star_low_chol = arma::chol(sigma_star, "lower");
  bool found = false;
  arma::vec proposed;
  for (int i = 0; i < 20000 && !found; i++) {
    proposed = sigma_star_low_chol * arma::vec(as<NumericVector>(rnorm(3)).begin(), 3) + auxtheta_star;
    if (-1+eps < proposed(0) && proposed(0) < 1-eps &&
        -1+eps < proposed(1) && proposed(1) < 1-eps &&
        0+eps < proposed(2)) {
      found = true;
    }
    
    if (i == 10000-1 && !found) {
      grad = grad*0;
      sigma_star_low_chol = arma::eye(3, 3) * 10;
      sigma_star = sigma_star_low_chol * sigma_star_low_chol.t();
      auxtheta_star = auxtheta_hat + sigma_star*grad;
      //Rcout << "&";
    }
  }
  if (!found) {
    Rcout << phi << " " << rho << " " << sigma2 << " " << mu << std::endl;
    print(optim_result);
    print(_deriv_result);
    print(wrap(hessian));
    ::Rf_error("could not propose new value");
  }
  const NumericVector _proposed(proposed.begin(), proposed.end());
  const double log_posterior_old = auxtheta_log_posterior(_auxtheta_old,
                                                          mixing_a, mixing_b, mixing_m, mixing_v,
                                                          d, y_star,
                                                          prior_phi, prior_rho, prior_sigma2, prior_mu);
  const double log_posterior_new = auxtheta_log_posterior(_proposed,
                                                          mixing_a, mixing_b, mixing_m, mixing_v,
                                                          d, y_star,
                                                          prior_phi, prior_rho, prior_sigma2, prior_mu);
  const double log_proposal_old = arma::dot(auxtheta_old-auxtheta_hat, grad) - 0.5*arma::as_scalar((auxtheta_old-auxtheta_hat).t()*arma::inv_sympd(sigma_star)*(auxtheta_old-auxtheta_hat));  // the hessian might be different already!
  const double log_proposal_new = arma::dot(proposed-auxtheta_hat, grad) - 0.5*arma::as_scalar((proposed-auxtheta_hat).t()*arma::inv_sympd(sigma_star)*(proposed-auxtheta_hat));
  
  const double log_acceptance_rate = (log_posterior_new - log_proposal_new) - (log_posterior_old - log_proposal_old);
  
  //Rcout << phi << " " << rho << " " << sigma2 << " " << mu << std::endl;
  //print(wrap(proposed));
  //Rcout << exp(log_acceptance_rate) << " " << log_posterior_new << " " << log_proposal_new-.5*arma::as_scalar(grad.t()*arma::inv_sympd(sigma_star)*grad) << " " << log_posterior_old << " " << log_proposal_old << std::endl;
  NumericVector auxtheta_new;
  if (log_acceptance_rate > 0 || exp(log_acceptance_rate) > R::runif(0, 1)) {
    auxtheta_new = _proposed;
  } else {
    auxtheta_new = _auxtheta_old;
  }
  
  // Draw mu
  const double phi_new = auxtheta_new(0), rho_new = auxtheta_new(1), sigma2_new = auxtheta_new(2);
  const List filter_result = aug_kalman_filter_c(phi_new, rho_new, sigma2_new,
                                                 mixing_a, mixing_b, mixing_m, mixing_v,
                                                 d, y_star,
                                                 prior_mu(0), pow(prior_mu(1), 2));
  const double mu_new = R::rnorm(as<double>(filter_result["q"]) / as<double>(filter_result["Q"]), 1/sqrt(as<double>(filter_result["Q"])));
  //Rcout << "===============" << std::endl;
  return NumericVector::create(phi_new, rho_new, sigma2_new, mu_new);
}

