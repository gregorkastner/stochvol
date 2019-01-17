#include <Rcpp.h>
#include <cmath>
#include "theta-utils.h"
#include "parameterization.hpp"

using namespace Rcpp;

double theta_log_likelihood(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const NumericVector& y,
    const NumericVector& h,
    const NumericVector& ht,
    const Parameterization centering) {
  double result;
  switch (centering) {
    case Parameterization::CENTERED:
    result = theta_log_likelihood_c(phi, rho, sigma2, mu, y, h);
    break;
    case Parameterization::NONCENTERED:
    result = theta_log_likelihood_nc(phi, rho, sigma2, mu, y, ht);
    break;
  }
  return result;
}

double theta_log_likelihood_c(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericVector& h) {
  const int n = y.size();
  const double sigma = sqrt(sigma2);
  double log_lik = 0;
  for (int i = 0; i < n; i++) {
    double h_mean, h_sd, y_mean, y_sd;
    if (i == 0) {
      h_mean = mu;
      h_sd = sigma/sqrt(1-phi*phi);
    } else {
      h_mean = mu+phi*(h[i-1]-mu) + rho*sigma*exp(-h[i-1]/2)*y[i-1];
      h_sd = sigma*sqrt(1-rho*rho);
    }
    y_mean = 0;
    y_sd = exp(h[i]/2);
    log_lik += R::dnorm(y[i], y_mean, y_sd, true) + R::dnorm(h[i], h_mean, h_sd, true);
  }
  return log_lik;
}

double theta_log_likelihood_nc(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericVector& ht) {
  const int n = y.size();
  const double sigma = sqrt(sigma2);
  double log_lik = 0;
  for (int i = 0; i < n; i++) {
    double h_mean, h_sd, y_mean, y_sd;
    if (i == 0) {
      h_mean = 0;
      h_sd = 1/sqrt(1-phi*phi);
    } else {
      h_mean = phi*ht[i-1];
      h_sd = 1;
    }
    if (i < n-1) {
      y_mean = exp((sigma*ht[i]+mu)/2)*rho*(ht[i+1]-phi*ht[i]);
      y_sd = exp((sigma*ht[i]+mu)/2)*sqrt(1-rho*rho);
    } else {
      y_mean = 0;
      y_sd = exp((sigma*ht[i] + mu)/2);
    }
    log_lik += R::dnorm(y[i], y_mean, y_sd, true) + R::dnorm(ht[i], h_mean, h_sd, true);
  }
  return log_lik;
}

double theta_log_prior(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const NumericVector& prior_phi,
    const NumericVector& prior_rho,
    const NumericVector& prior_sigma2,
    const NumericVector& prior_mu,
    const bool gammaprior) {
  // use variable names to clear the confusion about Gamma and InvGamma
  //const double gammarate = prior_sigma2(1);
  //const double gammascale = 1/gammarate;
  //const double invgammascale = gammarate;
  // Wikipedia notatio: prior_sigma2(1) is beta in both the case of Gamma and of InvGamma
  return (-log(2)) + R::dbeta((phi+1)/2, prior_phi[0], prior_phi[1], true) +
    (-log(2)) + R::dbeta((rho+1)/2, prior_rho[0], prior_rho[1], true) +
    R::dnorm(mu, prior_mu[0], prior_mu[1], true) +
    (gammaprior ?
      R::dgamma(sigma2, prior_sigma2[0], 1/prior_sigma2[1], true) :  // Gamma(shape, scale)
      1/R::dgamma(sigma2, prior_sigma2[0]+2, prior_sigma2[1]/(prior_sigma2[0]*(prior_sigma2[0]+1)), true));  // moment matched InvGamma
}

NumericVector theta_transform(
    const double f,
    const double r,
    const double s,
    const double m) {
  return {1-2/(exp(2*f)+1), 1-2/(exp(2*r)+1), exp(s), m};
}

NumericVector theta_transform_inv(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu) {
  return {0.5*log(2/(1-phi)-1), 0.5*log(2/(1-rho)-1), log(sigma2), mu};
}

double theta_transform_log_det_jac(
    const double f,
    const double r,
    const double s,
    const double m) {
  return 2*(log(4) + f+r+s/2 - log(exp(2*f)+1) - log(exp(2*r)+1));
}

double theta_transform_inv_log_det_jac(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu) {
  return -(log(1-phi*phi)+log(1-rho*rho)+log(sigma2));
}

NumericVector theta_proposal_stdev(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const NumericVector& y,
    const NumericVector& h,
    const NumericVector& ht,
    const double stdev) {
  return rep(stdev, 4);
}

NumericVector theta_propose(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const NumericVector& y,
    const NumericVector& h,
    const NumericVector& ht,
    const double stdev) {
  const NumericVector theta_old_t = theta_transform_inv(phi, rho, sigma2, mu);
  
  const NumericVector &proposal_mean_old = theta_old_t;
  const NumericVector proposal_sd_old = theta_proposal_stdev(phi, rho, sigma2, mu, y, h, ht, stdev);
  const NumericVector theta_new_t = rnorm(4)*proposal_sd_old + proposal_mean_old;
  const NumericVector theta_new = theta_transform(theta_new_t[0], theta_new_t[1], theta_new_t[2], theta_new_t[3]);
  const double phi_new = theta_new[0], rho_new = theta_new[1], sigma2_new = theta_new[2], mu_new = theta_new[3];
  double theta_density_new = theta_transform_inv_log_det_jac(phi_new, rho_new, sigma2_new, mu_new);
  for (int i = 0; i < 4; i++) {
    theta_density_new += R::dnorm(theta_new_t[i], proposal_mean_old[i], proposal_sd_old[i], true);
  }
  
  const NumericVector &proposal_mean_new = theta_new_t;
  const NumericVector proposal_sd_new = theta_proposal_stdev(phi_new, rho_new, sigma2_new, mu_new, y, h, ht, stdev);
  double theta_density_old = theta_transform_inv_log_det_jac(phi, rho, sigma2, mu);
  for (int i = 0; i < 4; i++) {
    theta_density_old += R::dnorm(theta_old_t[i], proposal_mean_new[i], proposal_sd_new[i], true);
  }
  
  return {phi_new, rho_new, sigma2_new, mu_new, theta_density_old, theta_density_new};
}

