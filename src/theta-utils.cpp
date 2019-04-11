#include <RcppArmadillo.h>
#include <cmath>
#include "theta-utils.h"
#include "parameterization.h"

using namespace Rcpp;

double theta_log_likelihood(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const arma::vec& y,
    const arma::vec& h,
    const arma::vec& ht,
    const Parameterization centering) {
  double result = 0;
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
    const arma::vec& y,
    const arma::vec& h) {
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
    const arma::vec& y,
    const arma::vec& ht) {
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
    const arma::vec& prior_phi,
    const arma::vec& prior_rho,
    const arma::vec& prior_sigma2,
    const arma::vec& prior_mu,
    const bool gammaprior) {
  return R::dnorm(mu, prior_mu[0], prior_mu[1], true) +
    thetamu_log_prior(phi, rho, sigma2, prior_phi, prior_rho, prior_sigma2, gammaprior);
}

double thetamu_log_prior(
    const double phi,
    const double rho,
    const double sigma2,
    const arma::vec& prior_phi,
    const arma::vec& prior_rho,
    const arma::vec& prior_sigma2,
    const bool gammaprior) {
  // use variable names to clear the confusion about Gamma and InvGamma
  //const double gammarate = prior_sigma2(1);
  //const double gammascale = 1/gammarate;
  //const double invgammascale = gammarate;
  // Wikipedia notation: prior_sigma2(1) is beta in both the case of Gamma and of InvGamma
  return (-log(2.)) + R::dbeta((phi+1)/2, prior_phi[0], prior_phi[1], true) +
    (-log(2.)) + R::dbeta((rho+1)/2, prior_rho[0], prior_rho[1], true) +
    (gammaprior ?
      R::dgamma(sigma2, prior_sigma2[0], 1/prior_sigma2[1], true) :  // Gamma(shape, scale)
      1/R::dgamma(sigma2, prior_sigma2[0]+2, prior_sigma2[1]/(prior_sigma2[0]*(prior_sigma2[0]+1)), true));  // moment matched InvGamma
}

arma::vec theta_transform(
    const double f,
    const double r,
    const double s,
    const double m) {
  return {1-2/(exp(2*f)+1), 1-2/(exp(2*r)+1), exp(s), m};
}

arma::vec theta_transform_inv(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu) {
  return {0.5*log(2./(1-phi)-1), 0.5*log(2./(1-rho)-1), log(sigma2), mu};
}

double theta_transform_log_det_jac(
    const double f,
    const double r,
    const double s,
    const double m) {
  return 2*(log(4.) + f+r+s/2. - log(exp(2.*f)+1.) - log(exp(2.*r)+1.));
}

double theta_transform_inv_log_det_jac(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu) {
  return -(log(1.-phi*phi)+log(1.-rho*rho)+log(sigma2));
}

arma::vec rnorm4_arma() {
  return {R::rnorm(0, 1), R::rnorm(0, 1), R::rnorm(0, 1), R::rnorm(0, 1)};
}

arma::vec theta_propose(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const arma::vec& y,
    const arma::vec& h,
    const arma::vec& ht,
    const arma::mat& proposal_chol,
    const arma::mat& proposal_chol_inv) {
  const arma::vec theta_old_t = theta_transform_inv(phi, rho, sigma2, mu);
  
  const arma::vec &proposal_mean_old = theta_old_t;
  const arma::vec theta_new_t_standardized = rnorm4_arma();
  const arma::vec theta_new_t = proposal_chol*theta_new_t_standardized + proposal_mean_old;
  const arma::vec theta_new = theta_transform(theta_new_t[0], theta_new_t[1], theta_new_t[2], theta_new_t[3]);
  const double phi_new = theta_new[0], rho_new = theta_new[1], sigma2_new = theta_new[2], mu_new = theta_new[3];
  double theta_density_new = theta_transform_inv_log_det_jac(phi_new, rho_new, sigma2_new, mu_new);
  for (int i = 0; i < 4; i++) {
    theta_density_new += R::dnorm(theta_new_t_standardized[i], 0., 1., true);
  }
  
  const arma::vec &proposal_mean_new = theta_new_t;
  const arma::vec theta_old_t_standardized = proposal_chol_inv*(theta_old_t-proposal_mean_new);
  double theta_density_old = theta_transform_inv_log_det_jac(phi, rho, sigma2, mu);
  for (int i = 0; i < 4; i++) {
    theta_density_old += R::dnorm(theta_old_t_standardized[i], 0., 1., true);
  }
  
  return {phi_new, rho_new, sigma2_new, mu_new, theta_density_old, theta_density_new};
}

arma::vec rnorm3_arma() {
  return {R::rnorm(0, 1), R::rnorm(0, 1), R::rnorm(0, 1)};
}

arma::vec thetamu_propose(
    const double phi,
    const double rho,
    const double sigma2,
    const arma::vec& y,
    const arma::vec& h,
    const arma::vec& ht,
    const arma::mat& proposal_chol,
    const arma::mat& proposal_chol_inv) {
  const double mu = 0;
  const arma::vec theta_old_t = theta_transform_inv(phi, rho, sigma2, mu);
  
  const arma::vec &proposal_mean_old = theta_old_t.head(3);
  const arma::vec theta_new_t_standardized = rnorm3_arma();
  const arma::vec theta_new_t = proposal_chol*theta_new_t_standardized + proposal_mean_old;
  const arma::vec theta_new = theta_transform(theta_new_t[0], theta_new_t[1], theta_new_t[2], mu);
  const double phi_new = theta_new[0], rho_new = theta_new[1], sigma2_new = theta_new[2];
  double theta_density_new = theta_transform_inv_log_det_jac(phi_new, rho_new, sigma2_new, mu);
  for (int i = 0; i < 3; i++) {
    theta_density_new += R::dnorm(theta_new_t_standardized[i], 0., 1., true);
  }
  
  const arma::vec &proposal_mean_new = theta_new_t;
  const arma::vec theta_old_t_standardized = proposal_chol_inv*(theta_old_t-proposal_mean_new);
  double theta_density_old = theta_transform_inv_log_det_jac(phi, rho, sigma2, mu);
  for (int i = 0; i < 3; i++) {
    theta_density_old += R::dnorm(theta_old_t_standardized[i], 0., 1., true);
  }
  
  return {phi_new, rho_new, sigma2_new, theta_density_old, theta_density_new};
}

