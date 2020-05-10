#include <RcppArmadillo.h>
#include <cmath>
#include "theta-utils.h"
#include "parameterization.h"

using namespace Rcpp;

// Helper declarations

double theta_log_likelihood_c(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const arma::vec& y,
    const arma::vec& h,
    const arma::vec& exp_h_half);

double theta_log_likelihood_nc(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const arma::vec& y,
    const arma::vec& ht,
    const arma::vec& exp_h_half_tilde);

arma::vec4 grad_theta_log_posterior(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const arma::vec& y,
    const arma::vec& h,
    const arma::vec2& prior_phi,
    const arma::vec2& prior_rho,
    const arma::vec2& prior_sigma2,
    const arma::vec2& prior_mu);

// Definitions

double theta_log_likelihood(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const arma::vec& y,
    const arma::vec& h,
    const arma::vec& ht,
    const arma::vec& exp_h_half,
    const Parameterization centering) {
  double result = 0;
  switch (centering) {
    case Parameterization::CENTERED:
    result = theta_log_likelihood_c(phi, rho, sigma2, mu, y, h, exp_h_half);
    break;
    case Parameterization::NONCENTERED:
    result = theta_log_likelihood_nc(phi, rho, sigma2, mu, y, ht, exp_h_half);
    break;
  }
  return result;
}

double dnorm_log(
    const double x,
    const double mu,
    const double sd,
    const double log_sd) {
  const double z = (x - mu) / sd;
  return -.5 * z * z - log_sd;
}

double theta_log_likelihood_c(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const arma::vec& y,
    const arma::vec& h,
    const arma::vec& exp_h_half) {
  const int n = y.size();
  const double sigma = std::sqrt(sigma2);
  const double h_sd_t = sigma * std::sqrt(1 - rho * rho);
  const double log_h_sd_t = std::log(h_sd_t);
  double log_lik = 0;
  for (int i = 0; i < n; i++) {
    double h_mean, h_sd, y_mean, y_sd, log_h_sd, log_y_sd;
    if (i == 0) {
      h_mean = mu;
      h_sd = sigma / std::sqrt(1 - phi * phi);
      log_h_sd = std::log(h_sd);
    } else {
      h_mean = mu + phi * (h[i-1] - mu) + rho * sigma / exp_h_half[i-1] * y[i-1];
      h_sd = h_sd_t;
      log_h_sd = log_h_sd_t;
    }
    y_mean = 0;
    y_sd = exp_h_half[i];
    log_y_sd = .5 * h[i];
    log_lik += dnorm_log(y[i], y_mean, y_sd, log_y_sd) + dnorm_log(h[i], h_mean, h_sd, log_h_sd);
  }

  return log_lik;
}

double theta_log_likelihood_nc(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const arma::vec& y,
    const arma::vec& ht,
    const arma::vec& exp_h_half) {
  const int n = y.size();
  const double sigma = std::sqrt(sigma2);
  const double rho_const = std::sqrt(1 - rho * rho);
  const double log_rho_const = .5 * std::log(1 - rho * rho);
  double log_lik = 0;
  for (int i = 0; i < n; i++) {
    double h_mean, h_sd, y_mean, y_sd, log_h_sd, log_y_sd;
    if (i == 0) {
      h_mean = 0;
      h_sd = 1 / std::sqrt(1 - phi * phi);
      log_h_sd = -.5 * std::log(1 - phi * phi);
    } else {
      h_mean = phi * ht[i-1];
      h_sd = 1;
      log_h_sd = 0;
    }
    if (i < n-1) {
      y_mean = exp_h_half[i] * rho * (ht[i+1] - phi * ht[i]);
      y_sd = exp_h_half[i] * rho_const;
      log_y_sd = .5 * (sigma * ht[i] + mu) + log_rho_const;
    } else {
      y_mean = 0;
      y_sd = exp_h_half[i];
      log_y_sd = .5 * (sigma * ht[i] + mu);
    }
    log_lik += dnorm_log(y[i], y_mean, y_sd, log_y_sd) + dnorm_log(ht[i], h_mean, h_sd, log_h_sd);
  }

  return log_lik;
}

double dlogpriornorm(const double x, const double mean, const double stdev) {  // only necessary terms are added
  return -.5 * std::pow((x - mean) / stdev, 2);
}

double dlogpriorbeta(const double x, const double alpha, const double beta) {
  return (alpha - 1.) * std::log(x) + (beta - 1.) * std::log(1. - x);
}

double dlogpriorgamma(const double x, const double alpha, const double beta) {
  return (alpha - 1.) * std::log(x) - beta * x;
}

double dlogpriorinvgamma(const double x, const double alpha, const double beta) {
  return -(alpha - 1.) * std::log(x) - beta / x;
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
  return dlogpriornorm(mu, prior_mu[0], prior_mu[1]) +
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
  return dlogpriorbeta(.5 * (phi + 1.), prior_phi[0], prior_phi[1]) +
    dlogpriorbeta(.5 * (rho + 1.), prior_rho[0], prior_rho[1]) +
    (gammaprior ?
      dlogpriorgamma(sigma2, prior_sigma2[0], prior_sigma2[1]) :
      // moment matched InvGamma
      (-2. * std::log(sigma2) + dlogpriorinvgamma(1 / sigma2, prior_sigma2[0] + 2, prior_sigma2[1] / (prior_sigma2[0] * (prior_sigma2[0] + 1.)))));
}

arma::vec4 theta_transform(
    const double f,
    const double r,
    const double s,
    const double m) {
  return {1 - 2/(std::exp(2 * f) + 1), 1 - 2/(std::exp(2 * r) + 1), std::exp(s), m};
}

arma::vec4 theta_transform_inv(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu) {
  return {0.5 * std::log(2. / (1 - phi) - 1), 0.5 * std::log(2. / (1 - rho) - 1), std::log(sigma2), mu};
}

double theta_transform_log_det_jac(
    const double f,
    const double r,
    const double s,
    const double m) {
  static const double log4 = std::log(4.);
  return 2 * (log4 + f + r + s / 2. - std::log(std::exp(2. * f) + 1.) - std::log(std::exp(2. * r) + 1.));
}

double theta_transform_inv_log_det_jac(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu) {
  return -(std::log((1. - phi * phi) * (1. - rho * rho) * sigma2));
}

arma::vec rnorm4_arma() {
  return {R::rnorm(0, 1), R::rnorm(0, 1), R::rnorm(0, 1), R::rnorm(0, 1)};
}

arma::vec6 theta_propose_rwmh(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const arma::vec& y,
    const arma::vec& h,
    const arma::vec& ht,
    const stochvol::Adaptation<4>::Result& adaptation_proposal) {
  const arma::vec4 theta_old_t = theta_transform_inv(phi, rho, sigma2, mu);
  
  const arma::vec4 &proposal_mean_old = theta_old_t;
  const arma::vec4 theta_new_t_standardized = rnorm4_arma();
  const arma::vec4 theta_new_t =
    std::sqrt(adaptation_proposal.scale) * adaptation_proposal.covariance_chol * theta_new_t_standardized +
    proposal_mean_old;
  const arma::vec4 theta_new = theta_transform(theta_new_t[0], theta_new_t[1], theta_new_t[2], theta_new_t[3]);
  const double phi_new = theta_new[0], rho_new = theta_new[1], sigma2_new = theta_new[2], mu_new = theta_new[3];
  double theta_density_new = theta_transform_inv_log_det_jac(phi_new, rho_new, sigma2_new, mu_new);
  for (int i = 0; i < 4; i++) {
    theta_density_new += -.5 * std::pow(theta_new_t_standardized[i], 2);  // dnorm(x, 0, 1, true)
  }
  
  const arma::vec4 &proposal_mean_new = theta_new_t;
  const arma::vec4 theta_old_t_standardized =
    1 / std::sqrt(adaptation_proposal.scale) * adaptation_proposal.covariance_chol_inv *
    (theta_old_t - proposal_mean_new);
  double theta_density_old = theta_transform_inv_log_det_jac(phi, rho, sigma2, mu);
  for (int i = 0; i < 4; i++) {
    theta_density_old += -.5 * std::pow(theta_old_t_standardized[i], 2);  // dnorm(x, 0, 1, true)
  }
  
  return {phi_new, rho_new, sigma2_new, mu_new, theta_density_old, theta_density_new};
}

double dmvnorm_mala(
    const arma::vec4& x_prime,
    const arma::vec4& x,
    const arma::vec4& grad_log_posterior,
    const arma::mat44& A,  // preconditioning matrix
    const arma::mat44& A_inv,
    const double tau,
    const arma::vec y,
    const arma::vec h,
    const arma::vec2& prior_phi,
    const arma::vec2& prior_rho,
    const arma::vec2& prior_sigma2,
    const arma::vec2& prior_mu,
    const bool log = false) {
  const auto mean = x + tau * A * grad_log_posterior;
  const arma::vec4 demeaned = x_prime - mean;
  const double log_density = -0.25 * arma::as_scalar(demeaned.t() * A_inv * demeaned) / tau;
  return log ? log_density : std::exp(log_density);
}

arma::mat44 crossprod(const arma::mat44 m) { return m.t() * m; }
arma::mat44 tcrossprod(const arma::mat44 m) { return m * m.t(); }

arma::vec6 theta_propose_mala(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const arma::vec& y,
    const arma::vec& h,
    const arma::vec2& prior_phi,
    const arma::vec2& prior_rho,
    const arma::vec2& prior_sigma2,
    const arma::vec2& prior_mu,
    const stochvol::Adaptation<4>::Result& adaptation_proposal) {
  const arma::vec4 theta_old_t = theta_transform_inv(phi, rho, sigma2, mu);
  
  const arma::vec4 grad_old = grad_theta_log_posterior(phi, rho, sigma2, mu, y, h, prior_phi, prior_rho, prior_sigma2, prior_mu) %
    arma::vec4({1 - std::pow(phi, 2), 1 - std::pow(rho, 2), sigma2, 1});
  const arma::vec4 proposal_mean_old =
    theta_old_t +
    adaptation_proposal.scale * adaptation_proposal.covariance * grad_old;
  const arma::vec4 theta_new_t =
    std::sqrt(2.) * std::sqrt(adaptation_proposal.scale) * adaptation_proposal.covariance_chol * rnorm4_arma() +
    proposal_mean_old;
  const arma::vec4 theta_new = theta_transform(theta_new_t(0), theta_new_t(1), theta_new_t(2), theta_new_t(3));
  const double phi_new = theta_new(0), rho_new = theta_new(1), sigma2_new = theta_new(2), mu_new = theta_new(3);
  const arma::vec4 grad_new = grad_theta_log_posterior(phi_new, rho_new, sigma2_new, mu_new, y, h, prior_phi, prior_rho, prior_sigma2, prior_mu) %
    arma::vec4({1 - std::pow(phi_new, 2), 1 - std::pow(rho_new, 2), sigma2_new, 1});

  double theta_density_new = theta_transform_inv_log_det_jac(phi_new, rho_new, sigma2_new, mu_new) +
    dmvnorm_mala(theta_new_t, theta_old_t, grad_old,
        adaptation_proposal.covariance, adaptation_proposal.precision, adaptation_proposal.scale,
        y, h, prior_phi, prior_rho, prior_sigma2, prior_mu, true);
  double theta_density_old = theta_transform_inv_log_det_jac(phi, rho, sigma2, mu) +
    dmvnorm_mala(theta_old_t, theta_new_t, grad_new,
        adaptation_proposal.covariance, adaptation_proposal.precision, adaptation_proposal.scale,
        y, h, prior_phi, prior_rho, prior_sigma2, prior_mu, true);

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
  const arma::vec theta_old_t = theta_transform_inv(phi, rho, sigma2, mu).head(3);
  
  const arma::vec &proposal_mean_old = theta_old_t;
  const arma::vec theta_new_t_standardized = rnorm3_arma();
  const arma::vec theta_new_t = proposal_chol*theta_new_t_standardized + proposal_mean_old;
  const arma::vec theta_new = theta_transform(theta_new_t[0], theta_new_t[1], theta_new_t[2], mu).head(3);
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

arma::vec4 grad_theta_log_posterior(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const arma::vec& y,
    const arma::vec& h,
    const arma::vec2& prior_phi,
    const arma::vec2& prior_rho,
    const arma::vec2& prior_sigma2,
    const arma::vec2& prior_mu) {
  const int n = y.n_elem;
  const double sigma = sqrt(sigma2);
  double h_tilde = (h(0) - mu) / sigma;
  // Grad of log likelihood
  double d_phi = phi * (std::pow(h_tilde, 2) - 1 / (1 - std::pow(phi, 2)));
  double d_rho = 0;
  double d_sigma2 = .5 / sigma2 * ((1 - std::pow(phi, 2)) * std::pow(h_tilde, 2) - 1);
  double d_mu = (1 - std::pow(phi, 2)) * h_tilde / sigma;
  const double rho_const = 1 / (1 - std::pow(rho, 2));
  for (int t = 0; t < n-1; t++) {
    const double y_h_const = y[t] * std::exp(-h[t] / 2);
    h_tilde = (h[t] - mu) / sigma;
    const double Delta = (h[t+1] - mu) / sigma - phi * h_tilde;
    const double dryh_const = Delta - rho * y_h_const;
    d_phi += rho_const * h_tilde * dryh_const;
    d_rho += rho_const * (rho - rho_const * rho * (std::pow(y_h_const, 2) - 2 * rho * y_h_const * Delta + Delta * Delta) + y_h_const * Delta);
    d_sigma2 += .5 / sigma2 * (rho_const * Delta * dryh_const - 1);
    d_mu += rho_const / sigma * (1 - phi) * dryh_const;
  }
  // Grad of log prior
  const double phi_beta = .5 * (phi + 1), rho_beta = .5 * (rho + 1);
  d_phi += .5 * ((prior_phi[0] - 1) / phi_beta - (prior_phi[1] - 1) / (1 - phi_beta));
  d_rho += .5 * ((prior_rho[0] - 1) / rho_beta - (prior_rho[1] - 1) / (1 - rho_beta));
  d_sigma2 += (prior_sigma2[0] - 1) / sigma2 - prior_sigma2[1];
  d_mu += -(mu - prior_mu[0]) / std::pow(prior_mu[1], 2);
  return {d_phi, d_rho, d_sigma2, d_mu};
}

