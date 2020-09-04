#include <RcppArmadillo.h>
#include <cmath>
#include "utils_parameters.h"
#include "type_definitions.h"
#include "densities.h"

using namespace Rcpp;

namespace stochvol {

double newton_raphson(
    const double startval,
    const double sumtau,
    const int n,
    const double lambda,
    const double tol,
    const int maxiter) {
  double x = startval;
  double error = R_PosInf;
  double xnew;
  bool converged = false;
  for (int i = 0; !converged && i < maxiter; i++) {
    xnew = x - dlogdnu(x, sumtau, lambda, n) / ddlogdnu(x, n);
    error = std::abs(xnew - x);
    x = xnew;
    converged = error < tol;
  }
  if (!converged) x = NA_REAL;
  return x;
}

double theta_log_likelihood_c(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const arma::vec& y,
    const double h0,
    const arma::vec& h,
    const arma::vec& exp_h_half);

double theta_log_likelihood_nc(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const arma::vec& y,
    const double h0,
    const arma::vec& ht,
    const arma::vec& exp_h_half_tilde);

arma::vec4 grad_theta_log_posterior(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const arma::vec& y,
    const double h0,
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
    const double h0,
    const arma::vec& h,
    const arma::vec& ht,
    const arma::vec& exp_h_half,
    const Parameterization centering) {
  double result;
  switch (centering) {
    case Parameterization::CENTERED:
      result = theta_log_likelihood_c(phi, rho, sigma2, mu, y, h0, h, exp_h_half);
      break;
    case Parameterization::NONCENTERED:
      result = theta_log_likelihood_nc(phi, rho, sigma2, mu, y, h0, ht, exp_h_half);  // TODO non-center h0
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
    const double h0,  // TODO test
    const arma::vec& h,
    const arma::vec& exp_h_half) {
  const int n = y.size();
  const double sigma = std::sqrt(sigma2),
               h_sd_t = sigma * std::sqrt(1 - rho * rho),
               log_h_sd_t = std::log(h_sd_t),
               B0 = std::sqrt(sigma2 / (1 - std::pow(phi, 2)));  // TODO stationary
  double log_lik = logdnorm(h0, mu, B0);
  for (int i = 0; i < n; i++) {
    double h_mean, h_sd, y_mean, y_sd, log_h_sd, log_y_sd;
    if (i == 0) {
      h_mean = mu + phi * (h0 - mu);
      h_sd = sigma;
      log_h_sd = std::log(h_sd);
    } else {
      h_mean = mu + phi * (h[i-1] - mu) + rho * sigma / exp_h_half[i-1] * y[i-1];
      h_sd = h_sd_t;
      log_h_sd = log_h_sd_t;
    }
    y_mean = 0;
    y_sd = exp_h_half[i];
    log_y_sd = .5 * h[i];
    log_lik += logdnorm2(y[i], y_mean, y_sd, log_y_sd) + logdnorm2(h[i], h_mean, h_sd, log_h_sd);
  }

  return log_lik;
}

double theta_log_likelihood_nc(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const arma::vec& y,
    const double h0,  // TODO implement
    const arma::vec& ht,
    const arma::vec& exp_h_half) {
  const int n = y.size();
  const double sigma = std::sqrt(sigma2),
               rho_const = std::sqrt(1 - rho * rho),
               log_rho_const = .5 * std::log(1 - rho * rho),
               B0 = std::sqrt(sigma2 / (1 - std::pow(phi, 2)));  // TODO stationary
  double log_lik = logdnorm(h0, mu, B0);
  for (int i = 0; i < n; i++) {
    double h_mean, y_mean, y_sd, log_y_sd;
    const double log_h_sd = 0,
                 h_sd = 1;
    if (i == 0) {
      h_mean = phi * (h0 - mu);  // because h0 is centered
    } else {
      h_mean = phi * ht[i-1];
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
    log_lik += logdnorm2(y[i], y_mean, y_sd, log_y_sd) + logdnorm2(ht[i], h_mean, h_sd, log_h_sd);
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
  return logdnorm2(mu, prior_mu[0], prior_mu[1]) +
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
  return logdbeta(.5 * (phi + 1.), prior_phi[0], prior_phi[1]) +
    logdbeta(.5 * (rho + 1.), prior_rho[0], prior_rho[1]) +
    (gammaprior ?
     logdgamma(sigma2, prior_sigma2[0], prior_sigma2[1]) :
     // moment matched InvGamma
     (-2. * std::log(sigma2) + logdinvgamma(1 / sigma2, prior_sigma2[0] + 2, prior_sigma2[1] / (prior_sigma2[0] * (prior_sigma2[0] + 1.)))));  // TODO I think this has been changed too fast; should be simply logdinvgamma(sigma2, ...)
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
  return 2 * (log4 + f + r + 0.5 * s - std::log((std::exp(2. * f) + 1.) * (std::exp(2. * r) + 1.)));
}

double theta_transform_inv_log_det_jac(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu) {
  return -(std::log((1. - phi * phi) * (1. - rho * rho) * sigma2));
}

template<unsigned int size>
arma::vec::fixed<size> rnorm_arma () {
  arma::vec::fixed<size> array;
  array.imbue([]() -> double { return R::norm_rand(); });
  return array;
}

arma::vec6 theta_propose_rwmh(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const ProposalDiffusionKen& diffusion_ken) {
  const arma::vec4 theta_old_t = theta_transform_inv(phi, rho, sigma2, mu);

  const arma::vec4 &proposal_mean_old = theta_old_t;
  const arma::vec4 theta_new_t_standardized = rnorm_arma<4>();
  const arma::vec4 theta_new_t =
    std::sqrt(diffusion_ken.get_scale()) * diffusion_ken.get_covariance_chol() * theta_new_t_standardized +
    proposal_mean_old;
  const arma::vec4 theta_new = theta_transform(theta_new_t[0], theta_new_t[1], theta_new_t[2], theta_new_t[3]);
  const double phi_new = theta_new[0], rho_new = theta_new[1], sigma2_new = theta_new[2], mu_new = theta_new[3];
  double theta_density_new = theta_transform_inv_log_det_jac(phi_new, rho_new, sigma2_new, mu_new);
  for (int i = 0; i < 4; i++) {
    theta_density_new += -.5 * std::pow(theta_new_t_standardized[i], 2);  // dnorm(x, 0, 1, true)
  }

  const arma::vec4 &proposal_mean_new = theta_new_t;
  const arma::vec4 theta_old_t_standardized =
    1 / std::sqrt(diffusion_ken.get_scale()) * diffusion_ken.get_covariance_chol_inv() *
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
    const double h0,  // TODO implement
    const arma::vec& h,
    const arma::vec2& prior_phi,
    const arma::vec2& prior_rho,
    const arma::vec2& prior_sigma2,
    const arma::vec2& prior_mu,
    const ProposalDiffusionKen& diffusion_ken) {
  const arma::vec4 theta_old_t = theta_transform_inv(phi, rho, sigma2, mu);

  const arma::vec4 grad_old = grad_theta_log_posterior(phi, rho, sigma2, mu, y, h0, h, prior_phi, prior_rho, prior_sigma2, prior_mu) %
    arma::vec4({1 - std::pow(phi, 2), 1 - std::pow(rho, 2), sigma2, 1});
  const arma::vec4 proposal_mean_old =
    theta_old_t +
    diffusion_ken.get_scale() * diffusion_ken.get_covariance() * grad_old;
  const arma::vec4 theta_new_t =
    std::sqrt(2.) * std::sqrt(diffusion_ken.get_scale()) * diffusion_ken.get_covariance_chol() * rnorm_arma<4>() +
    proposal_mean_old;
  const arma::vec4 theta_new = theta_transform(theta_new_t(0), theta_new_t(1), theta_new_t(2), theta_new_t(3));
  const double phi_new = theta_new(0), rho_new = theta_new(1), sigma2_new = theta_new(2), mu_new = theta_new(3);
  const arma::vec4 grad_new = grad_theta_log_posterior(phi_new, rho_new, sigma2_new, mu_new, y, h0, h, prior_phi, prior_rho, prior_sigma2, prior_mu) %
    arma::vec4({1 - std::pow(phi_new, 2), 1 - std::pow(rho_new, 2), sigma2_new, 1});

  double theta_density_new = theta_transform_inv_log_det_jac(phi_new, rho_new, sigma2_new, mu_new) +
    dmvnorm_mala(theta_new_t, theta_old_t, grad_old,
        diffusion_ken.get_covariance(), diffusion_ken.get_precision(), diffusion_ken.get_scale(), true);
  double theta_density_old = theta_transform_inv_log_det_jac(phi, rho, sigma2, mu) +
    dmvnorm_mala(theta_old_t, theta_new_t, grad_new,
        diffusion_ken.get_covariance(), diffusion_ken.get_precision(), diffusion_ken.get_scale(), true);

  return {phi_new, rho_new, sigma2_new, mu_new, theta_density_old, theta_density_new};
}

arma::vec thetamu_propose(
    const double phi,
    const double rho,
    const double sigma2,
    const ProposalDiffusionKen& adaptation_proposal) {
  const double mu = 0;
  const arma::vec3 theta_old_t = theta_transform_inv(phi, rho, sigma2, mu).head(3);

  const arma::vec3 &proposal_mean_old = theta_old_t;
  const arma::vec3 theta_new_t_standardized = rnorm_arma<3>();
  const arma::vec3 theta_new_t =
    std::sqrt(adaptation_proposal.get_scale()) * adaptation_proposal.get_covariance_chol() * theta_new_t_standardized +
    proposal_mean_old;
  const arma::vec3 theta_new = theta_transform(theta_new_t[0], theta_new_t[1], theta_new_t[2], mu).head(3);
  const double phi_new = theta_new[0], rho_new = theta_new[1], sigma2_new = theta_new[2];
  double theta_density_new = theta_transform_inv_log_det_jac(phi_new, rho_new, sigma2_new, mu);
  for (int i = 0; i < 3; i++) {
    theta_density_new += -.5 * std::pow(theta_new_t_standardized[i], 2);  // dnorm(x, 0, 1, true)
  }

  const arma::vec &proposal_mean_new = theta_new_t;
  const arma::vec theta_old_t_standardized =
    1 / std::sqrt(adaptation_proposal.get_scale()) * adaptation_proposal.get_covariance_chol_inv() *
    (theta_old_t - proposal_mean_new);
  double theta_density_old = theta_transform_inv_log_det_jac(phi, rho, sigma2, mu);
  for (int i = 0; i < 3; i++) {
    theta_density_old += -.5 * std::pow(theta_old_t_standardized[i], 2);  // dnorm(x, 0, 1, true)
  }

  return {phi_new, rho_new, sigma2_new, theta_density_old, theta_density_new};
}

arma::vec4 grad_theta_log_posterior(
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const arma::vec& y,
    const double h0,  // TODO implement
    const arma::vec& h,
    const arma::vec2& prior_phi,
    const arma::vec2& prior_rho,
    const arma::vec2& prior_sigma2,
    const arma::vec2& prior_mu) {
  const int n = y.n_elem;
  const double sigma = sqrt(sigma2);
  // Grad of log likelihood
  double d_phi = 0,
         d_rho = 0,
         d_sigma2 = 0,
         d_mu = 0;
  {  // p(h_0 | theta) and p(h_1 | h_0, theta)
    const double h_tilde = (h0 - mu) / sigma;
    d_phi += phi * (std::pow(h_tilde, 2) - 1 / (1 - std::pow(phi, 2)));
    d_rho += 0;
    d_sigma2 += .5 / sigma2 * ((1 - std::pow(phi, 2)) * std::pow(h_tilde, 2) - 1);
    d_mu += (1 - std::pow(phi, 2)) * h_tilde / sigma;
    const double delta = (h[0] - mu) / sigma - phi * h_tilde;
    d_phi += h_tilde * delta;
    d_rho += 0;
    d_sigma2 += .5 / sigma2 * (std::pow(delta, 2) - 1);
    d_mu += (1 - phi) * delta / sigma;
  }
  // p(h_{t+1} | h_t, theta) and p(y_t | h_{t+1}, h_t, theta)
  const double rho2 = std::pow(rho, 2);
  const double rho_const = 1 / (1 - rho2);
  for (int t = 0; t < n-1; t++) {
    const double y_h_const = y[t] * std::exp(-h[t] / 2),
                 h_tilde = (h[t] - mu) / sigma,
                 delta = (h[t+1] - mu) / sigma - phi * h_tilde,
                 delta2 = std::pow(delta, 2),
                 yhrd = -y_h_const + rho * delta;
    d_phi += h_tilde * delta;
    d_rho += 0;
    d_sigma2 += .5 * (delta2 - 1) / sigma2;
    d_mu += (1 - phi) * delta / sigma;
    d_phi += rho_const * h_tilde * rho * yhrd;
    d_rho += std::pow(rho_const, 2) * (-rho * std::pow(y_h_const, 2) + (1 + rho2) * y_h_const * delta - rho * delta2) - rho * rho_const;
    d_sigma2 += .5 * rho_const * (rho * delta * yhrd) / sigma2;
    d_mu += rho_const * (1 - phi) * rho * yhrd / sigma;
  }
  // Grad of log prior
  const double phi_beta = .5 * (phi + 1), rho_beta = .5 * (rho + 1);
  d_phi += .5 * ((prior_phi[0] - 1) / phi_beta - (prior_phi[1] - 1) / (1 - phi_beta));
  d_rho += .5 * ((prior_rho[0] - 1) / rho_beta - (prior_rho[1] - 1) / (1 - rho_beta));
  d_sigma2 += (prior_sigma2[0] - 1) / sigma2 - prior_sigma2[1];
  d_mu += -(mu - prior_mu[0]) / std::pow(prior_mu[1], 2);
  return {d_phi, d_rho, d_sigma2, d_mu};
}

}

