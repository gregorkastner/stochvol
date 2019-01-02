#include <RcppArmadillo.h>
#include <string>
#include <cmath>
#include "theta-utils.h"
#include "aug-kalman-filter.h"

using namespace Rcpp;

double theta_log_likelihood(const double phi, const double rho,
                            const double sigma2, const double mu,
                            const NumericVector y, const NumericVector h,
                            const CharacterVector centering) {
  double result;
  std::string scentering = as<std::string>(centering);
  if (scentering == "centered") {
    result = theta_log_likelihood_c(phi, rho, sigma2, mu, y, h);
  } else if (scentering == "non-centered") {
    result = theta_log_likelihood_nc(phi, rho, sigma2, mu, y, h);
  } else {
    ::Rf_error("invalid centering");
  }
  return result;
}

double theta_log_likelihood_c(const double phi, const double rho,
                              const double sigma2, const double mu,
                              const Rcpp::NumericVector y, const Rcpp::NumericVector h) {
  const int n = y.size();
  const double sigma = sqrt(sigma2);
  double log_lik = 0;
  for (int i = 0; i < n; i++) {
    double h_mean, h_sd, y_mean, y_sd;
    if (i == 0) {
      h_mean = mu;
      h_sd = sigma/sqrt(1-phi*phi);
    } else {
      h_mean = mu+phi*(h(i-1)-mu) + rho*sigma*exp(-h(i-1)/2)*y(i-1);
      h_sd = sigma*sqrt(1-rho*rho);
    }
    y_mean = 0;
    y_sd = exp(h(i)/2);
    log_lik += R::dnorm(y(i), y_mean, y_sd, true) + R::dnorm(h(i), h_mean, h_sd, true);
  }
  return log_lik;
}

double theta_log_likelihood_nc(const double phi, const double rho,
                               const double sigma2, const double mu,
                               const Rcpp::NumericVector y, const Rcpp::NumericVector h) {
  const int n = y.size();
  const double sigma = sqrt(sigma2);
  double log_lik = 0;
  for (int i = 0; i < n; i++) {
    double h_mean, h_sd, y_mean, y_sd;
    if (i == 0) {
      h_mean = 0;
      h_sd = 1/sqrt(1-phi*phi);
    } else {
      h_mean = phi*h(i-1);
      h_sd = 1;
    }
    if (i < n-1) {
      y_mean = exp((sigma*h(i)+mu)/2)*rho*(h(i+1)-phi*h(i));
      y_sd = exp((sigma*h(i)+mu)/2)*sqrt((1-rho*rho));
    } else {
      y_mean = 0;
      y_sd = exp((sigma*h(i) + mu)/2);
    }
    log_lik += R::dnorm(y(i), y_mean, y_sd, true) + R::dnorm(h(i), h_mean, h_sd, true);
  }
  return log_lik;
}

double auxtheta_log_prior(const double phi, const double rho,
                          const double sigma2,
                          const NumericVector prior_phi,
                          const NumericVector prior_rho,
                          const NumericVector prior_sigma2) {
  return (-log(2)) + R::dbeta((phi+1)/2, prior_phi(0), prior_phi(1), true) +
    (-log(2)) + R::dbeta((rho+1)/2, prior_rho(0), prior_rho(1), true) +
    R::dgamma(sigma2, prior_sigma2(0), 1/prior_sigma2(1), true);
}

double theta_log_prior(const double phi, const double rho,
                       const double sigma2, const double mu,
                       const NumericVector prior_phi,
                       const NumericVector prior_rho,
                       const NumericVector prior_sigma2,
                       const NumericVector prior_mu) {
  return auxtheta_log_prior(phi, rho, sigma2, prior_phi, prior_rho, prior_sigma2) +
    R::dnorm(mu, prior_mu(0), prior_mu(1), true);
}

NumericVector theta_transform(const double f, const double r,
                              const double s, const double m) {
  return NumericVector::create(1-2/(exp(2*f)+1), 1-2/(exp(2*r)+1), exp(s), m);
}

NumericVector theta_transform_inv(const double phi, const double rho,
                                  const double sigma2, const double mu) {
  return NumericVector::create(0.5*log(2/(1-phi)-1), 0.5*log(2/(1-rho)-1), log(sigma2), mu);
}

double theta_transform_log_det_jac(const double f, const double r,
                                   const double s, const double m) {
  return 2*(log(4) + f+r+s/2 - log(exp(2*f)+1) - log(exp(2*r)+1));
}

double theta_transform_inv_log_det_jac(const double phi, const double rho,
                                       const double sigma2, const double mu) {
  return -(log(1-phi*phi)+log(1-rho*rho)+log(sigma2));
}

NumericVector theta_proposal_stdev(const double phi, const double rho,
                                   const double sigma2, const double mu,
                                   const NumericVector y, const NumericVector h,
                                   const double stdev) {
  return NumericVector::create(1, 1/7, 2, 0.2)*0+stdev;
}

// source: http://gallery.rcpp.org/articles/dmvnorm_arma/
double dmvnrm(const arma::vec4 in,
              const arma::vec4 mean,
              const arma::mat44 sigma_chol_inv,  // inverse(chol(sigma)), lower tri
              const bool logd = false) {
  double out;

  const double constants = -0.5 * std::log(2.0 * M_PI);
  const double det_inv_root_sigma = arma::sum(arma::log(sigma_chol_inv.diag()));
  const arma::vec z = sigma_chol_inv * (in - mean);

  out = constants - 0.5*arma::sum(z%z) + det_inv_root_sigma;
  if (!logd) {
    out = std::exp(out);
  }

  return out;
}

arma::vec4 rwmh_logvars(const double phi, const double rho,
                        const double sigma2, const double mu,
                        const int T) {
  const arma::vec4 logVars = {-0.340073 - 0.000643900*T - 0.187647*log(1 - phi) - 1.528162*phi - 3.35185*sqrt(sigma2),
                              -0.787279 - 0.000691478*T + 0.170912*log(1 - phi) - 0.574605*phi - 3.97064*sqrt(sigma2),
                               3.388768 - 0.000817538*T + 0.499564*log(1 - phi) - 0.589691*phi - 6.62991*sqrt(sigma2),
                              -5.588485 - 0.000848774*T - 1.479299*log(1 - phi) - 2.430386*phi + 4.57432*sqrt(sigma2)};
  return logVars;
}

arma::vec4 rnorm4_arma() {
  const NumericVector res_rcpp = Rcpp::rnorm(4);
  return {res_rcpp(0), res_rcpp(1), res_rcpp(2), res_rcpp(3)};
}

NumericVector theta_propose(const double phi, const double rho,
                            const double sigma2, const double mu,
                            const NumericVector y, const NumericVector h,
                            const double stdev) {  // a factor for the covariance matrix
  const arma::mat44 cormat_chol = {{1, 0, 0, 0}, {0.01224282, 0.9999251, 0, 0}, {-0.54829, 0.0504385, 0.8347659, 0}, {0.05624419, -0.02407009, -0.1454298, 0.9874753}};

  const NumericVector theta_old_t_rcpp = theta_transform_inv(phi, rho, sigma2, mu);
  const arma::vec4 theta_old_t = {theta_old_t_rcpp(0), theta_old_t_rcpp(1), theta_old_t_rcpp(2), theta_old_t_rcpp(3)};

  const arma::vec4 &proposal_mean_old = theta_old_t;
  const arma::vec4 logvars_old = rwmh_logvars(phi, rho, sigma2, mu, h.size());
  const arma::mat44 sigma_chol_old = arma::diagmat(stdev*arma::exp(0.5*logvars_old)) * cormat_chol;
  const arma::mat44 sigma_chol_inv_old = arma::inv(arma::trimatl(sigma_chol_old));
  const arma::vec4 theta_new_t = sigma_chol_old*rnorm4_arma() + proposal_mean_old;
  const NumericVector theta_new = theta_transform(theta_new_t(0), theta_new_t(1), theta_new_t(2), theta_new_t(3));
  const double phi_new = theta_new(0), rho_new = theta_new(1), sigma2_new = theta_new(2), mu_new = theta_new(3);
  double theta_density_new = theta_transform_inv_log_det_jac(phi_new, rho_new, sigma2_new, mu_new);
  theta_density_new += dmvnrm(theta_new_t, proposal_mean_old, sigma_chol_inv_old, true);

  const arma::vec4 &proposal_mean_new = theta_new_t;
  const arma::vec4 logvars_new = rwmh_logvars(phi_new, rho_new, sigma2_new, mu_new, h.size());
  const arma::mat44 sigma_chol_new = arma::diagmat(stdev*arma::exp(0.5*logvars_new)) * cormat_chol;
  const arma::mat44 sigma_chol_inv_new = arma::inv(arma::trimatl(sigma_chol_new));
  double theta_density_old = theta_transform_inv_log_det_jac(phi, rho, sigma2, mu);
  theta_density_old += dmvnrm(theta_old_t, proposal_mean_new, sigma_chol_inv_new, true);

  return {phi_new, rho_new, sigma2_new, mu_new, theta_density_old, theta_density_new};
}

NumericVector theta_propose2(const double phi, const double rho,
                             const double sigma2, const double mu,
                             const NumericVector y, const NumericVector h,
                             const double stdev) {
  const NumericVector theta_old_t = theta_transform_inv(phi, rho, sigma2, mu);
  
  const NumericVector &proposal_mean_old = theta_old_t;
  const NumericVector proposal_sd_old = theta_proposal_stdev(phi, rho, sigma2, mu, y, h, stdev);
  const NumericVector theta_new_t = rnorm(4)*proposal_sd_old + proposal_mean_old;
  const NumericVector theta_new = theta_transform(theta_new_t(0), theta_new_t(1), theta_new_t(2), theta_new_t(3));
  const double phi_new = theta_new(0), rho_new = theta_new(1), sigma2_new = theta_new(2), mu_new = theta_new(3);
  double theta_density_new = theta_transform_inv_log_det_jac(phi_new, rho_new, sigma2_new, mu_new);
  for (int i = 0; i < 4; i++) {
    theta_density_new += R::dnorm(theta_new_t(i), proposal_mean_old(i), proposal_sd_old(i), true);
  }
  
  const NumericVector &proposal_mean_new = theta_new_t;
  const NumericVector proposal_sd_new = theta_proposal_stdev(phi_new, rho_new, sigma2_new, mu_new, y, h, stdev);
  double theta_density_old = theta_transform_inv_log_det_jac(phi, rho, sigma2, mu);
  for (int i = 0; i < 4; i++) {
    theta_density_old += R::dnorm(theta_old_t(i), proposal_mean_new(i), proposal_sd_new(i), true);
  }
  
  return NumericVector::create(phi_new, rho_new, sigma2_new, mu_new, theta_density_old, theta_density_new);
}

double auxtheta_log_posterior(const NumericVector auxtheta,
                              const NumericVector a,
                              const NumericVector b,
                              const NumericVector mm,
                              const NumericVector v,
                              const NumericVector d,
                              const NumericVector y_star,
                              const NumericVector prior_phi,
                              const NumericVector prior_rho,
                              const NumericVector prior_sigma2,
                              const NumericVector prior_mu) {
  //print(auxtheta);
  const double phi = auxtheta(0), rho = auxtheta(1), sigma2 = auxtheta(2);
  if (phi <= -1 || phi >= 1 || rho <= -1 || rho >= 1 || sigma2 <= 0) {
    return -1.0e100;
  }
  const NumericVector& m = mm;
  const List filter_results = aug_kalman_filter_c(phi, rho, sigma2,
                                                  a, b, m, v, d, y_star,
                                                  prior_mu(0), pow(prior_mu(1), 2));
  double likelihood = sum(log(as<NumericVector>(filter_results["D"]))) +
    log(as<double>(filter_results["Q"])) +
    sum(pow(as<NumericVector>(filter_results["f"]), 2) / as<NumericVector>(filter_results["D"])) -
    pow(as<double>(filter_results["q"]), 2) / as<double>(filter_results["Q"]);
  likelihood = -0.5 * likelihood;
  const double prior = auxtheta_log_prior(auxtheta(0), auxtheta(1), auxtheta(2),
                                          prior_phi, prior_rho, prior_sigma2);
  return likelihood+prior;
}
