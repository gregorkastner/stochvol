#include <RcppArmadillo.h>
#include "utils_latent_states.h"
#include "auxmix.h"
#include <cmath>

double h_log_posterior(
    const arma::vec& h,  // centered
    const arma::vec& y,
    const double phi,
    const double rho,
    const double sigma2,
    const double mu) {
  const double sigma = std::sqrt(sigma2);
  const double sigma2_inv = 1 / sigma2;
  const double rho_const = std::sqrt(1 - rho * rho);
  const int n = y.size();
  const arma::vec exp_h_half = arma::exp(0.5 * h);
  double result = -.5 * std::pow(h[0] - mu, 2) * sigma2_inv * (1 - phi * phi);  // log p(h_1 | theta)
  for (int t = 0; t < n - 1; t++) {
    result += -.5 * std::pow(h[t + 1] - (mu + phi * (h[t] - mu)), 2) * sigma2_inv;
    result += -.5 * std::pow((y[t] - (exp_h_half[t] * rho * (h[t + 1] - mu - phi * (h[t] - mu)) / sigma)) / (exp_h_half[t] * rho_const), 2) - .5 * h[t];
  }
  result += -.5 * std::pow(y[n - 1] / exp_h_half[n - 1], 2) - .5 * h[n - 1];
  return result;
}

double h_aux_log_posterior(
    const arma::vec& h,  // centered
    const arma::vec& y_star,
    const arma::ivec& d,
    const double phi,
    const double rho,
    const double sigma2,
    const double mu) {
  const int n = y_star.size();
  static const int mix_count = mix_a.n_elem;
  const double sigma = std::sqrt(sigma2);
  static const arma::vec::fixed<10> exp_m_half = arma::exp(mix_mean * .5);
  const arma::vec C = rho * sigma * exp_m_half;  // re-used constant

  double result = -.5 * std::pow(h[0] - mu, 2) / sigma2 * (1 - phi * phi);  // log p(h_1 | theta)
  for (int t = 0; t < n; t++) {
    double subresult = 0;
    if (t < n - 1) {
      for (int j = 0; j < mix_count; j++) {
        const double h_mean = mu + phi * (h[t] - mu) + d[t] * C[j] * mix_a[j];
        const double h_var = mix_var[j] * std::pow(C[j] * mix_b[j], 2) + sigma2 * (1 - rho * rho);
        const double yh_cov = d[t] * C[j] * mix_b[j] * mix_var[j];
        const double y_mean = h[t] + mix_mean[j] + yh_cov / h_var * (h[t + 1] - h_mean);
        const double y_var = (1 - std::pow(yh_cov, 2) / (mix_var[j] * h_var)) * mix_var[j];
        subresult += std::exp(-.5 * (std::pow(y_star[t] - y_mean, 2) / y_var + std::pow(h[t + 1] - h_mean, 2) / h_var)) / std::sqrt(y_var * h_var) * mix_prob[j];
      }
    } else {
      for (int j = 0; j < mix_count; j++) {
        subresult += std::exp(-.5 * std::pow(y_star[t] - (h[t] + mix_mean[j]), 2) / mix_var[j]) / std::sqrt(mix_var[j]) * mix_prob[j];
      }
    }
    result += std::log(subresult);
  }

  return result;
}

