/*
 * R package stochvol by
 *     Gregor Kastner Copyright (C) 2013-2018
 *     Gregor Kastner and Darjus Hosszejni Copyright (C) 2019-
 *  
 *  This file is part of the R package stochvol: Efficient Bayesian
 *  Inference for Stochastic Volatility Models.
 *  
 *  The R package stochvol is free software: you can redistribute it
 *  and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation, either version 2 or
 *  any later version of the License.
 *  
 *  The R package stochvol is distributed in the hope that it will be
 *  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 *  General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with the R package stochvol. If that is not the case, please
 *  refer to <http://www.gnu.org/licenses/>.
 */

/*
 * utils_parameters.cc
 * 
 * Definitions of the functions declared in utils_parameters.h.
 * Documentation can also be found in utils_parameters.h.
 */

#include <RcppArmadillo.h>
#include <expert.hpp>
#include <cmath>
#include "utils.h"
#include "utils_parameters.h"
#include "utils_latent_states.h"
#include "densities.h"
#include <array>
#include <utility>

using namespace Rcpp;

namespace stochvol {

double newton_raphson(
    const double startval,
    const double sum_tau,
    const int n,
    const double lambda,
    const double tol,
    const int maxiter) {
  double x = startval;
  double error = R_PosInf;
  double xnew;
  bool converged = false;
  for (int i = 0; not converged and i < maxiter; i++) {
    xnew = x - dlogdnu(x, sum_tau, lambda, n) / ddlogdnu(x, n);
    error = std::abs(xnew - x);
    x = xnew;
    converged = error < tol;
  }
  if (not converged) x = NA_REAL;
  return x;
}

namespace fast_sv {

namespace centered {

std::array<double, 9> regression_aggregates(
    const double h0,
    const arma::vec& h,
    const ExpertSpec_FastSV& expert) {
  const int T = h.size();

  double sum1 = h[0];
  double sum3 = h0*h[0];
  double sum4 = h[0]*h[0];
  for (int j = 1; j < T-1; j++) {
    sum1 += h[j];
    sum3 += h[j-1]*h[j];
    sum4 += h[j]*h[j];
  }
  double sum2 = sum1 + h[T-1];  // h_1 + h_2 + ... + h_T
  sum1 += h0;            // h_0 + h_1 + ... + h_{T-1}
  sum3 += h[T-2]*h[T-1]; // h_0*h_1 + h_1*h_2 + ... + h_{T-1}*h_T
  sum4 += h0*h0;         // h_0^2 + h_1^2 + ... + h_{T-1}^2

  const double tmp1 = 1/(((sum4 + expert.proposal_intercept_varinv) * (T + expert.proposal_phi_varinv) - std::pow(sum1, 2))),
               BT11 = (T + expert.proposal_phi_varinv)*tmp1,
               BT12 = -sum1*tmp1,
               BT22 = (sum4+expert.proposal_intercept_varinv)*tmp1,
               bT1 = BT11*sum3 + BT12*sum2,
               bT2 = BT12*sum3 + BT22*sum2;

  return {sum1, sum2, sum3, sum4, BT11, BT12, BT22, bT1, bT2};
}

std::pair<bool, double> sample_sigma(
    const double mu,
    const double phi,
    const double sigma,
    const double h0,
    const arma::vec& h,
    const PriorSpec& prior_spec,
    const ExpertSpec_FastSV& expert) {
  const int T = h.n_elem;

  const double Bh0inv = determine_Bh0inv(phi, prior_spec);
  const double cT = determine_cT(T, prior_spec, expert);
  const double C0 = determine_C0(prior_spec);

  double sigma_new = sigma;
  bool sigma_accepted = false;

  double z = std::pow(((h[0]-mu)-phi*(h0-mu)), 2);  // TODO: more efficiently via sum1, sum2, etc.
  for (int j = 0; j < (T-1); j++) {
    z += std::pow((h[j+1]-mu)-phi*(h[j]-mu), 2);
  }
  z += std::pow(h0-mu, 2) * Bh0inv;
  switch (expert.proposal_sigma2) {
    case ExpertSpec_FastSV::ProposalSigma2::LOG_RANDOM_WALK: {
      const double sigma2_prop = std::exp(R::rnorm(2 * std::log(sigma), expert.proposal_sigma2_rw_scale)),
            log_ar_prob = logacceptrateRW(sigma2_prop, std::pow(sigma, 2), 0.5 / prior_spec.sigma2.gamma.rate, T, z);

      if (log_ar_prob >= 0 or
          R::unif_rand() < std::exp(log_ar_prob)) {
        sigma_new = std::sqrt(sigma2_prop);
        sigma_accepted = true;
      }
    }
    break;
    case ExpertSpec_FastSV::ProposalSigma2::INDEPENDENCE: {
      switch (prior_spec.sigma2.distribution) {
        case PriorSpec::Sigma2::GAMMA: {
          const double CT = .5*z,
                       sigma2_prop = 1/R::rgamma(cT, 1/CT);
          if (std::log(R::unif_rand()) <
              logacceptrateGamma(sigma2_prop, std::pow(sigma, 2), 0.5 / prior_spec.sigma2.gamma.rate)) {
            sigma_new = std::sqrt(sigma2_prop);
            sigma_accepted = true;
          }
        }
        break;
        case PriorSpec::Sigma2::INVERSE_GAMMA: {
          const double CT = C0+.5*z;
          sigma_new = std::sqrt(1/R::rgamma(cT, 1/CT));
          sigma_accepted = true;
        }
        break;
        default:
        ::Rf_error("Constant prior for sigma not implemented in fast sv.");
      }
    }
  }

  return {sigma_accepted, sigma_new};
}

// propose gamma and phi given sigma
std::pair<bool, std::array<double, 2> > propose_beta(
    const double sigma,
    const std::array<double, 9>& agg,
    const ExpertSpec_FastSV& expert) {
  const double BT11 = agg[4],
               BT12 = agg[5],
               BT22 = agg[6],
               bT1 = agg[7],
               bT2 = agg[8];

  double chol11 = std::sqrt(BT11),
         chol12 = (BT12/chol11),
         chol22 = std::sqrt(BT22-chol12*chol12);
  chol11 *= sigma;
  chol12 *= sigma;
  chol22 *= sigma;

  double phi_prop,
         gamma_prop;
  switch (expert.proposal_phi) {
    case ExpertSpec_FastSV::ProposalPhi::REPEATED_ACCEPT_REJECT_NORMAL: { // draw from truncated normal via inversion method
      const double quant_low = R::pnorm(-1, bT1, chol11, true, false),
                   quant_high = R::pnorm(1, bT1, chol11, true, false);
      phi_prop = R::qnorm((quant_low + R::unif_rand()*(quant_high-quant_low)),
          bT1, chol11, true, false);
      gamma_prop = R::rnorm(bT2 + chol12*((phi_prop-bT1)/chol11),
          chol22);
    }
    break;
    case ExpertSpec_FastSV::ProposalPhi::IMMEDIATE_ACCEPT_REJECT_NORMAL: {  // draw from normal and reject (return) if outside
      const double innov = R::norm_rand();
      phi_prop = bT1 + chol11*innov;
      if (phi_prop >= 1 or
          phi_prop <= -1) { // outside the unit sphere
        return {false, {NAN, NAN}};
      }
      gamma_prop = bT2 + chol12*innov + chol22*R::norm_rand();
    }
    break;
  }

  return {true, {phi_prop, gamma_prop}};
}

// log acceptance rate for gamma and phi given sigma
double acceptance_rate_beta(
    const double mu,
    const double phi,
    const double sigma,
    const double gamma_prop,
    const double phi_prop,
    const double sigma_prop,
    const double h0,
    const PriorSpec& prior_spec,
    const ExpertSpec_FastSV& expert) {
  const double phi_prop_const = 1-phi_prop,  // some temps used for acceptance probability
               phi_const = 1-phi,
               proposal_intercept_sd = 1 / std::sqrt(expert.proposal_intercept_varinv);

  double ar_prob = 0;
  if (prior_spec.latent0.variance == PriorSpec::Latent0::STATIONARY) {
    ar_prob += logdnorm(h0, gamma_prop/phi_prop_const, sigma_prop/std::sqrt(1-phi_prop*phi_prop)) -
      logdnorm(h0, mu, sigma/std::sqrt(1-phi*phi));
  } else {
    ar_prob += logdnorm(h0, gamma_prop/phi_prop_const, std::sqrt(prior_spec.latent0.constant.value)*sigma_prop) -
      logdnorm(h0, mu, std::sqrt(prior_spec.latent0.constant.value)*sigma);
  }
  ar_prob += logdnorm(gamma_prop, prior_spec.mu.normal.mean*phi_prop_const, prior_spec.mu.normal.sd*phi_prop_const) -
    logdnorm(mu*phi_const, prior_spec.mu.normal.mean*phi_const, prior_spec.mu.normal.sd*phi_const) +
    logdbeta((phi_prop+1)/2, prior_spec.phi.beta.alpha, prior_spec.phi.beta.beta) -
    logdbeta((phi+1)/2, prior_spec.phi.beta.alpha, prior_spec.phi.beta.beta) +
    logdnorm(phi, 0, sigma * proposal_intercept_sd) -
    logdnorm(phi_prop, 0, sigma_prop * proposal_intercept_sd) +
    logdnorm(mu*phi_const, 0, sigma * proposal_intercept_sd) -
    logdnorm(gamma_prop, 0, sigma_prop * proposal_intercept_sd);
  
  return ar_prob;
}

SampledTheta draw_theta_1block(
    const double mu,
    const double phi,
    const double sigma,
    const double h0,
    const arma::vec& h,
    const PriorSpec& prior_spec,
    const ExpertSpec_FastSV& expert) {
  R_assert(prior_spec.sigma2.distribution == PriorSpec::Sigma2::GAMMA, "The centered 1-block sampler is only implemented with a gamma prior for sigma^2");
  R_assert(prior_spec.mu.distribution == PriorSpec::Mu::NORMAL, "The centered 1-block sampler is only implemented with a normal prior for mu");

  const int T = h.n_elem;
  const double cT = determine_cT(T, prior_spec, expert);

  const auto agg = regression_aggregates(h0, h, expert);
  const double sum2 = agg[1],
               sum3 = agg[2],
               sum4 = agg[3],
               bT1 = agg[7],
               bT2 = agg[8];

  // propose sigma
  const double CT = .5*((sum4 - h0*h0 + h[T-1]*h[T-1]) - bT1*sum3 - bT2*sum2);
  const double sigma_prop = std::sqrt(1/R::rgamma(cT, 1/CT));

  const auto proposed_beta = propose_beta(sigma, agg, expert);
  if (not proposed_beta.first) {  // rejected phi immediately
    return {mu, phi, sigma, false, false, false};
  }
  const double phi_prop = proposed_beta.second[0],
               gamma_prop = proposed_beta.second[1];

  // log acceptance probability
  const double ar_prob = logacceptrateGamma(std::pow(sigma_prop, 2), std::pow(sigma, 2), 0.5 / prior_spec.sigma2.gamma.rate) +
    acceptance_rate_beta(mu, phi, sigma, gamma_prop, phi_prop, sigma_prop, h0, prior_spec, expert);

  // accept/reject
  if (std::log(R::unif_rand()) < ar_prob) {
    return {gamma_prop/(1-phi_prop), phi_prop, sigma_prop, true, true, true};
  } else {
    return {mu, phi, sigma, false, false, false};
  }
}

SampledTheta draw_theta_2block(
    const double mu,
    const double phi,
    const double sigma,
    const double h0,
    const arma::vec& h,
    const PriorSpec& prior_spec,
    const ExpertSpec_FastSV& expert) {
  R_assert(prior_spec.mu.distribution == PriorSpec::Mu::NORMAL, "The centered 2-block sampler is only implemented with a normal prior for mu");

  const auto sigma_draw = sample_sigma(mu, phi, sigma, h0, h, prior_spec, expert);
  const bool sigma_accepted = sigma_draw.first;
  const double sigma_new = sigma_draw.second;

  const auto proposed_beta = propose_beta(sigma_new, regression_aggregates(h0, h, expert), expert);
  if (not proposed_beta.first) {  // rejected phi immediately
    return {mu, phi, sigma_new, false, false, sigma_accepted};
  }
  const double phi_prop = proposed_beta.second[0],
               gamma_prop = proposed_beta.second[1];

  // log acceptance probability
  const double ar_prob = acceptance_rate_beta(mu, phi, sigma_new, gamma_prop, phi_prop, sigma_new, h0, prior_spec, expert);

  // accept/reject
  if (std::log(R::unif_rand()) < ar_prob) {
    return {gamma_prop/(1-phi_prop), phi_prop, sigma_new, true, true, sigma_accepted};
  } else {
    return {mu, phi, sigma_new, false, false, sigma_accepted};
  }
}

SampledTheta draw_theta_3block(
    const double mu,
    const double phi,
    const double sigma,
    const double h0,
    const arma::vec& h,
    const PriorSpec& prior_spec,
    const ExpertSpec_FastSV& expert) {
  const int T = h.n_elem;

  SampledTheta result;
  result.mu = mu;
  result.phi = phi;
  result.sigma = sigma;
  result.mu_accepted = false;
  result.phi_accepted = false;
  result.sigma_accepted = false;

  const auto sigma_draw = sample_sigma(mu, phi, sigma, h0, h, prior_spec, expert);
  result.sigma_accepted = sigma_draw.first;
  result.sigma = sigma_draw.second;

  const auto agg = regression_aggregates(h0, h, expert);
  const double sum1 = agg[0],
               sum2 = agg[1],
               sum3 = agg[2],
               sum4 = agg[3];

  {
    // sampling of phi from full conditional:
    const double gamma = (1-phi)*mu,  // = 0 if mu = 0
                 sigma_new = result.sigma,
                 BTsqrt = sigma_new/std::sqrt(sum4+expert.proposal_intercept_varinv),
                 bT = (sum3-gamma*sum1)/(sum4+expert.proposal_intercept_varinv),
                 phi_prop = R::rnorm(bT, BTsqrt);

    double ar_prob = 0;
    if (prior_spec.latent0.variance == PriorSpec::Latent0::STATIONARY) { // needed only if prior of h0 depends on phi
      ar_prob += logdnorm(h0, mu, sigma_new/std::sqrt(1-phi_prop*phi_prop)) -
        logdnorm(h0, mu, sigma_new/std::sqrt(1-phi*phi));
    } 
    ar_prob += logdbeta((phi_prop+1)/2, prior_spec.phi.beta.alpha, prior_spec.phi.beta.beta) -
      logdbeta((phi+1)/2, prior_spec.phi.beta.alpha, prior_spec.phi.beta.beta) +
      logdnorm(phi, 0, sigma_new/std::sqrt(expert.proposal_intercept_varinv)) -
      logdnorm(phi_prop, 0, sigma_new/std::sqrt(expert.proposal_intercept_varinv));

    if (std::log(R::unif_rand()) < ar_prob) {
      result.phi = phi_prop;
      result.phi_accepted = true;
    }
  }

  if (prior_spec.mu.distribution != PriorSpec::Mu::CONSTANT) {
    // sampling of gamma from full conditional:
    const double phi_new = result.phi,
                 sigma_new = result.sigma,
                 phi_const = 1 - phi_new,
                 gamma = phi_const * mu,
                 BTsqrt = sigma_new/std::sqrt(T+expert.proposal_phi_varinv),
                 bT = (sum2-phi_new*sum1)/(T+expert.proposal_phi_varinv),
                 gamma_prop = R::rnorm(bT, BTsqrt);

    const double h0_sd = std::pow(determine_Bh0inv(phi_new, prior_spec), -0.5);
    const double ar_prob = logdnorm(h0, gamma_prop/phi_const, sigma_new * h0_sd) -
      logdnorm(h0, gamma/phi_const, sigma_new * h0_sd) +
      logdnorm(gamma_prop, prior_spec.mu.normal.mean*phi_const, prior_spec.mu.normal.sd*phi_const) -
      logdnorm(gamma, prior_spec.mu.normal.mean*phi_const, prior_spec.mu.normal.sd*phi_const) +
      logdnorm(gamma, 0, sigma_new/std::sqrt(expert.proposal_phi_varinv)) -
      logdnorm(gamma_prop, 0, sigma_new/std::sqrt(expert.proposal_phi_varinv));

    if (std::log(R::unif_rand()) < ar_prob) {
      result.mu = gamma_prop / phi_const;
      result.mu_accepted = true;
    }
  }

  return result;
}

}  // END namespace centered

namespace noncentered {

std::pair<bool, double> sample_phi(
    const double phi,
    const double ht0,
    const arma::vec& ht,
    const PriorSpec& prior_spec,
    const ExpertSpec_FastSV& expert) {
  const int T = ht.size();

  double sum_covar = ht0*ht[0],
         sum_var = ht0*ht0;
  double phi_prop;
  for (int j = 0; j < T-1; j++) {
    sum_covar += ht[j]*ht[j+1];
    sum_var += ht[j]*ht[j];
  }
  const double mean = sum_covar/sum_var,
               sd = 1/std::sqrt(sum_var);

  // Proposal
  switch (expert.proposal_phi) {
    case ExpertSpec_FastSV::ProposalPhi::REPEATED_ACCEPT_REJECT_NORMAL: { // draw from truncated normal via inversion method
      const double quant_low = R::pnorm(-1, mean, sd, true, false),
                   quant_high = R::pnorm(1, mean, sd, true, false);
      phi_prop = R::qnorm((quant_low + R::unif_rand()*(quant_high-quant_low)),
          mean, sd, true, false);
    }
    break;
    case ExpertSpec_FastSV::ProposalPhi::IMMEDIATE_ACCEPT_REJECT_NORMAL: {  // draw from normal and reject (return) if outside
      phi_prop = R::rnorm(mean, sd);
      if ((phi_prop >= 1) or
          (phi_prop <= -1)) { // outside the unit sphere
        return {false, phi};
      }
    }
    break;
    default:
      ::Rf_error("sample_phi: Mistake in the switch-case");
  }

  // MH step, acceptance prob ar_prob
  double ar_prob = 1;
  if (prior_spec.latent0.variance == PriorSpec::Latent0::STATIONARY) {
      ar_prob *= std::exp(logdnorm(ht0, 0, 1/std::sqrt(1-phi_prop*phi_prop))
          - logdnorm(ht0, 0, 1/std::sqrt(1-phi*phi)));
  }
  ar_prob *= propBeta((phi_prop+1)/2, (phi+1)/2, prior_spec.phi.beta.alpha, prior_spec.phi.beta.beta);
  //         ^^note that factor 1/2 from transformation of densities cancels

  // accept/reject
  if (R::unif_rand() < ar_prob) {
    return {true, phi_prop};
  } else {
    return {false, phi};
  }
}

SampledTheta draw_theta_2block(
    const arma::vec& log_data2,
    const double mu,
    const double phi,
    const double sigma,
    const double ht0,
    const arma::vec& ht,
    const arma::uvec& r,
    const PriorSpec& prior_spec,
    const ExpertSpec_FastSV& expert) {
  R_assert(prior_spec.mu.distribution == PriorSpec::Mu::NORMAL, "The non-centered 2-block sampler is only implemented with a normal prior for mu");

  // Draw mu and sigma
  const int T = ht.n_elem;
  double BT11 = std::pow(prior_spec.mu.normal.sd, -2),
         BT12 = 0,
         BT22 = 2*prior_spec.sigma2.gamma.rate,
         bT1 = 0,
         bT2 = prior_spec.mu.normal.mean * BT11;

  for (int j = 0; j < T; j++) {
    const double tmp1 = mix_varinv[r[j]],
                 tmp2 = (log_data2[j]-mix_mean[r[j]])*tmp1,
                 tmp3 = ht[j]*tmp1;
    BT11 += tmp1;
    BT12 -= tmp3;
    BT22 += ht[j]*tmp3;
    bT1 += ht[j]*tmp2;
    bT2 += tmp2;
  }

  {
    const double det = BT11*BT22-BT12*BT12;
    BT11 /= det;
    BT12 /= det;
    BT22 /= det;
  }

  {
    const double bT1_old = bT1;
    bT1 = BT11*bT1_old + BT12*bT2;
    bT2 = BT12*bT1_old + BT22*bT2;
  }

  const double chol11 = std::sqrt(BT11),
               chol12 = (BT12/chol11),
               chol22 = std::sqrt(BT22-chol12*chol12);

  const double innov = R::norm_rand(),
               sigma_new = bT1 + chol11*innov,
               mu_new = bT2 + chol12*innov + chol22*R::norm_rand();

  // Draw phi
  const auto phi_draw = sample_phi(phi, ht0, ht, prior_spec, expert);

  return {mu_new, phi_draw.second, std::fabs(sigma_new), true, phi_draw.first, true};
}

SampledTheta draw_theta_3block(
    const arma::vec& log_data2,
    const double mu,
    const double phi,
    const double sigma,
    const double ht0,
    const arma::vec& ht,
    const arma::uvec& r,
    const PriorSpec& prior_spec,
    const ExpertSpec_FastSV& expert) {
  const int T = ht.n_elem;
  // Draw sigma
  double sigma_new;
  {
    double tmp1 = 0,
           tmp2 = 0;
    for (int j = 0; j < T; j++) {
      tmp1 += ht[j]*ht[j]*mix_varinv[r[j]];
      tmp2 += ht[j]*(log_data2[j]-mix_mean[r[j]]-mu)*mix_varinv[r[j]];
    }
    const double BT11 = 1/(tmp1+2*prior_spec.sigma2.gamma.rate),
                 bT1 = BT11*tmp2;
    sigma_new = R::rnorm(bT1, std::sqrt(BT11));
  }

  // Draw mu
  double mu_new = mu;
  if (prior_spec.mu.distribution != PriorSpec::Mu::CONSTANT) {
    double tmp1 = 0,
           tmp2 = 0;
    for (int j = 0; j < T; j++) {
      tmp1 += mix_varinv[r[j]];
      tmp2 += (log_data2[j]-mix_mean[r[j]]-sigma_new*ht[j])*mix_varinv[r[j]];
    }
    const double Bmu_inv = std::pow(prior_spec.mu.normal.sd, -2);
    const int BT22 = 1/(tmp1+Bmu_inv),
              bT2 = BT22*(tmp2 + prior_spec.mu.normal.mean * Bmu_inv);
    mu_new = R::rnorm(bT2, std::sqrt(BT22));
  }

  // Draw phi
  const auto phi_draw = sample_phi(phi, ht0, ht, prior_spec, expert);

  return {mu_new, phi_draw.second, std::fabs(sigma_new), prior_spec.mu.distribution != PriorSpec::Mu::CONSTANT, phi_draw.first, true};
}

}  // END namespace noncentered

}  // END namespace fast_sv

namespace general_sv {

double theta_log_likelihood_c(
    const arma::vec& data,
    const double mu,
    const double phi,
    const double sigma,
    const double rho,
    const double h0,
    const arma::vec& h,
    const arma::vec& exp_h_half,
    const PriorSpec& prior_spec) {
  const arma::vec& y = data;
  const int n = y.size();
  const double h_sd_t = sigma * std::sqrt(1 - std::pow(rho, 2)),
               log_h_sd_t = std::log(h_sd_t),
               B0 = sigma * std::pow(determine_Bh0inv(phi, prior_spec), -0.5);
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
    const arma::vec& data,
    const double mu,
    const double phi,
    const double sigma,
    const double rho,
    const double ht0,
    const arma::vec& ht,
    const arma::vec& exp_h_half,
    const PriorSpec& prior_spec) {
  const arma::vec& y = data;
  const int n = y.size();
  const double rho_const = std::sqrt(1 - std::pow(rho, 2)),
               log_rho_const = std::log(rho_const),
               B0 = std::pow(determine_Bh0inv(phi, prior_spec), -0.5);
  double log_lik = logdnorm(ht0, 0, B0);
  for (int i = 0; i < n; i++) {
    double h_mean, y_mean, y_sd, log_y_sd;
    const double log_h_sd = 0,
                 h_sd = 1;
    if (i == 0) {
      h_mean = phi * ht0;
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

double theta_log_likelihood(
    const arma::vec& data,
    const double mu,
    const double phi,
    const double sigma,
    const double rho,
    const double h0,
    const double ht0,
    const arma::vec& h,
    const arma::vec& ht,
    const arma::vec& exp_h_half,
    const PriorSpec& prior_spec,
    const Parameterization centering) {
  double result;
  switch (centering) {
    case Parameterization::CENTERED:
      result = theta_log_likelihood_c(data, mu, phi, sigma, rho, h0, h, exp_h_half, prior_spec);
      break;
    case Parameterization::NONCENTERED:
      result = theta_log_likelihood_nc(data, mu, phi, sigma, rho, ht0, ht, exp_h_half, prior_spec);
      break;
    default:
      ::Rf_error("theta_log_likelihood: Mistake in the switch-case");
  }
  return result;
}

double theta_log_prior(
    const double mu,
    const double phi,
    const double sigma,
    const double rho,
    const PriorSpec& prior_spec) {
  double mu_density;
  switch (prior_spec.mu.distribution) {
    case PriorSpec::Mu::CONSTANT:
      mu_density = 0;
      break;
    case PriorSpec::Mu::NORMAL:
      mu_density = logdnorm2(mu, prior_spec.mu.normal.mean, prior_spec.mu.normal.sd);
      break;
    default:
      ::Rf_error("theta_log_prior: Mistake in the switch-case");
  }
  double phi_density;
  switch (prior_spec.phi.distribution) {
    case PriorSpec::Phi::CONSTANT:
      phi_density = 0;
      break;
    case PriorSpec::Phi::BETA:
      phi_density = logdbeta(.5 * (phi + 1.), prior_spec.phi.beta.alpha, prior_spec.phi.beta.beta);
      break;
    case PriorSpec::Phi::NORMAL:
      phi_density = logdnorm2(phi, prior_spec.phi.normal.mean, prior_spec.phi.normal.sd);
      break;
    default:
      ::Rf_error("theta_log_prior: Mistake in the switch-case");
  }
  double sigma_density;
  switch (prior_spec.sigma2.distribution) {
    case PriorSpec::Sigma2::CONSTANT:
      sigma_density = 0;
      break;
    case PriorSpec::Sigma2::GAMMA:
      sigma_density = std::log(sigma) + logdgamma(std::pow(sigma, 2), prior_spec.sigma2.gamma.shape, prior_spec.sigma2.gamma.rate);
      break;
    case PriorSpec::Sigma2::INVERSE_GAMMA:
      sigma_density = std::log(sigma) + logdinvgamma(std::pow(sigma, 2), prior_spec.sigma2.inverse_gamma.shape, prior_spec.sigma2.inverse_gamma.scale);
      break;
    default:
      ::Rf_error("theta_log_prior: Mistake in the switch-case");
  }
  double rho_density;
  switch (prior_spec.rho.distribution) {
    case PriorSpec::Rho::CONSTANT:
      rho_density = 0;
      break;
    case PriorSpec::Rho::BETA:
      rho_density = logdbeta(.5 * (rho + 1.), prior_spec.rho.beta.alpha, prior_spec.rho.beta.beta);
      break;
    default:
      ::Rf_error("theta_log_prior: Mistake in the switch-case");
  }
  return mu_density + phi_density + sigma_density + rho_density;
}

arma::vec4 theta_transform(
    const double m,
    const double f,
    const double s,
    const double r,
    const PriorSpec& prior_spec) {
  return {
    m,
    prior_spec.phi.distribution == PriorSpec::Phi::BETA ? 1 - 2/(std::exp(2 * f) + 1) : f,
    std::exp(s),
    prior_spec.rho.distribution == PriorSpec::Rho::BETA ? 1 - 2/(std::exp(2 * r) + 1) : r};
}

arma::vec4 theta_transform_inv(
    const double mu,
    const double phi,
    const double sigma,
    const double rho,
    const PriorSpec& prior_spec) {
  return {
    mu,
    prior_spec.phi.distribution == PriorSpec::Phi::BETA ? 0.5 * std::log(2. / (1 - phi) - 1) : phi,
    std::log(sigma),
    prior_spec.rho.distribution == PriorSpec::Rho::BETA ? 0.5 * std::log(2. / (1 - rho) - 1) : rho};
}

double theta_transform_log_det_jac(
    const double m,
    const double f,
    const double s,
    const double r,
    const PriorSpec& prior_spec) {
  static const double log4 = std::log(4.);
  const double phi_part1 = prior_spec.phi.distribution == PriorSpec::Phi::BETA ? f : 0,
               phi_part2 = prior_spec.phi.distribution == PriorSpec::Phi::BETA ? std::exp(2. * f) + 1. : 1;
  const double rho_part1 = prior_spec.rho.distribution == PriorSpec::Rho::BETA ? r : 0,
               rho_part2 = prior_spec.rho.distribution == PriorSpec::Rho::BETA ? std::exp(2. * r) + 1. : 1;
  return 2 * (log4 + phi_part1 + 0.5 * s + rho_part1 - std::log(phi_part2 * rho_part2));
}

double theta_transform_inv_log_det_jac(
    const double mu,
    const double phi,
    const double sigma,
    const double rho,
    const PriorSpec& prior_spec) {
  const double phi_part = prior_spec.phi.distribution == PriorSpec::Phi::BETA ? 1. - std::pow(phi, 2) : 1,
               rho_part = prior_spec.rho.distribution == PriorSpec::Rho::BETA ? 1. - std::pow(rho, 2) : 1;
  return -std::log(phi_part * sigma * rho_part);
}

arma::vec rnorm_arma (const size_t size) {
  arma::vec array(size);
  array.imbue([]() -> double { return R::norm_rand(); });
  return array;
}

std::array<double, 6> theta_propose_rwmh(
    const double mu,
    const double phi,
    const double sigma,
    const double rho,
    const PriorSpec& prior_spec,
    const ProposalDiffusionKen& diffusion_ken,
    const arma::uvec& update_indicator) {
  const arma::uvec update_index = arma::find(update_indicator);
  const arma::vec4 theta_old_full {mu, phi, sigma, rho};
  const arma::vec4 theta_old_t_full = theta_transform_inv(mu, phi, sigma, rho, prior_spec);

  const arma::vec theta_old_t = theta_old_t_full.elem(update_index);

  const arma::vec &proposal_mean_old = theta_old_t;  // rename
  const arma::vec theta_new_t_standardized = rnorm_arma(update_index.n_elem);
  const arma::vec theta_new_t =
    std::sqrt(diffusion_ken.get_scale()) * diffusion_ken.get_covariance_chol() * theta_new_t_standardized +
    proposal_mean_old;
  arma::vec4 theta_new_t_full {theta_old_t_full};
  theta_new_t_full.elem(update_index) = theta_new_t;
  const arma::vec4 theta_new_full = update_indicator % theta_transform(theta_new_t_full[0], theta_new_t_full[1], theta_new_t_full[2], theta_new_t_full[3], prior_spec) +
    (1u - update_indicator) % theta_old_full;
  
  const double mu_new = theta_new_full[0],
               phi_new = theta_new_full[1],
               sigma_new = theta_new_full[2],
               rho_new = theta_new_full[3];

  // Proposal density for theta_new given theta_old
  const double theta_density_new = theta_transform_inv_log_det_jac(mu_new, phi_new, sigma_new, rho_new, prior_spec) +
    (-0.5) * arma::sum(arma::square(theta_new_t_standardized));

  // Proposal density for theta_old given theta_new
  const arma::vec &proposal_mean_new = theta_new_t;
  const arma::vec theta_old_t_standardized =
    1 / std::sqrt(diffusion_ken.get_scale()) * diffusion_ken.get_covariance_chol_inv() *
    (theta_old_t - proposal_mean_new);
  const double theta_density_old = theta_transform_inv_log_det_jac(mu, phi, sigma, rho, prior_spec) +
    (-0.5) * arma::sum(arma::square(theta_old_t_standardized));

  return {mu_new, phi_new, sigma_new, rho_new, theta_density_old, theta_density_new};
}

}  // END namespace general_sv

}

