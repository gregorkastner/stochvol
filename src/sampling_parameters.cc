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
 * sampling_parameters.cc
 * 
 * Definitions of the functions declared in sampling_parameters.h.
 * Documentation can also be found in sampling_parameters.h.
 */

#include <RcppArmadillo.h>
#include <cmath>
#include <expert.hpp>
#include "sampling_latent_states.h"
#include "sampling_parameters.h"
#include "type_definitions.hpp"
#include "utils_parameters.h"
#include "utils_latent_states.h"
#include "densities.h"
#include "utils.h"

using namespace Rcpp;

namespace stochvol {

namespace fast_sv {

namespace centered {

// Fast SV step (b): sample mu, phi, sigma - __CENTERED__ version
SampledTheta regression(
    const double mu,
    const double phi,
    const double sigma,
    const double h0,
    const arma::vec& h,
    const PriorSpec& prior_spec,
    const ExpertSpec_FastSV& expert) {
  switch (expert.mh_blocking_steps) {
    case 1:
      return centered::draw_theta_1block(mu, phi, sigma, h0, h, prior_spec, expert);
    case 2:
      return centered::draw_theta_2block(mu, phi, sigma, h0, h, prior_spec, expert);
    case 3:
      return centered::draw_theta_3block(mu, phi, sigma, h0, h, prior_spec, expert);
    default:
      ::Rf_error("Parameter fast_sv$mh_blocking_steps should an integer between 1 and 3.");
  }
}

}  // END namespace centered

namespace noncentered {

// Fast SV step (b): sample mu, phi, sigma - __NONCENTERED__ version
SampledTheta regression(
    const arma::vec& log_data2,
    const double mu,
    const double phi,
    const double sigma,
    const double h0,
    const arma::vec& h,
    const arma::uvec& r,
    const PriorSpec& prior_spec,
    const ExpertSpec_FastSV& expert) {
  switch (expert.mh_blocking_steps) {
    case 1: //[[fallthrough]];
    case 2:
      return noncentered::draw_theta_2block(log_data2, mu, phi, sigma, h0, h, r, prior_spec, expert);
    case 3:
      return noncentered::draw_theta_3block(log_data2, mu, phi, sigma, h0, h, r, prior_spec, expert);
    default:
      ::Rf_error("Parameter fast_sv$mh_blocking_steps should an integer between 1 and 3.");
  }
}

}  // END namespace noncentered

// Main function for drawing the parameters.
// See documentation above the declaration
SampledTheta draw_theta(
    const arma::vec& log_data2,
    const double mu,
    const double phi,
    const double sigma,
    const double h0,
    const double ht0,
    const arma::vec& h,
    const arma::vec& ht,
    const arma::uvec& r,
    const PriorSpec& prior_spec,
    const ExpertSpec_FastSV& expert,
    const Parameterization parameterization) {
  switch (parameterization) {
    case Parameterization::CENTERED:
      return centered::regression(mu, phi, sigma, h0, h, prior_spec, expert);
    case Parameterization::NONCENTERED:
      return noncentered::regression(log_data2, mu, phi, sigma, ht0, ht, r, prior_spec, expert);
    default:
      ::Rf_error("draw_theta: Mistake in the switch-case");
  }
}

}  // END namespace fast_sv

namespace general_sv {

// Main function for drawing the parameters.
// See documentation above the declaration
SampledTheta draw_theta(
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
    arma::vec& exp_h_half_proposal_nc,
    const arma::uvec& update_indicator,
    const PriorSpec& prior_spec,
    const ExpertSpec_GeneralSV& expert,
    const ProposalDiffusionKen& diffusion_ken,
    const Parameterization parameterization) {
  const std::array<double, 6> proposed = theta_propose_rwmh(mu, phi, sigma, rho, prior_spec, diffusion_ken, update_indicator);
  const double mu_prop = proposed[0], phi_prop = proposed[1], sigma_prop = proposed[2],
    rho_prop = proposed[3], prop_old_logdens = proposed[4], prop_new_logdens = proposed[5];
  if (parameterization == Parameterization::NONCENTERED) {
    exp_h_half_proposal_nc = arma::exp(.5 * ::stochvol::noncentered_to_centered(mu_prop, sigma_prop, ht));
  }
  const arma::vec& exp_h_half_proposal = parameterization == Parameterization::CENTERED ? exp_h_half : exp_h_half_proposal_nc;

  const double log_acceptance =
    (theta_log_prior(mu_prop, phi_prop, sigma_prop, rho_prop, prior_spec) +
     theta_log_likelihood(data, mu_prop, phi_prop, sigma_prop, rho_prop, h0, ht0, h, ht, exp_h_half_proposal, prior_spec, parameterization)) -
    (theta_log_prior(mu, phi, sigma, rho, prior_spec) +
     theta_log_likelihood(data, mu, phi, sigma, rho, h0, ht0, h, ht, exp_h_half, prior_spec, parameterization)) -
    (prop_new_logdens - prop_old_logdens);

  const bool accepted = log_acceptance > 0 or
    std::exp(log_acceptance) > R::unif_rand();
  if (accepted) {
    return {mu_prop, phi_prop, sigma_prop, rho_prop, update_indicator[0] == 1u, update_indicator[1] == 1u, update_indicator[2] == 1u, update_indicator[3] == 1u};
  } else {
    return {mu, phi, sigma, rho, false, false, false, false};
  }
}

namespace centered {

// Main function for drawing the parameters.
// See documentation above the declaration
SampledTheta draw_theta(
    const double mu,
    const double phi,
    const double sigma,
    const double rho,
    const SufficientStatistic& sufficient_statistic,
    const arma::uvec& update_indicator,
    const PriorSpec& prior_spec,
    const ExpertSpec_GeneralSV& expert,
    const ProposalDiffusionKen& diffusion_ken) {
  const std::array<double, 6> proposed = theta_propose_rwmh(mu, phi, sigma, rho, prior_spec, diffusion_ken, update_indicator);
  const double mu_prop = proposed[0], phi_prop = proposed[1], sigma_prop = proposed[2],
    rho_prop = proposed[3], prop_old_logdens = proposed[4], prop_new_logdens = proposed[5];

  const double log_acceptance =
    (theta_log_prior(mu_prop, phi_prop, sigma_prop, rho_prop, prior_spec) +
     theta_log_likelihood(mu_prop, phi_prop, sigma_prop, rho_prop, sufficient_statistic, prior_spec)) -
    (theta_log_prior(mu, phi, sigma, rho, prior_spec) +
     theta_log_likelihood(mu, phi, sigma, rho, sufficient_statistic, prior_spec)) -
    (prop_new_logdens - prop_old_logdens);

  const bool accepted = log_acceptance > 0 or
    std::exp(log_acceptance) > R::unif_rand();
  if (accepted) {
    return {mu_prop, phi_prop, sigma_prop, rho_prop, update_indicator[0] == 1u, update_indicator[1] == 1u, update_indicator[2] == 1u, update_indicator[3] == 1u};
  } else {
    return {mu, phi, sigma, rho, false, false, false, false};
  }
}

#ifndef NDEBUG

double test_function (
    const double d1,
    const double d2,
    const double d3,
    const double h0,
    const double h1,
    const double h2,
    const double h3,
    const double mu1,
    const double mu2,
    const double phi1,
    const double phi2,
    const double sigma1,
    const double sigma2,
    const double rho1,
    const double rho2) {
  const PriorSpec prior_spec;
  const arma::vec data{d1, d2, d3};
  const arma::vec h{h1, h2, h3};
  const arma::vec exp_h_half = arma::exp(.5 * h);
  const auto sufficient_statistic2 = compute_sufficient_statistic(data, h0, h);
  const double hmmmmm = 
    theta_log_likelihood(mu1, phi1, sigma1, rho1, sufficient_statistic2, prior_spec) -
    theta_log_likelihood(mu2, phi2, sigma2, rho2, sufficient_statistic2, prior_spec);
  const double hnnnnn =
    ::stochvol::general_sv::theta_log_likelihood(data, mu1, phi1, sigma1, rho1, h0, h0, h, h, exp_h_half, prior_spec, Parameterization::CENTERED) -
    ::stochvol::general_sv::theta_log_likelihood(data, mu2, phi2, sigma2, rho2, h0, h0, h, h, exp_h_half, prior_spec, Parameterization::CENTERED);
  return hmmmmm - hnnnnn;
}

#endif

}  // END namespace centered

#ifndef NDEBUG

double test_function (
    const double d1,
    const double d2,
    const double d3,
    const double h0,
    const double h1,
    const double h2,
    const double h3,
    const double mu1,
    const double mu2,
    const double phi1,
    const double phi2,
    const double sigma1,
    const double sigma2,
    const double rho1,
    const double rho2) {
  const PriorSpec prior_spec;
  const arma::vec data{d1, d2, d3};
  const arma::vec h{h1, h2, h3};
  const arma::vec exp_h_half = arma::exp(.5 * h);
  const arma::vec c1 = centered_to_noncentered(mu1, phi1, sigma1, rho1, data, h0, h, prior_spec);
  const arma::vec c2 = centered_to_noncentered(mu2, phi2, sigma2, rho2, data, h0, h, prior_spec);
  // compute back
  const LatentVector h_back1 = noncentered_to_centered(mu1, phi1, sigma1, rho1, data, c1, prior_spec);
  const LatentVector h_back2 = noncentered_to_centered(mu2, phi2, sigma2, rho2, data, c2, prior_spec);
  const auto sufficient_statistic2 = centered::compute_sufficient_statistic(data, h0, h);
  const double hmmmmm = 
    centered::theta_log_likelihood(mu1, phi1, sigma1, rho1, sufficient_statistic2, prior_spec) -
    centered::theta_log_likelihood(mu2, phi2, sigma2, rho2, sufficient_statistic2, prior_spec);
  const double hnnnnn =
    noncentered::theta_log_likelihood(mu1, phi1, sigma1, rho1, {data, c1}, prior_spec) -
    noncentered::theta_log_likelihood(mu2, phi2, sigma2, rho2, {data, c2}, prior_spec);
  return hmmmmm - hnnnnn;
}

#endif

namespace noncentered {

// Main function for drawing the parameters.
// See documentation above the declaration
SampledTheta draw_theta(
    const double mu,
    const double phi,
    const double sigma,
    const double rho,
    const SufficientStatistic& sufficient_statistic,
    const arma::uvec& update_indicator,
    const PriorSpec& prior_spec,
    const ExpertSpec_GeneralSV& expert,
    const ProposalDiffusionKen& diffusion_ken) {
  const std::array<double, 6> proposed = theta_propose_rwmh(mu, phi, sigma, rho, prior_spec, diffusion_ken, update_indicator);
  const double mu_prop = proposed[0], phi_prop = proposed[1], sigma_prop = proposed[2],
    rho_prop = proposed[3], prop_old_logdens = proposed[4], prop_new_logdens = proposed[5];

  const double log_acceptance =
    (theta_log_prior(mu_prop, phi_prop, sigma_prop, rho_prop, prior_spec) +
     theta_log_likelihood(mu_prop, phi_prop, sigma_prop, rho_prop, sufficient_statistic, prior_spec)) -
    (theta_log_prior(mu, phi, sigma, rho, prior_spec) +
     theta_log_likelihood(mu, phi, sigma, rho, sufficient_statistic, prior_spec)) -
    (prop_new_logdens - prop_old_logdens);

  const bool accepted = log_acceptance > 0 or
    std::exp(log_acceptance) > R::unif_rand();
  if (accepted) {
    return {mu_prop, phi_prop, sigma_prop, rho_prop, update_indicator[0] == 1u, update_indicator[1] == 1u, update_indicator[2] == 1u, update_indicator[3] == 1u};
  } else {
    return {mu, phi, sigma, rho, false, false, false, false};
  }
}

}  // END namespace noncentered

}  // END namespace general_sv

}

