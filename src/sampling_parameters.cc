/*
 * R package stochvol by
 *     Gregor Kastner Copyright (C) 2013-2020
 *     Darjus Hosszejni Copyright (C) 2019-2020
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
#include "sampling_parameters.h"
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
    case 1: [[fallthrough]];
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
    exp_h_half_proposal_nc = arma::exp(.5 * noncentered_to_centered(mu_prop, sigma_prop, ht));
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

}  // END namespace general_sv

}

