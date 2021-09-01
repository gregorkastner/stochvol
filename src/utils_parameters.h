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
 * utils_parameters.h
 * 
 * Utility functions and long implementations related to sampling_parameters.cc.
 * 
 * The functions are separated into the fast_sv and the general_sv
 * variants due to the varying expert options.
 */

#ifndef _UTILS_PARAMETERS_H_
#define _UTILS_PARAMETERS_H_

#include <RcppArmadillo.h>
#include <adaptation.hpp>
#include <expert.hpp>
#include <array>

namespace stochvol {

// Find the root of a function (Newton-Raphson)
double newton_raphson(
    const double startval,
    const double sum_tau,
    const int n,
    const double lambda,
    const double tol = 1e-03,
    const int maxiter = 50);

namespace fast_sv {

// Encapsulate the four parameters mu, phi, and sigma,
// as one returned object, including whether
// they were updated/accepted.
struct SampledTheta {
  double mu,
         phi,
         sigma;
  bool mu_accepted,
       phi_accepted,
       sigma_accepted;
};

namespace centered {

// Single-block update of mu, phi, and sigma -- centered.
// For more details, see Kastner and Frühwirth-Schnatter (2014)
SampledTheta draw_theta_1block(
    const double mu,
    const double phi,
    const double sigma,
    const double h0,
    const arma::vec& h,
    const PriorSpec& prior_spec,
    const ExpertSpec_FastSV& expert);

// Two-block update of mu, phi, and sigma -- centered.
// For more details, see Kastner and Frühwirth-Schnatter (2014)
SampledTheta draw_theta_2block(
    const double mu,
    const double phi,
    const double sigma,
    const double h0,
    const arma::vec& h,
    const PriorSpec& prior_spec,
    const ExpertSpec_FastSV& expert);

// Three-block update of mu, phi, and sigma -- centered.
// For more details, see Kastner and Frühwirth-Schnatter (2014)
SampledTheta draw_theta_3block(
    const double mu,
    const double phi,
    const double sigma,
    const double h0,
    const arma::vec& h,
    const PriorSpec& prior_spec,
    const ExpertSpec_FastSV& expert);

}  // END namespace centered

namespace noncentered {

// Two-block update of mu, phi, and sigma -- noncentered.
// For more details, see Kastner and Frühwirth-Schnatter (2014)
SampledTheta draw_theta_2block(
    const arma::vec& log_data2,
    const double mu,
    const double phi,
    const double sigma,
    const double ht0,
    const arma::vec& ht,
    const arma::uvec& r,
    const PriorSpec& prior_spec,
    const ExpertSpec_FastSV& expert);

// Three-block update of mu, phi, and sigma -- noncentered.
// For more details, see Kastner and Frühwirth-Schnatter (2014)
SampledTheta draw_theta_3block(
    const arma::vec& log_data2,
    const double mu,
    const double phi,
    const double sigma,
    const double ht0,
    const arma::vec& ht,
    const arma::uvec& r,
    const PriorSpec& prior_spec,
    const ExpertSpec_FastSV& expert);

}  // END namespace noncentered

// Determine c_T from the prior specification and the
// expert settings. This value is used in the update
// of sigma.
inline
double determine_cT(
    const int T,
    const PriorSpec& prior_spec,
    const ExpertSpec_FastSV& expert) {
  switch (prior_spec.sigma2.distribution) {
    case PriorSpec::Sigma2::GAMMA:
      if (expert.mh_blocking_steps == 1) {
        return 0.5 * (T - 1);
      } else {
        return 0.5 * T;
      }
    case PriorSpec::Sigma2::INVERSE_GAMMA:
      if (expert.mh_blocking_steps == 2) {
        return prior_spec.sigma2.inverse_gamma.shape + 0.5 * (T + 1);
      } else {
        return std::numeric_limits<double>::lowest();
      }
    default:
      return std::numeric_limits<double>::lowest();
  }
}

// Determine C_0 from the prior specification and the
// expert settings. This value is used in the update
// of sigma.
inline
double determine_C0(
    const PriorSpec& prior_spec) {
  switch (prior_spec.sigma2.distribution) {
    case PriorSpec::Sigma2::INVERSE_GAMMA:
      return prior_spec.sigma2.inverse_gamma.scale;
    default:
      return std::numeric_limits<double>::lowest();
  }
}

}  // END namespace fast_sv

namespace general_sv {

// Encapsulate the four parameters mu, phi, sigma,
// and rho as one returned object, including whether
// they were updated/accepted.
struct SampledTheta {
  double mu,
         phi,
         sigma,
         rho;
  bool mu_accepted,
       phi_accepted,
       sigma_accepted,
       rho_accepted;
};

namespace centered {  // sufficient augmentation

// Encapsulate the sufficient statistics for
// the sufficient augmentation 
struct SufficientStatistic {
  unsigned int length;  // T == length(y)
  double h_zero,  // h_0
         h_one,  // h_1
         sum_h,  // sum_{t=2}^{T-1} h_t
         sum_h_square,  // sum_{t=2}^{T-1} h_t^2
         h_last,  // h_T
         h_first_autocov,  // h_1 h_0
         sum_h_autocov,  // sum_{t=2}^T h_t h_{t-1}
         normalized_data,  // sum_{t=1}^{T-1} y_t / sqrt(tau_t) / exp(h_t / 2)
         normalized_data_square,  // sum_{t=1}^{T-1} y_t^2 / tau_t / exp(h_t)
         normalized_cov,  // sum_{t=1}^{T-1} y_t h_t / sqrt(tau_t) / exp(h_t / 2)
         normalized_leverage;  // sum_{t=1}^{T-1} y_t h_{t+1} / sqrt(tau_t) / exp(h_t / 2)
};

// Compute the log-likelihood for mu, phi, sigma,
// and rho in the exact SV model.
double theta_log_likelihood(
    const double mu,
    const double phi,
    const double sigma,
    const double rho,
    const SufficientStatistic& sufficient_statistic,
    const PriorSpec& prior_spec);

// Compute the sufficient statistic for mu, phi,
// sigma, and rho in the centered parameterization.
SufficientStatistic compute_sufficient_statistic(
    const arma::vec& data,
    const double h0,
    const arma::vec& h);

}  // END namespace centered

namespace noncentered {

// Encapsulate the sufficient statistics for
// the ancillary augmentation 
struct SufficientStatistic {
  const arma::vec& data;  // y_t / sqrt(tau_t)
  const arma::vec& c;  // ancillary latent states
};

// Compute the log-likelihood for mu, phi, sigma,
// and rho in the exact SV model.
double theta_log_likelihood(
    const double mu,
    const double phi,
    const double sigma,
    const double rho,
    const SufficientStatistic& sufficient_statistic,
    const PriorSpec& prior_spec);

}  // END namespace noncentered

// Compute the log-likelihood for mu, phi, sigma,
// and rho in the exact SV model.
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
    const Parameterization centering);

// Compute the log-prior density for mu, phi, sigma,
// and rho.
double theta_log_prior(
    const double mu,
    const double phi,
    const double sigma,
    const double rho,
    const PriorSpec& prior_spec);

// The transformation that is used to put mu, phi,
// sigma, and rho in an unbounded space for the
// random walk.
arma::vec4 theta_transform_inv(
    const double mu,
    const double phi,
    const double sigma,
    const double rho,
    const PriorSpec& prior_spec);

// Propose the next value for vector <mu, phi, sigma, rho>
// in the unbounded space. Also computes the log density
// of the old and new parameter values.
std::array<double, 6> theta_propose_rwmh(
    const double mu,
    const double phi,
    const double sigma,
    const double rho,
    const PriorSpec& prior_spec,
    const ProposalDiffusionKen& diffusion_ken,
    const arma::uvec& update_indicator);

}  // END namespace general_sv

}

#endif  // THETA_UTILS_H

