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
 * single_update.cc
 * 
 * Definitions of the functions declared in single_update.h.
 * Documentation can also be found in single_update.h.
 */

#include <RcppArmadillo.h>
#include "single_update.h"
#include "densities.h"
#include "type_definitions.h"
#include "sampling_latent_states.h"
#include "sampling_parameters.h"
#include "utils_latent_states.h"
#include "utils_parameters.h"
#include "utils.h"

using namespace Rcpp;

namespace stochvol {

namespace fast_sv {

// Determine ASIS from 'expert'
std::vector<Parameterization> expert_to_strategy(
    const ExpertSpec_FastSV& expert) {
  switch (expert.baseline) {
    case Parameterization::CENTERED:
      if (expert.interweave) {
        return {Parameterization::CENTERED, Parameterization::NONCENTERED};
      } else {
        return {Parameterization::CENTERED};
      }
    case Parameterization::NONCENTERED:
      if (expert.interweave) {
        return {Parameterization::NONCENTERED, Parameterization::CENTERED};
      } else {
        return {Parameterization::NONCENTERED};
      }
  }
}

}

// Single update accoring to fast SV.
// See documentation above the declaration
void update_fast_sv(
    const arma::vec& log_data2,
    double& mu,
    double& phi,
    double& sigma,
    double& h0,
    arma::vec& h,
    arma::uvec& r,
    const PriorSpec& prior_spec,  // C0, cT, Bsigma, a0, b0, bmu, Bmu, Gammaprior, dontupdatemu, priorlatent0 feed into this
    const ExpertSpec_FastSV& expert) {  // parameterization, centered_baseline, B011inv, B022inv, truncnormal, Mhcontrol, MHsteps feed into this
  // TODO setup validation (put in a "debug" environment in the end)
  // phi is beta
  // sigma2 is either inverse_gamma with shape == 2.5 or gamma with shape == 0.5
  // inverse_gamma prior and mh_blocking_steps != 2 is nyi
  // constant mu implies mh_blocking_steps == 3 (mh_blocking_steps == 1 is nyi and 2 does not make sense because mu is drawn jointly with either phi or sigma2 in the 2-block)

  double ht0 = centered_to_noncentered(mu, sigma, h0);
  arma::vec ht = centered_to_noncentered(mu, sigma, h);

  if (prior_spec.mu.distribution == PriorSpec::Mu::CONSTANT) {
    mu = 0;
  }

  if (expert.update.mixture_indicators) {  // Step (c): sample indicators
    r = fast_sv::draw_mixture_indicators(log_data2, mu, phi, sigma, h);
  }

  if (expert.update.latent_vector) {  // Step (a): sample the latent volatilities h:
    const auto latent_new = fast_sv::draw_latent(log_data2, mu, phi, sigma, r, prior_spec, expert);

    switch (expert.baseline) {
      case Parameterization::CENTERED:
        h = latent_new.h;
        h0 = latent_new.h0;
        ht = centered_to_noncentered(mu, sigma, h);
        ht0 = centered_to_noncentered(mu, sigma, h0);
        break;
      case Parameterization::NONCENTERED:
        ht = latent_new.h;
        ht0 = latent_new.h0;
        h = noncentered_to_centered(mu, sigma, ht);
        h0 = noncentered_to_centered(mu, sigma, ht0);
        break;
    }
  }

  if (expert.update.parameters) {  // Step (b): sample mu, phi, sigma
    const auto strategy = fast_sv::expert_to_strategy(expert);
    for (const auto parameterization : strategy) {
      // Draw theta
      const auto parameter_draw = fast_sv::draw_theta(log_data2, mu, phi, sigma, h0, ht0, h, ht, r, prior_spec, expert, parameterization);
      mu = parameter_draw.mu;
      phi = parameter_draw.phi;
      sigma = parameter_draw.sigma;
      // Update latent vectors
      switch (parameterization) {
        case Parameterization::CENTERED:
          ht = centered_to_noncentered(mu, sigma, h);
          ht0 = centered_to_noncentered(mu, sigma, h0);
          break;
        case Parameterization::NONCENTERED:
          h = noncentered_to_centered(mu, sigma, ht);
          h0 = noncentered_to_centered(mu, sigma, ht0);
          break;
      }
    }
  }
}

// Single update of a single tau: the marginal data augmentation for the degrees of freedom.
double update_single_tau(
    const double homosked_data_i,
    const double tau_i,
    const double mean_i,
    const double sd_i,
    const double nu,
    const bool do_tau_acceptance_rejection) {
  // Watch out, R::rgamma(shape, scale), not Rf_rgamma(shape, rate)
  const double proposal_shape = (nu + 1) / 2.,
               proposal_scale = (nu - 2 + std::pow(homosked_data_i, 2)) / 2.,
               tau_proposal_i = 1./R::rgamma(proposal_shape, 1 / proposal_scale);  // non-parallelizable
  if (do_tau_acceptance_rejection) {
    const double sqrt_tau_prop = std::sqrt(tau_proposal_i),
                 sqrt_tau = std::sqrt(tau_i),
                 log_ar =
                   (logdnorm(homosked_data_i, sqrt_tau_prop * mean_i, sqrt_tau_prop * sd_i) +
                    logdinvgamma(tau_proposal_i, .5 * nu, .5 * (nu - 2)) -
                    logdinvgamma(tau_proposal_i, proposal_shape, proposal_scale)) -
                   (logdnorm(homosked_data_i, sqrt_tau * mean_i, sqrt_tau * sd_i) +
                    logdinvgamma(tau_i, .5 * nu, .5 * (nu - 2)) -
                    logdinvgamma(tau_i, proposal_shape, proposal_scale));
    if (log_ar >= 0 or R::unif_rand() < std::exp(log_ar)) {
      return tau_proposal_i;
    } else {
      return tau_i;
    }
  } else {
    return tau_proposal_i;
  }
}

// Single update of the degrees of freedom.
// See documentation above the declaration
void update_t_error(
    const arma::vec& homosked_data,
    arma::vec& tau,
    const arma::vec& mean,
    const arma::vec& sd,
    double& nu,
    const PriorSpec& prior_spec,
    const bool do_tau_acceptance_rejection) {
  R_assert(prior_spec.nu.distribution != PriorSpec::Nu::INFINITE, "Call to update_t_error: Non-matching model specification. Prior for nu should not be infinity.");

  const int T = homosked_data.n_elem;

  double sum_tau = 0;
  for (int i = 0; i < T; i++) {
    tau[i] = update_single_tau(homosked_data[i], tau[i], mean[i], sd[i], nu, do_tau_acceptance_rejection);
    sum_tau += std::log(tau[i]) + 1. / tau[i];
  }

  if (prior_spec.nu.distribution == PriorSpec::Nu::EXPONENTIAL) {
    const double lambda = prior_spec.nu.exponential.rate;
    const double numean = newton_raphson(nu, sum_tau, T, lambda);
    const double auxsd = std::sqrt(-1/ddlogdnu(numean, T)); 
    const double nuprop = R::rnorm(numean, auxsd);
    const double log_ar =
      logdnu(nuprop, sum_tau, lambda, T) - logdnorm(nuprop, numean, auxsd) -
      (logdnu(nu, sum_tau, lambda, T) - logdnorm(nu, numean, auxsd));

    if (log_ar >= 0 or R::unif_rand() < std::exp(log_ar)) {
      nu = nuprop;
    }
  }
}

// Single update of the coefficients in a Bayesian linear regression.
// See documentation above the declaration
void update_regressors(
    const arma::vec& dependent_variable,
    const arma::mat& independent_variables,
    arma::vec& beta,
    const PriorSpec& prior_spec) {
  const int p = independent_variables.n_cols;
  arma::mat postprecchol;
  arma::mat postpreccholinv;
  arma::mat postcov;
  arma::vec postmean;
  arma::vec armadraw(p);

  bool success = true;
  // cholesky factor of posterior precision matrix
  success = success and arma::chol(postprecchol, independent_variables.t() * independent_variables + prior_spec.beta.multivariate_normal.precision);
  // inverse cholesky factor of posterior precision matrix 
  success = success and arma::inv(postpreccholinv, arma::trimatu(postprecchol));
  if (!success) {
    Rcpp::stop("Cholesky or its inverse failed");
  }

  // posterior covariance matrix and posterior mean vector
  postcov = postpreccholinv * postpreccholinv.t();
  postmean = postcov * (independent_variables.t() * dependent_variable + prior_spec.beta.multivariate_normal.precision * prior_spec.beta.multivariate_normal.mean);

  armadraw.imbue([]() -> double { return R::norm_rand(); });

  // posterior betas
  beta = postmean + postpreccholinv * armadraw;
}

// Single update accoring to general SV.
// See documentation above the declaration
void update_general_sv(
    const arma::vec& data,
    const arma::vec& log_data2,
    const arma::ivec& sign_data,
    double& mu,
    double& phi,
    double& sigma,
    double& rho,
    double& h0,
    arma::vec& h,
    AdaptationCollection& adaptation_collection,
    const PriorSpec& prior_spec,  // prior_mu, prior_phi, prior_sigma2, prior_rho, gammaprior, dontupdatemu feed into this (plus priorlatent0, truncnormal nyi)
    const ExpertSpec_GeneralSV& expert) {  // strategy, correct, use_mala feed into this
  double ht0;
  arma::vec ht;
  ht0 = centered_to_noncentered(mu, sigma, h0);
  ht = centered_to_noncentered(mu, sigma, h);

  if (expert.update.latent_vector) {  // Sample the latent states
    //const  // not const to be able to std::move
    LatentVector h_full = general_sv::draw_latent(data, log_data2, sign_data, mu, phi, sigma, rho, h, ht, prior_spec, expert);
    h0 = h_full.h0;
    h = std::move(h_full.h);
    ht0 = centered_to_noncentered(mu, sigma, h0);
    ht = centered_to_noncentered(mu, sigma, h);
  }

  if (expert.update.parameters) {  // Sample the parameters
    const arma::uvec update_indicator {
          prior_spec.mu.distribution != PriorSpec::Mu::CONSTANT,
          prior_spec.phi.distribution != PriorSpec::Phi::CONSTANT,
          prior_spec.sigma2.distribution != PriorSpec::Sigma2::CONSTANT,
          prior_spec.rho.distribution != PriorSpec::Rho::CONSTANT};

    arma::vec exp_h_half {arma::exp(.5 * h)};  // cache exp() calculations
    arma::vec exp_h_half_proposal_nc;
    for (const Parameterization par : expert.strategy) {
      Adaptation& adaptation = adaptation_collection[par];
      const ProposalDiffusionKen& adapted_proposal = expert.adapt ? adaptation.get_proposal() : expert.proposal_diffusion_ken;
      const auto theta = general_sv::draw_theta(
          data,
          mu, phi, sigma, rho,
          h0, ht0, h, ht,
          exp_h_half, exp_h_half_proposal_nc,
          update_indicator,
          prior_spec, expert,
          adapted_proposal,
          par);
      mu = theta.mu, phi = theta.phi, sigma = theta.sigma, rho = theta.rho;
      const bool accepted = theta.mu_accepted or theta.phi_accepted or theta.sigma_accepted or theta.rho_accepted;  // was anything accepted?
      if (expert.adapt) {
        adaptation.register_sample(
            accepted,
            general_sv::theta_transform_inv(mu, phi, sigma, rho, prior_spec).elem(arma::find(update_indicator)));  // current sample
      }

      if (accepted) {
        switch (par) {
          case Parameterization::CENTERED:
            ht0 = centered_to_noncentered(mu, sigma, h0);
            ht = centered_to_noncentered(mu, sigma, h);
            break;
          case Parameterization::NONCENTERED:
            h = noncentered_to_centered(mu, sigma, ht);
            h0 = noncentered_to_centered(mu, sigma, ht0);
            exp_h_half = std::move(exp_h_half_proposal_nc);
            break;
        }
      }
    }
  }

  return;
}

}

