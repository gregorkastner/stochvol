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

// Single update of the marginal data augmentation for the degrees of freedom.
void update_tau(
    const arma::vec& homosked_data,
    arma::vec& tau,
    const double nu) {
  // Watch out, R::rgamma(shape, scale), not Rf_rgamma(shape, rate)
  std::transform(
      homosked_data.cbegin(), homosked_data.cend(),
      tau.begin(),
      [nu](const double homosked_data_i) -> double { return 1./R::rgamma((nu + 1.) / 2., 2. / (nu + std::pow(homosked_data_i, 2))); });
}

// Single update of the degrees of freedom.
// See documentation above the declaration
void update_t_error(
    const arma::vec& homosked_data,
    arma::vec& tau,
    double& nu,
    const PriorSpec& prior_spec) {
  R_assert(prior_spec.nu.distribution != PriorSpec::Nu::INFINITE, "Call to update_t_error: Non-matching model specification. Prior for nu should not be infinity.");

  update_tau(homosked_data, tau, nu);
  // TODO since C++17 this can be done with transform_reduce (in parallel)
  const double sumtau = std::accumulate(
      tau.cbegin(), tau.cend(),
      0.,
      [](const double partial_sum, const double tau_i) -> double { return partial_sum + std::log(tau_i) + 1. / tau_i;});

  const int T = homosked_data.size();

  if (prior_spec.nu.distribution == PriorSpec::Nu::EXPONENTIAL) {
    const double lambda = prior_spec.nu.exponential.rate;
    const double numean = newton_raphson(nu, sumtau, T, lambda);
    const double auxsd = std::sqrt(-1/ddlogdnu(numean, T)); 
    const double nuprop = R::rnorm(numean, auxsd);
    const double logR =
      logdnu(nuprop, sumtau, lambda, T) - logdnorm(nuprop, numean, auxsd) -
      (logdnu(nu, sumtau, lambda, T) - logdnorm(nu, numean, auxsd));

    if (logR >= 0 or std::log(R::unif_rand()) < logR) {
      nu = nuprop;
    }
  }
}

// Single update of the coefficients in a Bayesian linear regression.
// See documentation above the declaration
void update_regressors(
    arma::vec dependent_variable,  // by value on purpose
    arma::mat independent_variables,  // by value on purpose
    arma::vec& beta,
    arma::vec& tau,
    const arma::vec& normalizer,  // inverse of heteroskedastic scales
    const double df,
    const PriorSpec& prior_spec) {
  dependent_variable %= normalizer;
  independent_variables.each_col() %= normalizer;
  if (df < 10000) {  // heavy-tailed
    update_tau(dependent_variable - independent_variables * beta, tau, df);
    dependent_variable /= tau;
    independent_variables.each_col() /= tau;
  }

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
    //const Proposal proposal = use_mala ? Proposal::MALA : Proposal::RWMH;  TODO delete
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
            general_sv::theta_transform_inv(mu, phi, sigma, rho).elem(arma::find(update_indicator)));  // current sample
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

