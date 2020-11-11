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
 * sampling_main.cc
 * 
 * Definitions of the functions declared in sampling_main.h.
 * Documentation can also be found in sampling_main.h.
 */

#include <RcppArmadillo.h>
#include <adaptation.hpp>
#include "sampling_main.h"
#include "single_update.h"
#include <type_definitions.h>
#include "utils_main.h"
#include "utils_latent_states.h"
#include "utils.h"

using namespace Rcpp;

namespace stochvol {

// Wrapper function around fast SV.
// See documentation above the declaration
List svsample_fast_cpp(
    const arma::vec& data,
    const int draws,
    const int burnin,
    const arma::mat& X,
    const Rcpp::List& priorspec_in,
    const int thinpara,
    const int thinlatent,
    const Rcpp::CharacterVector& keeptime_in,
    const Rcpp::List& startpara,
    const arma::vec& startlatent,
    const bool keep_tau,
    const Rcpp::List print_settings,
    const bool correct_model_specification,
    const bool interweave,
    const double offset,
    const Rcpp::List& expert_in) {
  const int T = data.size();
  const int p = X.n_cols;

  const PriorSpec prior_spec = list_to_priorspec(priorspec_in);
  const ExpertSpec_FastSV expert = list_to_fast_sv(expert_in, interweave);

  const bool is_regression = !::ISNA(X.at(0,0)),
             is_heavy_tail = prior_spec.nu.distribution != PriorSpec::Nu::INFINITE;

  if (prior_spec.mu.distribution == PriorSpec::Mu::CONSTANT && expert.mh_blocking_steps == 1) { // not implemented (would be easy, though)
    ::Rf_error("Single block update leaving mu constant is not yet implemented");
  }

  // shortcuts / precomputations
  const int thintime = determine_thintime(T, keeptime_in);

  // number of MCMC draws
  const int N = burnin + draws;

  // verbosity control
  const int chain = print_settings["chain"],
            n_chains = print_settings["n_chains"];
  const bool quiet = print_settings["quiet"],
             verbose = !quiet,
             single_chain = n_chains == 1;

  // initialize the variables:
  double mu = startpara["mu"],
         phi = startpara["phi"],
         sigma = startpara["sigma"],
         nu = startpara["nu"];
  arma::vec h = startlatent;  // contains h1 to hT, but not h0!
  double h0 = startpara["latent0"];
  arma::vec beta = startpara["beta"];
  arma::vec tau = expert_in["init_tau"];
  if (tau.n_elem == 1) {
    const double tau_elem = tau[0];
    tau.set_size(T);
    tau.fill(tau_elem);
  } else if (tau.n_elem != T) {
    ::Rf_error("Bad initialization for the vector tau. Should have length %d, received length %d, first element %f", T, tau.n_elem, tau[0]);
  }
  arma::uvec r = expert_in["init_indicators"];  // mixture indicators
  if (r.n_elem == 1) {
    const double r_elem = r[0];
    r.set_size(T);
    r.fill(r_elem);
  } else if (r.n_elem != T) {
    ::Rf_error("Bad initialization for the vector of mixture indicators. Should have length %d, received length %d, first element %f", T, r.n_elem, r[0]);
  }
  r -= 1u;
  R_assert(arma::all(r <= 9u), "Initial values of the mixture indicators need to be between 1 and 10 inclusive.");

  // keep some arrays cached
  arma::vec data_demean = is_regression ? data - X * beta : data,
            log_data2_normal = arma::log(arma::square(data_demean) / tau),  // standardized "data" (different for t-errors, "normal data")
            exp_h_half_inv;
  clamp_log_data2(log_data2_normal);
  const arma::vec conditional_mean(data.n_elem, arma::fill::zeros),
                  conditional_sd(data.n_elem, arma::fill::ones);

  // storage
  const bool keep_r = expert_in["store_indicators"];
  const int para_draws = draws / thinpara;
  Rcpp::NumericMatrix para_store(5, para_draws);
  Rcpp::NumericMatrix beta_store(p, is_regression * para_draws);
  const int latent_length = T / thintime;  // thintime must be either 1 or T
  const int latent_draws = draws / thinlatent;
  Rcpp::NumericVector latent0_store(latent_draws);
  Rcpp::NumericMatrix latent_store(latent_length, latent_draws);
  Rcpp::NumericMatrix tau_store(latent_length, keep_tau * latent_draws);
  Rcpp::IntegerMatrix r_store(latent_length, keep_r * latent_draws);
  Rcpp::NumericVector correction_weight_latent(correct_model_specification * latent_draws);
  Rcpp::NumericVector correction_weight_para(correct_model_specification * para_draws);

  // a warning about intended use
  if (correct_model_specification and (is_regression or is_heavy_tail)) {
    Rcpp::warning("Function 'svsample_fast_cpp' is not meant to do correction for model misspecification along with regression/heavy tails. Function 'svsample_general_cpp' can do this correction.");
  }

  // initializes the progress bar
  // "show" holds the number of iterations per progress sign
  int show = 0;
  if (verbose) {
    if (single_chain) {
      show = progressbar_init(N);
    } else {
      show = chain_print_init(chain, burnin, draws);
    }
  }

  for (int i = -burnin + 1; i < draws + 1; i++) {  // BEGIN main MCMC loop
    if (i % 20 == 0) {
      ::R_CheckUserInterrupt();
    }

    const bool parasave_round = i % thinpara == thinpara - 1,  // is this a parameter saving round?
               latentsave_round = i % thinlatent == thinlatent - 1;  // is this a latent saving round?

    // print a progress sign every "show" iterations
    if (verbose and i % show == 0) {
      if (single_chain) {
        progressbar_print();
      } else {
        show = chain_print(chain, i, burnin, draws);
      }
    }

    // a single MCMC update: update indicators, latent volatilities,
    // and parameters ONCE
    update_fast_sv(log_data2_normal, mu, phi, sigma, h0, h, r, prior_spec, expert);
    if (correct_model_specification or is_regression or is_heavy_tail) {
      exp_h_half_inv = arma::exp(-.5 * h);
    }

    // update tau and nu
    if (is_heavy_tail) {
      update_t_error(data_demean % exp_h_half_inv, tau, conditional_mean, conditional_sd, nu, prior_spec, false);
      // update cached data arrays after the regression
    }

    // update beta
    if (is_regression) {
      const arma::vec normalizer = exp_h_half_inv / arma::sqrt(tau);
      update_regressors(
          data % normalizer,
          X.each_col() % normalizer,
          beta, prior_spec);
      // update cached data arrays
      data_demean = data - X * beta;
      if (is_heavy_tail) {
        log_data2_normal = arma::log(arma::square(data_demean) / tau);
      } else {
        log_data2_normal = arma::log(arma::square(data_demean));
      }
      clamp_log_data2(log_data2_normal);
    } else if (is_heavy_tail) {
      log_data2_normal = arma::log(arma::square(data_demean) / tau);
      clamp_log_data2(log_data2_normal);
    } else {
      ;  // no update needed
    }


    // store draws
    double correction_weight_i = 1;
    if ((parasave_round or latentsave_round) and correct_model_specification) {
      correction_weight_i = fast_sv::compute_correction_weight(data_demean, log_data2_normal, h, 1/exp_h_half_inv);
    }
    if (i >= 1 and parasave_round) {
      const unsigned int index = i / thinpara - 1;
      save_para_sample(index, mu, phi, sigma, nu, beta, para_store, beta_store, is_regression);
      if (correct_model_specification) {
        correction_weight_para[index] = correction_weight_i;
      }
    }
    if (i >= 1 and latentsave_round) {
      const unsigned int index = i / thinlatent - 1;
      save_latent_sample(index, h0, h, tau, r, thintime, latent_length, latent0_store, latent_store, tau_store, r_store, keep_tau and is_heavy_tail, keep_r);
      if (correct_model_specification) {
        correction_weight_latent[index] = correction_weight_i;
      }
    }
  }  // END main MCMC loop

  if (verbose) {
    if (single_chain) {
      progressbar_finish(N);  // finalize progress bar
    } else {
      chain_print_finish(chain);
    }
  }

  Rcpp::NumericMatrix latent0_store_mat(latent0_store.length(), 1);
  latent0_store_mat(_, 0) = latent0_store;
  return cleanup(T, para_store, latent0_store_mat, latent_store, tau_store, beta_store, r_store, correction_weight_para, correction_weight_latent);
}

arma::ivec arma_sign(
    const arma::vec& vec) {
  arma::ivec d(vec.n_elem);
  std::transform(vec.cbegin(), vec.cend(), d.begin(), [](const double y_elem) -> int { return y_elem > 0 ? 1 : -1; });
  return d;
}

// Wrapper function around general SV.
// See documentation above the declaration
List svsample_general_cpp(
    const arma::vec& data,
    const int draws,
    const int burnin,
    const arma::mat& X,
    const Rcpp::List& priorspec_in,
    const int thinpara,
    const int thinlatent,
    const Rcpp::CharacterVector& keeptime_in,
    const Rcpp::List& startpara,
    const arma::vec& startlatent,
    const bool keep_tau,
    const Rcpp::List print_settings,
    const bool correct_model_specification,
    const bool interweave,
    const double offset,
    const Rcpp::List& expert_in) {

  const int N = burnin + draws;
  const int T = data.size();
  const int p = X.n_cols;

  const PriorSpec prior_spec = list_to_priorspec(priorspec_in);
  const ExpertSpec_GeneralSV expert = list_to_general_sv(expert_in, correct_model_specification, interweave);

  const bool is_regression = !::ISNA(X.at(0,0)),
             is_heavy_tail = prior_spec.nu.distribution != PriorSpec::Nu::INFINITE,
             is_leverage = not(prior_spec.rho.distribution == PriorSpec::Rho::CONSTANT and prior_spec.rho.constant.value == 0);

  const int thintime = determine_thintime(T, keeptime_in);

  // verbosity control
  const int chain = print_settings["chain"],
            n_chains = print_settings["n_chains"];
  const bool quiet = print_settings["quiet"],
             verbose = !quiet,
             single_chain = n_chains == 1;

  // starting values
  double mu = startpara["mu"];
  double phi = startpara["phi"];
  double sigma = startpara["sigma"];
  double nu = startpara["nu"];
  double rho = startpara["rho"];
  double h0 = startpara["latent0"];
  arma::vec h = startlatent;
  arma::vec beta = startpara["beta"];
  arma::vec tau = expert_in["init_tau"];
  if (tau.n_elem == 1) {
    const double tau_elem = tau[0];
    tau.set_size(T);
    tau.fill(tau_elem);
  } else if (tau.n_elem != T) {
    ::Rf_error("Bad initialization for the vector tau. Should have length %d, received length %d, first element %f", T, tau.n_elem, tau[0]);
  }

  // keep some arrays cached
  arma::vec data_demean = is_regression ? data - X * beta : data,
            data_demean_normal = data_demean / arma::sqrt(tau),  // standardized "data" (different for t-errors, "normal data")
            log_data2_normal = arma::log(arma::square(data_demean_normal)),
            exp_h_half_inv;
  clamp_log_data2(log_data2_normal);
  arma::ivec d = arma_sign(data_demean);
  auto conditional_moments = decorrelate(mu, phi, sigma, rho, h);

  // storage
  const int para_draws = draws / thinpara;
  Rcpp::NumericMatrix para_store(5, para_draws);
  Rcpp::NumericMatrix beta_store(p, is_regression * para_draws);
  const int latent_length = T / thintime;  // thintime must be either 1 or T
  const int latent_draws = draws / thinlatent;
  Rcpp::NumericVector latent0_store(latent_draws);
  Rcpp::NumericMatrix latent_store(latent_length, latent_draws);
  Rcpp::NumericMatrix tau_store(latent_length, keep_tau * latent_draws);

  // adaptive MH
  const int batch_size = 200,
            memory_size = expert.adapt ? expert.strategy.size() * (draws + burnin) / batch_size + 1 : 1;
  const double target_acceptance = 0.234,
               lambda = 0.1,
               init_scale = 0.001;
  AdaptationCollection adaptation_collection(
      (prior_spec.mu.distribution != PriorSpec::Mu::CONSTANT) +
      (prior_spec.phi.distribution != PriorSpec::Phi::CONSTANT) +
      (prior_spec.sigma2.distribution != PriorSpec::Sigma2::CONSTANT) +
      (prior_spec.rho.distribution != PriorSpec::Rho::CONSTANT),
      memory_size,
      batch_size,
      target_acceptance,
      lambda,  // between 0 and 1: the larger the value the stronger and longer the adaptation
      init_scale);

  // initializes the progress bar
  // "show" holds the number of iterations per progress sign
  int show = 0;
  if (verbose) {
    if (single_chain) {
      show = progressbar_init(N);
    } else {
      show = chain_print_init(chain, burnin, draws);
    }
  }

  for (int i = -burnin+1; i < draws+1; i++) {
    if (i % 20 == 0) {
      ::R_CheckUserInterrupt();
    }

    const bool parasave_round = i % thinpara == thinpara - 1,  // is this a parameter saving round?
               latentsave_round = i % thinlatent == thinlatent - 1;  // is this a latent saving round?

    // print a progress sign every "show" iterations
    if (verbose and i % show == 0) {
      if (single_chain) {
        progressbar_print();
      } else {
        show = chain_print(chain, i, burnin, draws);
      }
    }

    // update theta and h
    update_general_sv(data_demean_normal, log_data2_normal, d, mu, phi, sigma, rho, h0, h, adaptation_collection, prior_spec, expert);
    if (is_regression or is_heavy_tail) {
      exp_h_half_inv = arma::exp(-.5 * h);
      if (is_leverage) {
        conditional_moments = decorrelate(mu, phi, sigma, rho, h);
      }
    }

    // update tau and nu
    if (is_heavy_tail) {
      update_t_error(data_demean % exp_h_half_inv, tau, conditional_moments.conditional_mean, conditional_moments.conditional_sd, nu, prior_spec);
      // update cached data arrays after the regression
    }

    // update beta
    if (is_regression) {
      const arma::vec normalizer = exp_h_half_inv / (arma::sqrt(tau) % conditional_moments.conditional_sd);
      update_regressors(
          data % normalizer - conditional_moments.conditional_mean / conditional_moments.conditional_sd,
          X.each_col() % normalizer,
          beta, prior_spec);
      // update cached data arrays
      data_demean = data - X * beta;
      d = arma_sign(data_demean);
      if (is_heavy_tail) {
        log_data2_normal = arma::log(arma::square(data_demean) / tau);
      } else {
        log_data2_normal = arma::log(arma::square(data_demean));
      }
      clamp_log_data2(log_data2_normal);
    } else if (is_heavy_tail) {
      data_demean_normal = data_demean / arma::sqrt(tau);
      log_data2_normal = arma::log(arma::square(data_demean_normal));
      clamp_log_data2(log_data2_normal);
    } else {
      ;  // no update needed
    }

    // store draws
    if (i >= 1 and parasave_round) {
      save_para_sample(i / thinpara - 1, mu, phi, sigma, nu, rho, beta, para_store, beta_store, is_regression);
    }
    if (i >= 1 and latentsave_round) {
      save_latent_sample(i / thinlatent - 1, h0, h, tau, thintime, latent_length, latent0_store, latent_store, tau_store, keep_tau and is_heavy_tail);
    }
  }

  if (verbose) {
    if (single_chain) {
      progressbar_finish(N);  // finalize progress bar
    } else {
      chain_print_finish(chain);
    }
  }

  Rcpp::NumericMatrix latent0_store_mat(latent0_store.length(), 1);
  latent0_store_mat(_, 0) = latent0_store;
  return cleanup(T, para_store, latent0_store_mat, latent_store, tau_store, beta_store, adaptation_collection);
}

}

