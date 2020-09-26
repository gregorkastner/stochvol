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
#include "type_definitions.h"
#include "utils_main.h"
#include "utils.h"

using namespace Rcpp;

namespace stochvol {

int determine_thintime(
    const int T,
    const Rcpp::CharacterVector& keeptime_in) {
  const std::string keeptime = as<std::string>(keeptime_in);
  if (keeptime == "all") {
    return 1;
  } else if (keeptime == "last") {
    return T;
  } else {
    Rf_error("Unknown value for 'keeptime'; got \"%s\"", keeptime.c_str());
  }
}

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
    const bool quiet,
    const bool correct_model_specification,
    const bool interweave,
    const double offset,
    const Rcpp::List& expert_in) {
  const int T = data.size();
  const int p = X.n_cols;

  // should we model the mean as well?
  const bool regression = !::ISNA(X.at(0,0));
  const PriorSpec prior_spec = list_to_priorspec(priorspec_in);
  const ExpertSpec_FastSV expert = list_to_fast_sv(expert_in, interweave);

  if (prior_spec.mu.distribution == PriorSpec::Mu::CONSTANT && expert.mh_blocking_steps == 1) { // not implemented (would be easy, though)
    ::Rf_error("Single block update leaving mu constant is not yet implemented");
  }

  // shortcuts / precomputations
  const int thintime = determine_thintime(T, keeptime_in);

  // number of MCMC draws
  const int N = burnin + draws;

  // verbosity control
  const bool verbose = !quiet;

  // t-errors:
  const bool terr = prior_spec.nu.distribution != PriorSpec::Nu::INFINITE;

  // initialize the variables:
  double mu = startpara["mu"],
         phi = startpara["phi"],
         sigma = startpara["sigma"],
         nu = startpara["nu"];
  arma::vec h = startlatent;  // contains h1 to hT, but not h0!
  double h0 = startpara["latent0"];
  arma::vec beta = startpara["beta"];

  const bool keep_r = expert_in["store_indicators"];
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

  // other forms of data
  arma::vec data_demean {data};
  arma::vec log_data2 = arma::log(arma::square(data) + offset);  // commonly used transformation
  arma::vec log_data2_normal = log_data2;  // standardized "data" (different for t-errors, "normal data")

  // some stuff for the t-errors
  arma::vec tau(T, arma::fill::ones);

  // storage
  const int para_draws = draws / thinpara + 1;
  Rcpp::NumericMatrix para_store(5, para_draws);
  Rcpp::NumericMatrix beta_store(p, regression * para_draws);
  const int latent_length = T / thintime;  // thintime must be either 1 or T
  const int latent_draws = draws / thinlatent + 1;
  Rcpp::NumericVector latent0_store(latent_draws);
  Rcpp::NumericMatrix latent_store(latent_length, latent_draws);
  Rcpp::NumericMatrix tau_store(latent_length, keep_tau * latent_draws);
  Rcpp::IntegerMatrix r_store(latent_length, keep_r * latent_draws);

  // initializes the progress bar
  // "show" holds the number of iterations per progress sign
  int show = 0;
  if (verbose) show = progressbar_init(N);

  for (int i = -burnin + 1; i < draws + 1; i++) {  // BEGIN main MCMC loop
    if (i % 20 == 0) {
      ::R_CheckUserInterrupt();
    }

    const bool thinpara_round = (thinpara > 1) && (i % thinpara != 0);  // is this a parameter thinning round?
    const bool thinlatent_round = (thinlatent > 1) && (i % thinlatent != 0);  // is this a latent thinning round?

    // print a progress sign every "show" iterations
    if (verbose && (i % show == 0)) progressbar_print();

    if (regression) {
      data_demean = data - X*beta;
      log_data2_normal = log_data2 = arma::log(arma::square(data_demean));
    }

    if (terr) {
      update_t_error(data_demean % arma::exp(-0.5 * h), tau, nu, prior_spec);
      log_data2_normal = log_data2 - arma::log(tau);
    }

    // a single MCMC update: update indicators, latent volatilities,
    // and parameters ONCE
    update_fast_sv(log_data2_normal, mu, phi, sigma, h0, h, r, prior_spec, expert);

    if (regression) { // update betas (regression)
      update_regressors(data, X, beta, tau, arma::exp(-0.5 * h), nu, prior_spec);
    }

    // store draws
    if ((i >= 1) && !thinpara_round) {
      save_para_sample(i / thinpara - 1, mu, phi, sigma, nu, beta, para_store, beta_store, regression);
    }
    if ((i >= 1) && !thinlatent_round) {
      save_latent_sample(i / thinlatent - 1, h0, h, tau, r, thintime, latent_length, latent0_store, latent_store, tau_store, r_store, keep_tau and terr, keep_r);
    }
  }  // END main MCMC loop

  if (verbose) progressbar_finish(N);  // finalize progress bar

  Rcpp::NumericMatrix latent0_store_mat(latent0_store.length(), 1);
  latent0_store_mat(_, 0) = latent0_store;
  return cleanup(T, para_store, latent0_store_mat, latent_store, tau_store, beta_store, r_store);
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
    const bool quiet,
    const bool correct_model_specification,
    const bool interweave,
    const double offset,
    const Rcpp::List& expert_in) {

  const int N = burnin + draws;
  const bool regression = !::ISNA(X.at(0,0));
  const int T = data.size();
  const int p = X.n_cols;

  const PriorSpec prior_spec = list_to_priorspec(priorspec_in);
  const ExpertSpec_GeneralSV expert = list_to_general_sv(expert_in, correct_model_specification, interweave);

  const bool verbose = !quiet;
  const int thintime = determine_thintime(T, keeptime_in);

  arma::vec data_demean {data};
  arma::vec log_data2 = arma::log(arma::square(data) + offset);  // commonly used transformation
  arma::vec log_data2_normal = log_data2;  // standardized "data" (different for t-errors, "normal data")
  arma::ivec d = arma_sign(data);

  // t-errors:
  const bool terr = prior_spec.nu.distribution != PriorSpec::Nu::INFINITE;
  arma::vec tau(T, arma::fill::ones);

  // starting values
  double mu = startpara["mu"];
  double phi = startpara["phi"];
  double sigma = startpara["sigma"];
  double nu = startpara["nu"];
  double rho = startpara["rho"];
  double h0 = startpara["latent0"];
  arma::vec h = startlatent, ht = centered_to_noncentered(mu, sigma, h);
  arma::vec beta = startpara["beta"];

  // storage
  const int para_draws = draws / thinpara + 1;
  Rcpp::NumericMatrix para_store(5, para_draws);
  Rcpp::NumericMatrix beta_store(p, regression * para_draws);
  const int latent_length = T / thintime;  // thintime must be either 1 or T
  const int latent_draws = draws / thinlatent + 1;
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
  const int show = verbose ? progressbar_init(N) : 0;

  for (int i = -burnin+1; i < draws+1; i++) {
    if (i % 20 == 0) {
      ::R_CheckUserInterrupt();
    }

    const bool thinpara_round = (thinpara > 1) and (i % thinpara != 0);  // is this a parameter thinning round?
    const bool thinlatent_round = (thinlatent > 1) and (i % thinlatent != 0);  // is this a latent thinning round?

    // print a progress sign every "show" iterations
    if (verbose && (i % show == 0)) progressbar_print();

    if (regression) {
      data_demean = data - X*beta;
      log_data2_normal = log_data2 = arma::log(arma::square(data_demean));
      d = arma_sign(data_demean);
    }

    if (terr) {
      update_t_error(data_demean % arma::exp(-0.5 * h), tau, nu, prior_spec);  // TODO correct for rho
      log_data2_normal = log_data2 - arma::log(tau);
    }

    // update theta and h
    update_general_sv(data_demean, log_data2_normal, d, mu, phi, sigma, rho, h0, h, adaptation_collection, prior_spec, expert);
    ht = centered_to_noncentered(mu, sigma, h);

    // update beta
    if (regression) {
      arma::vec normalizer = arma::exp(-h/2);
      normalizer.head(T-1) /= std::sqrt(1 - std::pow(rho, 2));
      arma::vec data_reg = data;
      data_reg.head(T-1) -= rho * (arma::exp(h.head(T-1)/2) % (ht.tail(T-1) - phi*ht.head(T-1)));

      update_regressors(data_reg, X, beta, tau, normalizer, nu, prior_spec);
    }

    // store draws
    if ((i >= 1) && !thinpara_round) {
      save_para_sample(i / thinpara - 1, mu, phi, sigma, nu, rho, beta, para_store, beta_store, regression);
    }
    if ((i >= 1) && !thinlatent_round) {
      save_latent_sample(i / thinlatent - 1, h0, h, tau, thintime, latent_length, latent0_store, latent_store, tau_store, keep_tau and terr);
    }
  }

  if (verbose) progressbar_finish(N);  // finalize progress bar

  Rcpp::NumericMatrix latent0_store_mat(latent0_store.length(), 1);
  latent0_store_mat(_, 0) = latent0_store;
  return cleanup(T, para_store, latent0_store_mat, latent_store, tau_store, beta_store, adaptation_collection);
}

}

