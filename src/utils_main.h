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
 * utils_main.h
 * 
 * Utility functions related to the sampling_main.cc.
 */

#ifndef _UTILS_MAIN_H_
#define _UTILS_MAIN_H_

// Helper functions for the main sampler

#include <RcppArmadillo.h>
#include <expert.hpp>
#include <adaptation.hpp>

namespace stochvol {

// De-correlate the observations
//
// This function is relevant when rho is non-zero. Any regression
// conditional on the asymmetric SV model, e.g. when asymmetric SV
// is embedded in a Bayesian linear regression, needs to correct
// both the observations and the error weights. That is so because
// conditional on rho and h, the observation errors can be decomposed
// and thus the conditional errors have non-zero mean and modified
// variance (1-rho^2).
//
// The function 'decorrelate' stands for "removing" rho from
// the model. It returns the conditional mean and the conditional
// standard deviation.
struct AsymmetricConditionalMoments {
  arma::vec conditional_mean;
  arma::vec conditional_sd;
};
inline
AsymmetricConditionalMoments decorrelate (
    const double mu,
    const double phi,
    const double sigma,
    const double rho,
    const arma::vec& h) {
  arma::vec mean(h.n_elem, arma::fill::zeros);
  arma::vec sd(h.n_elem, arma::fill::ones);
  sd.head(sd.n_elem - 1) *= std::sqrt(1 - std::pow(rho, 2));
  
  const double rho_div_sigma = rho / sigma;
  auto it_mean = mean.begin();
  auto it_h = h.cbegin();
  for (unsigned int t = 0; t < h.n_elem-1; t++) {
    *it_mean = rho_div_sigma * ((*(it_h+1)) - mu - phi * ((*it_h) - mu));
    // Step
    it_mean++;
    it_h++;
  }
  return {std::move(mean), std::move(sd)};
}

// Avoid values of -inf
inline
void clamp_log_data2(
    arma::vec& log_data2) {
  // -100 ~= log(4e-44)
  std::for_each(log_data2.begin(), log_data2.end(), [](double& value) { value = std::max(value, -100.0); });
  if (not arma::is_finite(log_data2)) {
    ::Rf_error("Non-finite (+-inf or NaN) elements in the data set. This should not happen. It would help us if you could contact the maintainer with a reproducible example.");
  }
}

namespace fast_sv {

// Compute log-weight for fast SV re-sampling.
// Relevant when correcting for model mis-specification.
double compute_correction_weight(
    const arma::vec& data,
    const arma::vec& log_data2,
    const arma::vec& h,
    const arma::vec& exp_h_half);

}  // END namespace fast_sv

// Save a single sample in the storage objects
inline
void save_para_sample(
    const int index_para,
    const double mu,
    const double phi,
    const double sigma,
    const double nu,
    const double rho,
    const arma::vec& beta,
    Rcpp::NumericMatrix& para_store,
    Rcpp::NumericMatrix& beta_store,
    const bool save_beta) {
  para_store(0, index_para) = mu;
  para_store(1, index_para) = phi;
  para_store(2, index_para) = sigma;
  para_store(3, index_para) = nu;
  para_store(4, index_para) = rho;
  if (save_beta) {
    std::copy(beta.cbegin(), beta.cend(), beta_store(Rcpp::_, index_para).begin());
  }
}

// Save a single sample in the storage objects
inline
void save_para_sample(
    const int index_para,
    const double mu,
    const double phi,
    const double sigma,
    const double nu,
    const arma::vec& beta,
    Rcpp::NumericMatrix& para_store,
    Rcpp::NumericMatrix& beta_store,
    const bool save_beta) {
  para_store(0, index_para) = mu;
  para_store(1, index_para) = phi;
  para_store(2, index_para) = sigma;
  para_store(3, index_para) = nu;
  para_store(4, index_para) = 0;
  if (save_beta) {
    std::copy(beta.cbegin(), beta.cend(), beta_store(Rcpp::_, index_para).begin());
  }
}

// Save a single sample in the storage objects
inline
void save_latent_sample(
    const int index_latent,
    const double h0,
    const arma::vec& h,
    const arma::vec& tau,
    const int thintime,
    const int latent_length,
    Rcpp::NumericVector& latent0_store,
    Rcpp::NumericMatrix& latent_store,
    Rcpp::NumericMatrix& tau_store,
    const bool save_tau) {
  latent0_store[index_latent] = h0;
  for (int volind = 0, thincol = 0; thincol < latent_length; volind++, thincol++) {
    latent_store(volind, index_latent) = h[thintime * (thincol + 1) - 1];
  }
  if (save_tau) {
    for (int volind = 0, thincol = 0; thincol < latent_length; volind++, thincol++) {
      tau_store.at(volind, index_latent) = tau[thintime * (thincol + 1) - 1];
    }
  }
}

// Save a single sample in the storage objects
inline
void save_latent_sample(
    const int index_latent,
    const double h0,
    const arma::vec& h,
    const arma::vec& tau,
    const arma::uvec& r,
    const int thintime,
    const int latent_length,
    Rcpp::NumericVector& latent0_store,
    Rcpp::NumericMatrix& latent_store,
    Rcpp::NumericMatrix& tau_store,
    Rcpp::IntegerMatrix& r_store,
    const bool save_tau,
    const bool save_r) {
  save_latent_sample(index_latent, h0, h, tau, thintime, latent_length, latent0_store, latent_store, tau_store, save_tau);
  if (save_r) {
    for (int volind = 0, thincol = 0; thincol < latent_length; volind++, thincol++) {
      r_store.at(volind, 0) = r[thintime * (thincol + 1) - 1];
    }
  }
}

// Sums up results and prepares return value
Rcpp::List cleanup(
    const int T,
    Rcpp::NumericMatrix& para,
    Rcpp::NumericMatrix& latent0,
    Rcpp::NumericMatrix& latent,
    Rcpp::NumericMatrix& tau,
    Rcpp::NumericMatrix& betas,
    Rcpp::IntegerMatrix& mixture_indicators,
    Rcpp::NumericVector& correction_weight_para,
    Rcpp::NumericVector& correction_weight_latent);

// Sums up results and prepares return value
Rcpp::List cleanup(
    const int T,
    Rcpp::NumericMatrix& para,
    Rcpp::NumericMatrix& latent0,
    Rcpp::NumericMatrix& latent,
    Rcpp::NumericMatrix& tau,
    Rcpp::NumericMatrix& betas,
    AdaptationCollection& adaptation_collection);

// Sets up the progress bar
int progressbar_init(
    const int N);

// Adds one '+' to progress bar
inline
void progressbar_print() {
  ::REprintf("+");
  ::R_FlushConsole();
}

// Finalizes progress bar
void progressbar_finish(
    const int N);

// Sets up parallel print
inline
int chain_print_init(
    const int chain,
    const int burnin,
    const int draws) {
  ::REprintf("Chain %d starting\n", chain);
  ::R_FlushConsole();
  const int next_big_checkpoint = burnin <= 0 ? draws : burnin;
  if (next_big_checkpoint < 50) {
    return next_big_checkpoint;
  } else if (next_big_checkpoint < 200) {
    return next_big_checkpoint / 2;
  } else if (next_big_checkpoint < 500) {
    return next_big_checkpoint / 5;
  } else {
    return next_big_checkpoint / 10;
  }
}

// Prints progress for parallel chains
inline
int chain_print(
    const int chain,
    const int i,
    const int burnin,
    const int draws) {
  if (i < 0) {
    ::REprintf("Chain %d at iteration %d / %d (warmup)\n", chain, i+burnin, burnin+draws);
  } else {
    ::REprintf("Chain %d at iteration %d / %d (sampling)\n", chain, i+burnin, burnin+draws);
  }
  ::R_FlushConsole();
  const int next_big_checkpoint = i < 0 ? burnin : draws,
            remaining = std::abs(next_big_checkpoint - i);
  int step;
  if (next_big_checkpoint < 50) {
    step = next_big_checkpoint;
  } else if (next_big_checkpoint < 200) {
    step = next_big_checkpoint / 2;
  } else if (next_big_checkpoint < 500) {
    step = next_big_checkpoint / 5;
  } else {
    step = next_big_checkpoint / 10;
  }
  if (remaining <= step) {
    return remaining;
  } else {
    return step;
  }
}

// Finalizes parallel print
inline
void chain_print_finish(
    const int chain) {
  ::REprintf("Chain %d done\n", chain);
  ::R_FlushConsole();
}

// Transform an R priorspec list to its corresponding
// C++ object
PriorSpec list_to_priorspec(
    const Rcpp::List& list);

// Transform an R general expert list to its corresponding
// C++ object
ExpertSpec_GeneralSV list_to_general_sv(
    const Rcpp::List& list,
    const bool correct_model_specification,
    const bool interweave);

// Transform an R fast expert list to its corresponding
// C++ object
ExpertSpec_FastSV list_to_fast_sv(
    const Rcpp::List& list,
    const bool interweave);

// Transform an R adaptation list into its corresponding
// C++ object
Adaptation list_to_adaptation (
    const Rcpp::List& list);

// Transform an R adaptation collection list into
// its corresponding C++ object
AdaptationCollection list_to_adaptationcollection (
    const Rcpp::List& list);

}

#endif  // _UTILS_MAIN_H_

