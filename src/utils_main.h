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
#include <type_definitions.h>
#include <adaptation.hpp>

namespace stochvol {

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
    Rcpp::IntegerMatrix& mixture_indicators);

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

}

#endif  // _UTILS_MAIN_H_

