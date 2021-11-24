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
 * utils_latent_states.h
 *
 * Constants and utility functions related to the sampling of the
 * latent states or the mixture indicators.
 *
 * The functions are separated into the fast_sv and the general_sv
 * variants due to the varying expert options.
 */

#ifndef _UTILS_LATENT_STATES_H_
#define _UTILS_LATENT_STATES_H_

#include <RcppArmadillo.h>

namespace stochvol {

// Constants and their transforms from Omori et al (2007).
// These define the auxiliary mixture approximation used both
// in the fast SV and in the general SV samplers.
const arma::vec::fixed<10> mix_prob {.00609, .04775, .13057, .20674, .22715, .18842, .12047, .05591, .01575, .00115};
const arma::vec::fixed<10> mix_mean {1.92677, 1.34744, .73504, .02266, -.85173, -1.97278, -3.46788, -5.55246, -8.68384, -14.65000};
const arma::vec::fixed<10> mix_var {.11265, .17788, .26768, .40611, .62699, .98583, 1.57469, 2.54498, 4.16591, 7.33342};
const arma::vec::fixed<10> mix_a {1.01418, 1.02248, 1.03403, 1.05207, 1.08153, 1.13114, 1.21754, 1.37454, 1.68327, 2.50097};
const arma::vec::fixed<10> mix_b {0.50710, 0.51124, 0.51701, 0.52604, 0.54076, 0.56557, 0.60877, 0.68728, 0.84163, 1.25049};
const arma::vec::fixed<10> mix_sd {arma::sqrt(mix_var)};
const arma::vec::fixed<10> mix_varinv {1 / mix_var};
const arma::vec::fixed<10> mix_2varinv {0.5 * mix_varinv};
const arma::vec::fixed<10> mix_pre {  // TODO what are these numbers?
-4.0093723912083900628999799664597958326339721679687500000,
-2.1784531553855770447114537091692909598350524902343750000,
-1.3768642766903782526100030736415646970272064208984375000,
-1.1257277037836319610875079888501204550266265869140625000,
-1.2487323430568648685579091761610470712184906005859375000,
-1.6619460888428292388852014482836239039897918701171875000,
-2.3433837334574310062862423365004360675811767578125000000,
-3.3510734196563021214387845247983932495117187500000000000,
-4.8643822832849297199686589010525494813919067382812500000,
-7.7642143280080739842219372803810983896255493164062500000};

// Export the Omori constants
inline
Rcpp::List get_omori_constants () {
  return Rcpp::List::create(
      Rcpp::_["prob"] = mix_prob,
      Rcpp::_["mean"] = mix_mean,
      Rcpp::_["var"] = mix_var,
      Rcpp::_["a"] = mix_a,
      Rcpp::_["b"] = mix_b);
}

namespace fast_sv {

// Encapsulation of chol_diag and chol_offdiag
struct CholeskyTridiagonal {
  arma::vec chol_diag,
            chol_offdiag;
};

// Cholesky factor for a tridiagonal matrix with constant off-diagonal
CholeskyTridiagonal cholesky_tridiagonal(
    const arma::vec& omega_diag,
    const double omega_offdiag);

// Solves Chol*x = covector ("forward algorithm")
arma::vec forward_algorithm(
    const arma::vec& chol_diag,
    const arma::vec& chol_offdiag,
    const arma::vec& covector);

// Solves (Chol')*x = htmp ("backward algorithm")
arma::vec backward_algorithm(
    const arma::vec& chol_diag,
    const arma::vec& chol_offdiag,
    const arma::vec& htmp);

// Draws length(r) RVs, expects the non-normalized CDF mixprob
arma::uvec inverse_transform_sampling(
    const arma::vec& mixprob,
    const int T);

// Computes the CDF of the mixture indicators
arma::vec find_mixture_indicator_cdf(
    const arma::vec& datanorm);

}  // END namespace fast_sv

namespace general_sv {

// Computes the log-posterior of the latent
// vector in the exact SV model
double h_log_posterior(
    const arma::vec& h,
    const arma::vec& y,
    const double phi,
    const double rho,
    const double sigma,
    const double mu,
    const double h0);

// Computes the log-posterior of the latent
// vector in the auxiliary SV model
double h_aux_log_posterior(
    const arma::vec& h,
    const arma::vec& y_star,
    const arma::ivec& d,
    const double phi,
    const double rho,
    const double sigma,
    const double mu,
    const double h0);

}  // END namespace general_sv

}

#endif  // H_UTILS_H
