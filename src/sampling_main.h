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
 * sampling_main.h
 * 
 * Functions handling the data transformation from and to R,
 * the main Markov chain Monte Carlo (MCMC) loop, data storage additionally
 * to calling the "single update" functions.
 * 
 * The functions are separated into the fast_sv and the general_sv
 * variants due to the varying expert options.
 * 
 * These functions get exposed from the shared library and, most importantly,
 * to the R session. Detailed documentation of the parameters can be found
 * in the package manuals (which are generated from R/exports.R for these
 * functions). Type help("svsample_cpp") into your R session.
 */

#ifndef _SAMPLING_MAIN_H_
#define _SAMPLING_MAIN_H_

#include <RcppArmadillo.h>

namespace stochvol {

// Sample the latent states and the model parameters by using the fast_sv
// sampler. This is an extended version of the original algorithm by
// Kastner and Frühwirth-Schnatter (2014), and as such this covers the
// case of SV with heavy tails and without asymmetry (leverage).
// Additionally, the level of log variance can be fixed at zero.
Rcpp::List svsample_fast_cpp(
    const arma::vec& y_in,
    const int draws,
    const int burnin,
    const arma::mat& X_in,
    const Rcpp::List& priorspec_in,
    const int thinpara,
    const int thinlatent,
    const Rcpp::CharacterVector& keeptime_in,
    const Rcpp::List& startpara,
    const arma::vec& startlatent,
    const bool keeptau,
    const bool quiet,
    const bool correct_model_specification,
    const bool interweave,
    const double offset,
    const Rcpp::List& expert);

// Sample the latent states and the model parameters by using the general_sv
// sampler. This algorithm can sample from SV models with asymmetry (leverage),
// and handles a larger class of prior distributions. The sampler employs
// auxiliary mixture sampling for the latent states just like the fast_sv sampler
// but can correct the draws using a Metropolis-Hastings acceptance-rejection step.
// In our experience, this correction is more of a theoretical subtlety than an
// actual practical need.
// For the model parameters, the sampler employs an adaptive random walk
// Metropolis-Hastings (RWMH) scheme or an adaptive Metropolis adjusted Langevin (MALA)
// scheme.
Rcpp::List svsample_general_cpp(
    const arma::vec& y_in,
    const int draws,
    const int burnin,
    const arma::mat& X_in,
    const Rcpp::List& priorspec_in,
    const int thinpara,
    const int thinlatent,
    const Rcpp::CharacterVector& keeptime_in,
    const Rcpp::List& startpara,
    const arma::vec& startlatent,
    const bool keeptau,
    const bool quiet,
    const bool correct_model_specification,
    const bool interweave,
    const double offset,
    const Rcpp::List& expert);

}

#endif  // _SAMPLING_MAIN_H_

