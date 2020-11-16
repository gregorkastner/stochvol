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
 * sampling_latent_states.h
 * 
 * Functions related to the sampling of the latent states.
 * This includes the all-without-a-loop (AWOL, McCausland et al., 2011)
 * sampling of the latent vector in the auxiliary model from Omori et al. (2007)
 * and the corresponding correction step thereafter in the general SV sampler.
 */

#ifndef _SAMPLING_LATENT_STATES_H_
#define _SAMPLING_LATENT_STATES_H_

#include <RcppArmadillo.h>
#include <expert.hpp>

namespace stochvol {

// Store and encapsulate an entire latent vector in one structure.
struct LatentVector {
  double h0;
  arma::vec h;
};

namespace fast_sv {

// Sample from the conditional posterior distribution
// p(latent_states | data, parameters). Employ auxiliary mixture
// sampling (by Omori et al 2007) in an efficient manner all-without-a-loop
// (AWOL, by McCausland 2011).
LatentVector draw_latent(
    const arma::vec& data,
    const double mu,
    const double phi,
    const double sigma,
    const arma::uvec& r,
    const PriorSpec& prior_spec,
    const ExpertSpec_FastSV& expert);

// Sample the vector of mixture indicators. For more
// information see Omori et al (2007).
arma::uvec draw_mixture_indicators(
    const arma::vec& data,
    const double mu,
    const double phi,
    const double sigma,
    const arma::vec& h);

}  // END namespace fast_sv

namespace general_sv {

// Sample from the conditional posterior distribution
// p(latent_states | data, parameters). Employ auxiliary mixture
// sampling (by Omori et al 2007) as the proposal distribution and
// accept or reject the proposed vector according to the
// Metropolis-Hastings acceptance-rejection rate. The very first
// element h0 is drawn from its conditional posterior
// p(h0 | h1, parameters), which is a normal distribution.
// 
// If parameter 'correct' is false then the exact conditional
// posterior is replaced by the proposal distribution.
LatentVector draw_latent(
    const arma::vec& y,
    const arma::vec& y_star,
    const arma::ivec& d,
    const double mu,
    const double phi,
    const double sigma,
    const double rho,
    const arma::vec& h,
    const arma::vec& ht,
    const PriorSpec& prior_spec,
    const ExpertSpec_GeneralSV& expert);

}  // END namespace general_sv

}

#endif  // _SAMPLING_LATENT_STATES_H_

