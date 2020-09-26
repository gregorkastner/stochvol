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
 * sampling_parameters.h
 * 
 * Functions for sampling the time independent model parameters: mu, phi,
 * sigma, and rho.
 * 
 * The functions are separated into the fast_sv and the general_sv
 * variants due to the varying expert options.
 */

#ifndef _SAMPLING_PARAMETERS_H_
#define _SAMPLING_PARAMETERS_H_

#include <RcppArmadillo.h>
#include <adaptation.hpp>
#include "utils_parameters.h"
#include "type_definitions.h"

namespace stochvol {

namespace fast_sv {

// Draw mu, phi, and sigma.
// In the centered parameterization, all three parameters are
// sampled conditional on the latent states (they're independent
// the data given the latent vector).
// In the noncentered parameterization, mu and sigma are sampled
// conditional on the data and the latents while phi is sampled
// again only conditional on the latent vector.
// 
// For a detailed explanation, see Kastner and Frühwirth-Schnatter (2014).
SampledTheta draw_theta(
    const arma::vec& log_data2,
    const double mu,
    const double phi,
    const double sigma,
    const double h0,
    const double ht0,
    const arma::vec& h,
    const arma::vec& ht,
    const arma::uvec& r,
    const PriorSpec& prior_spec,
    const ExpertSpec_FastSV& expert,
    const Parameterization parameterization);

}  // END namespace fast_sv

namespace general_sv {

// Draw mu, phi, sigma, and rho.
// Both in the centered and the noncentered parameterization
// all parameters are sampled together according to a random
// walk Metropolis-Hastings algorithm. The tuning parameters
// of the random walk are passed in the diffusion_ken object.
// Input parameter update_indicator controls which of the four
// model parameters are to be updated.
// Finally, exp_h_half and exp_h_half_proposal_nc serve
// computational efficiency purposes, they are cached versions
// the current and the proposed exp(h/2).
SampledTheta draw_theta(
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
    arma::vec& exp_h_half_proposal_nc,
    const arma::uvec& update_indicator,
    const PriorSpec& prior_spec,
    const ExpertSpec_GeneralSV& expert,
    const ProposalDiffusionKen& diffusion_ken,
    const Parameterization parameterization);

}  // END namespace general_sv

}

#endif  // _SAMPLING_PARAMETERS_H_

