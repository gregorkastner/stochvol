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
 * single_update.h
 * 
 * Functions that implement single updates of the "parameter groups."
 * The parameters are grouped into
 *   - SV parameters in
 *     . Fast SV (mu, phi, sigma, r, and h), and
 *     . General SV (mu, phi, sigma, rho, and h),
 *   - t-error parameter (nu), and
 *   - regression parameters (betas).
 * 
 * The functions for the SV parameters are separated into the fast_sv
 * and the general_sv variants due to the varying expert options.
 * 
 * For backwards compatibility, update_sv is defined and exported. The
 * function is deprecated and will be removed in a future release.
 * 
 * This documentation is reachable from R by typing help("stochvol_cpp").
 */

#ifndef _SINGLE_UPDATES_H_
#define _SINGLE_UPDATES_H_

#include <RcppArmadillo.h>
#include <adaptation.hpp>
#include "type_definitions.h"

namespace stochvol {

// Single MCMC update using fast SV
// 
// Samples the mixture indicators, the latent variables, and the model independent
// parameters mu, phi, and sigma. The input is the logarithm of the squared de-meaned
// observations. An approximate SV model is estimated instead of the exact SV model
// by the use of auxiliary mixture sampling.
// Depending on the prior specification, mu might not be updated.
// Depending on the expert settings, the function might follow the ancillarity-sufficiency
// interweaving strategy (ASIS, Yu and Meng, 2011) for sampling mu, phi, and sigma.
// Furthermore, the user can turn off the sampling of the parameters, the latents, or the
// mixture indicators in the expert settings.
// 
// @param log_data2: log(data^2), where data is the vector of de-meaned observations
// @param mu: parameter mu. Level of the latent process h. Updated in place
// @param phi: parameter phi, persistence of the latent process h. Updated in place
// @param sigma: parameter sigma, volatility of the latent process h, also called volvol. Updated in place
// @param h0: parameter h0, the initial value of the latent process h. Updated in place
// @param h: the vector of the latent process. Updated in place
// @param r: the vector of the mixture indicators. Updated in place
// @param prior_spec: prior specification object. See type_definitions.h
// @param expert: expert settings for this function. See type_definitions.h
void update_fast_sv(
    const arma::vec& log_data2,
    double& mu,
    double& phi,
    double& sigma,
    double& h0,
    arma::vec& h,
    arma::uvec& r,
    const PriorSpec& prior_spec,  // old parameters: C0, cT, Bsigma, a0, b0, bmu, Bmu, Gammaprior, truncnormal, dontupdatemu, priorlatent0 feed into this
    const ExpertSpec_FastSV& expert);  // old parameters: parameterization, centered_baseline, B011inv, B022inv, Mhcontrol, MHsteps feed into this

// Single MCMC update Student's t-distribution
// 
// Samples the degrees of freedom parameter of de-meaned and homoskedastic
// t-distributed input variates. Marginal data augmentation (MDA) is applied, tau
// is the vector of auxiliary latent states.
// Depending on the prior specification, nu might not be updated, just tau.
// 
// @param homosked_data: de-meaned and homoskedastic observations
// @param tau: the vector of the latent states used in MDA. Updated in place
// @param nu: parameter nu. The degrees of freedom for the t-distribution. Updated in place
// @param prior_spec: prior specification object. See type_definitions.h
void update_t_error(
    const arma::vec& homosked_data,
    arma::vec& tau,
    double& nu,
    const PriorSpec& prior_spec);

// Single MCMC update using general SV
// 
// Samples the latent variables and the model independent parameters mu, phi, sigma,
// and rho. The observations need to be provided in different formats for efficiency.
// An approximate SV model is as the default posterior distribution for the latent vector; however,
// there is the option to correct for model misspecification through the expert settings.
// Depending on the prior specification, some of mu, phi, sigma, and rho might not be updated.
// Depending on the expert settings, the function might follow the ancillarity-sufficiency
// interweaving strategy (ASIS, Yu and Meng, 2011) for sampling mu, phi, sigma, and rho.
// Also controlled by the expert settings, 
// Furthermore, the user can turn off the sampling of the parameters, the latents, or the
// mixture indicators in the expert settings.
// 
// @param data: the vector of de-meaned observations
// @param log_data2: log(data^2), where data is the vector of de-meaned observations
// @param sign_data: the sign of the data
// @param mu: parameter mu. Level of the latent process h. Updated in place
// @param phi: parameter phi, persistence of the latent process h. Updated in place
// @param sigma: parameter sigma, volatility of the latent process h, also called volvol. Updated in place
// @param rho: parameter rho. Accounts for asymmetry/the leverage effect. Updated in place
// @param h0: parameter h0, the initial value of the latent process h. Updated in place
// @param h: the vector of the latent process. Updated in place
// @param adaptation: object implementing the adaptive Metropolis-Hastings scheme. Updated in place. See adaptation.hpp
// @param prior_spec: prior specification object. See type_definitions.h
// @param expert: expert settings for this function. See type_definitions.h
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
    AdaptationCollection& adaptation,
    const PriorSpec& prior_spec,  // old parameters: prior_mu, prior_phi, prior_sigma2, prior_rho, gammaprior, dontupdatemu feed into this (plus priorlatent0, truncnormal nyi)
    const ExpertSpec_GeneralSV& expert);  // old parameters: strategy, correct, use_mala feed into this

// Single MCMC update of weighted robust (t-distributed) Bayesian regression
// 
// Samples the coefficients of a linear regression with known weights for the errors.
// The errors might be t-distributed with known degrees of freedom. Marginal data
// augmentation (MDA) is applied in that case, and tau is the vector of auxiliary
// latent states.
// 
// @param dependent_variable: the left hand side
// @param independent_variables: the matrix of the independent variables. Has to be of same height as the length of the dependent variable
// @param beta: the vector of the latent states used in MDA. Updated in place
// @param tau: the vector of the latent states used in MDA. Updated in place
// @param normalizer: the reciprocal of the error weights. E.g. the reciprocal of the known standard deviations given df.
// @param df: the degrees of freedom for the t-distribution. If large, then a standard/normal linear regression is conducted
// @param prior_spec: prior specification object. See type_definitions.h
void update_regressors(
    arma::vec dependent_variable,  // by value on purpose
    arma::mat independent_variables,  // by value on purpose
    arma::vec& beta,
    arma::vec& tau,
    const arma::vec& normalizer,  // inverse of heteroskedastic scales
    const double df,
    const PriorSpec& prior_spec);

// OLD FUNCTIONS

// Deprecated
inline
void update_sv(
    const arma::vec& data,
    arma::vec& curpara,
    arma::vec& h,
    double& h0,
    arma::vec& mixprob,
    arma::ivec& r,
    const bool centered_baseline,
    const double C0,
    const double cT,
    const double Bsigma,
    const double a0,
    const double b0,
    const double bmu,
    const double Bmu,
    const double B011inv,
    const double B022inv,
    const bool Gammaprior,
    const bool truncnormal,
    const double MHcontrol,
    const int MHsteps,
    const int parameterization,
    const bool dontupdatemu,
    const double priorlatent0) {
  double mu = curpara[0],
         phi = curpara[1],
         sigma2 = std::pow(curpara[2], 2);
  const PriorSpec prior_spec {
    (priorlatent0 <= 0) ? PriorSpec::Latent0() : PriorSpec::Latent0(PriorSpec::Constant(priorlatent0)),
    dontupdatemu ? PriorSpec::Mu(PriorSpec::Constant(0)) : PriorSpec::Mu(PriorSpec::Normal(bmu, std::sqrt(Bmu))),
    PriorSpec::Phi(PriorSpec::Beta(a0, b0)),
    Gammaprior ? PriorSpec::Sigma2(PriorSpec::Gamma(0.5, 0.5 / Bsigma)) : PriorSpec::Sigma2(PriorSpec::InverseGamma(2.5, C0))
  };
  const ExpertSpec_FastSV expert {
    parameterization > 2,  // interweave
    parameterization % 2 ? Parameterization::CENTERED : Parameterization::NONCENTERED,  // centered_baseline
    B011inv,
    B022inv,
    MHsteps,
    MHcontrol < 0 ? ExpertSpec_FastSV::ProposalSigma2::INDEPENDENCE : ExpertSpec_FastSV::ProposalSigma2::LOG_RANDOM_WALK,
    MHcontrol,
    truncnormal ? ExpertSpec_FastSV::ProposalPhi::REPEATED_ACCEPT_REJECT_NORMAL : ExpertSpec_FastSV::ProposalPhi::IMMEDIATE_ACCEPT_REJECT_NORMAL
  };
  arma::uvec r_u;
  std::copy(r.cbegin(), r.cend(), r_u.begin());
  update_fast_sv(data, mu, phi, sigma2, h0, h, r_u, prior_spec, expert);
  curpara = {mu, phi, std::sqrt(sigma2)};
}

}

#endif  // _SINGLE_UPDATES_H_
