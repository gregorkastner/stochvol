#  #####################################################################################
#  R package stochvol by
#     Gregor Kastner Copyright (C) 2013-2020
#     Darjus Hosszejni Copyright (C) 2019-2020
#  
#  This file is part of the R package stochvol: Efficient Bayesian
#  Inference for Stochastic Volatility Models.
#  
#  The R package stochvol is free software: you can redistribute it
#  and/or modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation, either version 2 or
#  any later version of the License.
#  
#  The R package stochvol is distributed in the hope that it will be
#  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#  General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with the R package stochvol. If that is not the case, please
#  refer to <http://www.gnu.org/licenses/>.
#  #####################################################################################

#' Single MCMC Update Using Fast SV
#' 
#' Samples the mixture indicators, the latent variables, and the model independent
#' parameters mu, phi, and sigma. The input is the logarithm of the squared de-meaned
#' observations. An approximate SV model is estimated instead of the exact SV model
#' by the use of auxiliary mixture sampling.
#' Depending on the prior specification, mu might not be updated.
#' Depending on the expert settings, the function might follow the ancillarity-sufficiency
#' interweaving strategy (ASIS, Yu and Meng, 2011) for sampling mu, phi, and sigma.
#' Furthermore, the user can turn off the sampling of the parameters, the latents, or the
#' mixture indicators in the expert settings.
#' 
#' @param log_data2 log(data^2), where data is the vector of de-meaned observations
#' @param mu parameter mu. Level of the latent process h. Updated in place
#' @param phi parameter phi, persistence of the latent process h. Updated in place
#' @param sigma parameter sigma, volatility of the latent process h, also called volvol. Updated in place
#' @param h0 parameter h0, the initial value of the latent process h. Updated in place
#' @param h the vector of the latent process. Updated in place
#' @param r the vector of the mixture indicators. Updated in place
#' @param prior_spec prior specification object. See type_definitions.h
#' @param expert expert settings for this function. See type_definitions.h
#' @family stochvol_cpp
#' @keywords update
#' @export
update_fast_sv <- function(
  log_data2,
  mu,
  phi,
  sigma,
  h0,
  h,
  r,
  prior_spec,
  expert) {
  stop("This is a dummy R function showcasing a C++ function in stochvol.")
}

#' Single MCMC Update of Student's t-Distribution
#' 
#' Samples the degrees of freedom parameter of de-meaned and homoskedastic
#' t-distributed input variates. Marginal data augmentation (MDA) is applied, tau
#' is the vector of auxiliary latent states.
#' Depending on the prior specification, nu might not be updated, just tau.
#' 
#' @param homosked_data de-meaned and homoskedastic observations
#' @param tau the vector of the latent states used in MDA. Updated in place
#' @param nu parameter nu. The degrees of freedom for the t-distribution. Updated in place
#' @param prior_spec prior specification object. See type_definitions.h
#' @family stochvol_cpp
#' @keywords update
#' @export
update_t_error <- function(
  homosked_data,
  tau,
  nu,
  prior_spec) {
  stop("This is a dummy R function showcasing a C++ function in stochvol.")
}

#' Single MCMC Update Using General SV
#' 
#' Samples the latent variables and the model independent parameters mu, phi, sigma,
#' and rho. The observations need to be provided in different formats for efficiency.
#' An approximate SV model is as the default posterior distribution for the latent vector; however,
#' there is the option to correct for model misspecification through the expert settings.
#' Depending on the prior specification, some of mu, phi, sigma, and rho might not be updated.
#' Depending on the expert settings, the function might follow the ancillarity-sufficiency
#' interweaving strategy (ASIS, Yu and Meng, 2011) for sampling mu, phi, sigma, and rho.
#' Also controlled by the expert settings, 
#' Furthermore, the user can turn off the sampling of the parameters, the latents, or the
#' mixture indicators in the expert settings.
#' 
#' @param data the vector of de-meaned observations
#' @param log_data2 log(data^2), where data is the vector of de-meaned observations
#' @param sign_data the sign of the data
#' @param mu parameter mu. Level of the latent process h. Updated in place
#' @param phi parameter phi, persistence of the latent process h. Updated in place
#' @param sigma parameter sigma, volatility of the latent process h, also called volvol. Updated in place
#' @param rho parameter rho. Accounts for asymmetry/the leverage effect. Updated in place
#' @param h0 parameter h0, the initial value of the latent process h. Updated in place
#' @param h the vector of the latent process. Updated in place
#' @param adaptation object implementing the adaptive Metropolis-Hastings scheme. Updated in place. See adaptation.hpp
#' @param prior_spec prior specification object. See type_definitions.h
#' @param expert expert settings for this function. See type_definitions.h
#' @family stochvol_cpp
#' @keywords update
#' @export
update_general_sv <- function(
  data,
  log_data2,
  sign_data,
  mu,
  phi,
  sigma,
  rho,
  h0,
  h,
  adaptation,
  prior_spec,
  expert) {
  stop("This is a dummy R function showcasing a C++ function in stochvol.")
}

#' Single MCMC Update of Weighted Robust Bayesian Regression
#' 
#' Samples the coefficients of a linear regression with known weights for the errors.
#' The errors might be t-distributed with known degrees of freedom; hence robust. Marginal data
#' augmentation (MDA) is applied in that case, and tau is the vector of auxiliary
#' latent states.
#' 
#' @param dependent_variable the left hand side
#' @param independent_variables the matrix of the independent variables. Has to be of same height as the length of the dependent variable
#' @param beta the vector of the latent states used in MDA. Updated in place
#' @param tau the vector of the latent states used in MDA. Updated in place
#' @param normalizer the reciprocal of the error weights. E.g. the reciprocal of the known standard deviations given df.
#' @param df the degrees of freedom for the t-distribution. If large, then a standard/normal linear regression is conducted
#' @param prior_spec prior specification object. See type_definitions.h
#' @family stochvol_cpp
#' @keywords update
#' @export
update_regressors <- function(
  dependent_variable,
  independent_variables,
  beta,
  tau,
  normalizer,
  df,
  prior_spec) {
  stop("This is a dummy R function showcasing a C++ function in stochvol.")
}
