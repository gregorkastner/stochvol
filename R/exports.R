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

#' Bindings to \code{C++} Functions in \code{stochvol}
#' 
#' All the heavy lifting in \code{stochvol} is implemented in \code{C++}
#' with the help of \code{R} packages \code{Rcpp} and \code{RcppArmadillo}.
#' These functions call the MCMC samplers in \code{C++} directly without any
#' any validation and transformations, expert use only!
#' 
#' The sampling functions are separated into fast SV and general SV. See more details
#' in the sections below.
#' 
#' @param y numeric vector of the observations
#' @param draws single positive integer, the number of draws to
#' return (after the burn-in)
#' @param burnin single positive integer, length of warm-up
#' period, this number of draws are discarded from the beginning
#' @param designmatrix numeric matrix of covariates. Dimensions:
#' \code{length(y)} times the number of covariates. If there are
#' no covariates then this should be \code{matrix(NA)}
#' @param priorspec a \code{priorspec} object created by
#' \link{specify_priors}
#' @param thinpara single number greater or equal to 1, coercible to integer.
#' Every \code{thinpara}th parameter draw is kept and returned. The default
#' value is 1, corresponding to no thinning of the parameter draws i.e. every
#' draw is stored.
#' @param thinlatent single number greater or equal to 1, coercible to integer.
#' Every \code{thinlatent}th latent variable draw is kept and returned. The
#' default value is 1, corresponding to no thinning of the latent variable
#' draws, i.e. every draw is kept.
#' @param keeptime Either 'all' (the default) or 'last'. Indicates which latent
#' volatility draws should be stored.
#' @param startpara named list, containing the starting values
#' for the parameter draws. It must contain
#' elements
#' \itemize{
#'   \item{mu}{an arbitrary numerical value}
#'   \item{phi}{real number between \code{-1} and \code{1}}
#'   \item{sigma}{a positive real number}
#'   \item{nu}{a number larger than \code{2}; can be \code{Inf}}
#'   \item{rho}{real number between \code{-1} and \code{1}}
#'   \item{beta}{a numeric vector of the same length as the number of covariates}
#'   \item{latent0}{a single number, the initial value for \code{h0}}}
#' @param startlatent vector of length \code{length(y)},
#' containing the starting values for the latent log-volatility draws.
#' @param keeptau Logical value indicating whether the 'variance inflation
#' factors' should be stored (used for the sampler with conditional t
#' innovations only). This may be useful to check at what point(s) in time the
#' normal disturbance had to be 'upscaled' by a mixture factor and when the
#' series behaved 'normally'.
#' @param print_settings List of three elements:
#'  \itemize{
#'    \item{quiet}{logical value indicating whether the progress bar and other
#' informative output during sampling should be omitted}
#'    \item{n_chains}{number of independent MCMC chains}
#'    \item{chain}{index of this chain}}
#' Please note that this function does not run multiple independent chains
#' but \code{svsample} offers different printing functionality depending on
#' whether it is executed as part of several MCMC chains in parallel
#' (chain specific messages) or simply as a single chain (progress bar).
#' @param correct_model_misspecification Logical value. If \code{FALSE},
#' then auxiliary mixture sampling is used to sample the latent
#' states. If \code{TRUE}, extra computations are made to correct for model
#' misspecification either ex-post by reweighting or on-line using a
#' Metropolis-Hastings step.
#' @param interweave Logical value. If \code{TRUE},
#' then ancillarity-sufficiency interweaving strategy (ASIS) is applied
#' to improve on the sampling efficiency for the parameters.
#' Otherwise one parameterization is used.
#' @param myoffset Single non-negative number that is used in
#' \code{log(y^2 + myoffset)} to prevent \code{-Inf} values in the auxiliary
#' mixture sampling scheme.
#' @param fast_sv named list of expert settings. We recommend the use of \code{default_fast_sv}.
#' @param general_sv named list of expert settings. We recommend the use of \code{default_general_sv}.
#' @section Fast SV:
#' Fast SV was developed in Kastner and Fruehwirth-Schnatter (2014). Fast SV estimates an
#' approximate SV model without leverage, where the approximation comes in through
#' auxiliary mixture approximations to the exact SV model. The sampler uses
#' the ancillarity-sufficiency interweaving strategy (ASIS) to improve on the sampling
#' efficiency of the model parameters, and it employs all-without-a-loop (AWOL)
#' for computationally efficient Kalman filtering of the conditionally Gaussian state space.
#' Correction for model misspecification happens as a post-processing step.
#' 
#' Fast SV employs sampling strategies that have been fine-tuned and specified for
#' vanilla SV (no leverage), and hence it can be fast and efficient but also more limited
#' in its feature set. The conditions for the fast SV sampler: \code{rho == 0}; \code{mu}
#' has either a normal prior or it is also constant \code{0}; the prior for \code{phi}
#' is a beta distribution; the prior for \code{sigma^2} is either a gamma distribution
#' with shape \code{0.5} or a mean- and variance-matched inverse gamma distribution;
#' either \code{keeptime == 'all'} or \code{correct_model_misspecification == FALSE}.
#' These criteria are NOT VALIDATED by fast SV on the \code{C++} level!
#' @section General SV:
#' General SV also estimates an
#' approximate SV model without leverage, where the approximation comes in through
#' auxiliary mixture approximations to the exact SV model. The sampler uses
#' both ASIS and AWOL.
#' 
#' General SV employs adapted random walk Metropolis-Hastings as the proposal for
#' the parameters \code{mu}, \code{phi}, \code{sigma}, and \code{rho}. Therefore,
#' more general prior distributions are allowed in this case.
#' @example inst/examples/svsample_cpp.R
#' @rdname svsample_cpp
#' @export
svsample_fast_cpp <- function(y, draws = 1, burnin = 0, designmatrix = matrix(NA), priorspec = specify_priors(), thinpara = 1, thinlatent = 1, keeptime = "all", startpara, startlatent, keeptau = !inherits(priorspec$nu, "sv_infinity"), print_settings = list(quiet = TRUE, n_chains = 1, chain = 1), correct_model_misspecification = FALSE, interweave = TRUE, myoffset = 0, fast_sv = default_fast_sv) {
    .Call(`_stochvol_svsample_fast_cpp`, y, draws, burnin, designmatrix, priorspec, thinpara, thinlatent, keeptime, startpara, startlatent, keeptau, print_settings, correct_model_misspecification, interweave, myoffset, fast_sv, PACKAGE = "stochvol")
}
#' @rdname svsample_cpp
#' @export
svsample_general_cpp <- function(y, draws = 1, burnin = 0, designmatrix = matrix(NA), priorspec = specify_priors(), thinpara = 1, thinlatent = 1, keeptime = "all", startpara, startlatent, keeptau = !inherits(priorspec$nu, "sv_infinity"), print_settings = list(quiet = TRUE, n_chains = 1, chain = 1), correct_model_misspecification = FALSE, interweave = TRUE, myoffset = 0, general_sv = default_general_sv) {
    .Call(`_stochvol_svsample_general_cpp`, y, draws, burnin, designmatrix, priorspec, thinpara, thinlatent, keeptime, startpara, startlatent, keeptau, print_settings, correct_model_misspecification, interweave, myoffset, general_sv, PACKAGE = "stochvol")
}

get_omori_constants <- function() {
    .Call(`_stochvol_get_omori_constants`, PACKAGE = "stochvol")
}

# Register entry points for exported C++ functions
methods::setLoadAction(function(ns) {
    .Call('_stochvol_Export_registerCCallable', PACKAGE = 'stochvol')
})

