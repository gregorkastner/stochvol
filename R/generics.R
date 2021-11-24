#  #####################################################################################
#  R package stochvol by
#     Gregor Kastner Copyright (C) 2013-2018
#     Gregor Kastner and Darjus Hosszejni Copyright (C) 2019-
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

#' Computes the Log Returns of a Time Series
#'
#' \code{logret} computes the log returns of a time
#' series, with optional de-meaning and/or standardization.
#'
#' @param dat The raw data.
#' @param demean Logical value indicating whether the data should
#' be de-meaned.
#' @param standardize Logical value indicating whether the data should
#' be standardized (in the sense that each component series has an empirical
#' variance equal to one).
#' @param ... Ignored.
#' @return Log returns of the (de-meaned / standardized) data.
#' @family utilities
#' @export
logret <- function(dat, demean = FALSE, standardize = FALSE, ...) {
 UseMethod("logret")
}


#' Trace Plot of MCMC Draws from the Parameter Posteriors
#'
#' Generic function for plotting iterations vs. sampled parameter values.
#' A detailed help for the method implemented in \pkg{stochvol} can be found in
#' \code{\link{paratraceplot.svdraws}}.
#'
#' @param x An object used to select a method.
#' @param ... Further arguments passed to or from other methods.
#' @return Called for its side effects. Returns argument \code{x} invisibly.
#' @keywords hplot
#' @family plotting
#' @export
paratraceplot <- function(x, ...) {
 UseMethod("paratraceplot")
}

