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
#'
#' @return Log returns of the (de-meaned / standardized) data.
#'
#' @family utilities
#'
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
#'
#' @return Called for its side effects. Returns argument \code{x} invisibly.
#'
#' @keywords hplot
#'
#' @family plotting
#'
#' @export

paratraceplot <- function(x, ...) {
 UseMethod("paratraceplot")
}
