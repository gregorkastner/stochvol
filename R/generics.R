#' Computes the log returns of a time series
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
