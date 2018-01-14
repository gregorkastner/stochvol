

#' Euro exchange rate data
#' 
#' The data set contains the daily bilateral prices of one Euro in 23
#' currencies from January 3, 2000, until April 4, 2012. Conversions to New
#' Turkish Lira and Fourth Romanian Leu have been incorporated.
#' 
#' 
#' @name exrates
#' @docType data
#' @seealso \code{\link{svsample}}
#' @source ECB Statistical Data Warehouse (\url{http://sdw.ecb.europa.eu})
#' @keywords datasets
#' @examples
#' 
#' \dontrun{
#' data(exrates)
#' dat <- logret(exrates$USD, demean = TRUE)  ## de-meaned log-returns
#' res <- svsample(dat)                       ## run MCMC sampler
#' plot(res, forecast = 100)                  ## display results
#' }
#' 
NULL





#' Common Extractors for 'svdraws' Objects
#' 
#' Some simple extractors returning the corresponding element of an
#' \code{svdraws} object.
#' 
#' 
#' @name extractors
#' @aliases para latent latent0 priors thinning runtime
#' @param x \code{svdraws} object.
#' @return The return value depends on the actual funtion:
#' \item{para(x)}{extracts the parameter draws and returns them as an
#' \code{mcmc} object.} \item{latent(x)}{extracts the latent contemporaneous
#' log-volatility draws and returns them as an \code{mcmc} object.}
#' \item{latent0(x)}{extracts the latent initial log-volatility draws and
#' returns as an \code{mcmc} object.} \item{priors(x)}{extracts the prior
#' parameters used and returns them in a \code{list}.}
#' \item{thinning(x)}{extracts the thinning parameters used and returns them in
#' a \code{list}.} \item{runtime(x)}{extracts the runtime and returns it as a
#' \code{proc_time} object.}
#' @author Gregor Kastner \email{gregor.kastner@@wu.ac.at}
#' @keywords utilities
NULL





#' Efficient Bayesian Inference for Stochastic Volatility (SV) Models
#' 
#' This package provides an efficient algorithm for fully Bayesian estimation
#' of stochastic volatility (SV) models via Markov chain Monte Carlo (MCMC)
#' methods. Methodological details are given in Kastner and Frühwirth-Schnatter
#' (2014) <doi:10.1016/j.csda.2013.01.002>; the most common use cases are
#' described in Kastner (2016) <doi:10.18637/jss.v069.i05>.
#' 
#' Bayesian inference for stochastic volatility models using MCMC methods
#' highly depends on actual parameter values in terms of sampling efficiency.
#' While draws from the posterior utilizing the standard centered
#' parameterization break down when the volatility of volatility parameter in
#' the latent state equation is small, non-centered versions of the model show
#' deficiencies for highly persistent latent variable series. The novel
#' approach of ancillarity-sufficiency interweaving (Yu and Meng, 2011) has
#' recently been shown to aid in overcoming these issues for a broad class of
#' multilevel models. This package provides software for ``combining best of
#' different worlds'' which allows for inference for parameter constellations
#' that have previously been infeasible to estimate without the need to select
#' a particular parameterization beforehand.
#' 
#' @name stochvol-package
#' @aliases stochvol-package stochvol
#' @docType package
#' @useDynLib stochvol, .registration = TRUE
#' @note This package is currently in active development. Your comments,
#' suggestions and requests are warmly welcome!
#' @author Gregor Kastner \email{gregor.kastner@@wu.ac.at}
#' @references Kastner, G. and Frühwirth-Schnatter, S. (2014).
#' Ancillarity-sufficiency interweaving strategy (ASIS) for boosting MCMC
#' estimation of stochastic volatility models. \emph{Computational Statistics &
#' Data Analysis}, \bold{76}, 408--423,
#' \url{http://dx.doi.org/10.1016/j.csda.2013.01.002}.
#' 
#' Yu, Y. and Meng, X.-L. (2011). To Center or Not to Center: That is Not the
#' Question---An Ancillarity-Suffiency Interweaving Strategy (ASIS) for
#' Boosting MCMC Efficiency. \emph{Journal of Computational and Graphical
#' Statistics}, \bold{20}, 531--570,
#' \url{http://dx.doi.org/10.1198/jcgs.2011.203main}.
#' @keywords package models ts
#' @examples
#' 
#' ## Simulate a highly persistent SV process 
#' sim <- svsim(500, mu = -10, phi = 0.99, sigma = 0.2)
#' 
#' ## Obtain 4000 draws from the sampler (that's too little!)
#' draws <- svsample(sim$y, draws = 4000, burnin = 100, priormu = c(-10, 1),
#'                   priorphi = c(20, 1.2), priorsigma = 0.2)
#' 
#' ## Predict 20 days ahead
#' fore <- predict(draws, 20)
#' 
#' ## plot the results
#' plot(draws, forecast = fore)
#' 
NULL



