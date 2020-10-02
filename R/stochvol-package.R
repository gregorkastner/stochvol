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
#' \code{mcmc} object.}
#' \item{latent(x)}{extracts the latent contemporaneous
#' log-volatility draws and returns them as an \code{mcmc} object.}
#' \item{latent0(x)}{extracts the latent initial log-volatility draws and
#' returns as an \code{mcmc} object.}
#' \item{priors(x)}{extracts the prior
#' parameters used and returns them in a \code{list}.}
#' \item{thinning(x)}{extracts the thinning parameters used and returns them in
#' a \code{list}.}
#' \item{runtime(x)}{extracts the runtime and returns it as a
#' \code{proc_time} object.}
#' @keywords utilities
NULL

#' Efficient Bayesian Inference for Stochastic Volatility (SV) Models
#' 
#' This package provides an efficient algorithm for fully Bayesian estimation
#' of stochastic volatility (SV) models via Markov chain Monte Carlo (MCMC)
#' methods. Methodological details are given in Kastner and Frühwirth-Schnatter
#' (2014); the most common use cases are described in Kastner (2016). Recently,
#' the package has been extended to allow for the leverage effect.
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
#' @importFrom utils tail head flush.console
#' @importFrom graphics plot par hist mtext lines title matplot points abline layout plot.default axis boxplot
#' @importFrom stats cov rt rgamma rnorm sd IQR density time lowess dnorm dbeta dgamma dexp qnorm qt ppoints ts.plot median quantile predict plot.ts qqline qqnorm qqplot
#' @importFrom coda mcmc nvar niter varnames traceplot mcmc.list nvar nchain effectiveSize mcpar
#' @importFrom Rcpp sourceCpp
#' @note This package is currently in active development. Your comments,
#' suggestions and requests are warmly welcome!
#' @author Gregor Kastner \email{gregor.kastner@@wu.ac.at}, Darjus Hosszejni \email{darjus.hosszejni@@wu.ac.at}
#' @references Kastner, G. and Frühwirth-Schnatter, S. (2014).
#' Ancillarity-Sufficiency Interweaving Strategy (ASIS) for Boosting MCMC
#' Estimation of Stochastic Volatility Models. \emph{Computational Statistics &
#' Data Analysis}, \bold{76}, 408--423,
#' \url{http://dx.doi.org/10.1016/j.csda.2013.01.002}.
#'
#' Kastner, G. (2016). Dealing with Stochastic Volatility in Time Series Using the R Package stochvol.
#' \emph{Journal of Statistical Software}, \bold{69}(5), 1--30,
#' \url{http://dx.doi.org/10.18637/jss.v069.i05}.
#' 
#' Yu, Y. and Meng, X.-L. (2011). To Center or Not to Center: That is Not the
#' Question---An Ancillarity-Suffiency Interweaving Strategy (ASIS) for
#' Boosting MCMC Efficiency. \emph{Journal of Computational and Graphical
#' Statistics}, \bold{20}(3), 531--570,
#' \url{http://dx.doi.org/10.1198/jcgs.2011.203main}.
#' @keywords package models ts
#' @example inst/examples/stochvol-package.R
NULL

