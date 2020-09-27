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

# R wrapper functions for the main MCMC loop

#' Markov Chain Monte Carlo (MCMC) Sampling for the Stochastic Volatility (SV)
#' Model
#' 
#' \code{svsample} simulates from the joint posterior distribution of the SV
#' parameters \code{mu}, \code{phi}, \code{sigma} (and potentially \code{nu}),
#' along with the latent log-volatilities \code{h_0,...,h_n} and returns the
#' MCMC draws. If a design matrix is provided, simple Bayesian regression can
#' also be conducted.
#' 
#' For details concerning the algorithm please see the paper by Kastner and
#' Frühwirth-Schnatter (2014).
#' 
#' @param y numeric vector containing the data (usually log-returns), which
#' must not contain zeros. Alternatively, \code{y} can be an \code{svsim}
#' object. In this case, the returns will be extracted and a message is signalled.
#' @param draws single number greater or equal to 1, indicating the number of
#' draws after burn-in (see below). Will be automatically coerced to integer.
#' The default value is 10000.
#' @param burnin single number greater or equal to 0, indicating the number of
#' draws discarded as burn-in. Will be automatically coerced to integer. The
#' default value is 1000.
#' @param designmatrix regression design matrix for modeling the mean. Must
#' have \code{length(y)} rows. Alternatively, \code{designmatrix} may be a
#' string of the form \code{"arX"}, where \code{X} is a nonnegative integer. To
#' fit a constant mean model, use \code{designmatrix = "ar0"} (which is
#' equivalent to \code{designmatrix = matrix(1, nrow = length(y))}). To fit an
#' AR(1) model, use \code{designmatrix = "ar1"}, and so on. If some elements of
#' \code{designmatrix} are \code{NA}, the mean is fixed to zero (pre-1.2.0
#' behavior of \pkg{stochvol}).
#' @param priormu numeric vector of length 2, indicating mean and standard
#' deviation for the Gaussian prior distribution of the parameter \code{mu},
#' the level of the log-volatility. The default value is \code{c(0, 100)},
#' which constitutes a practically uninformative prior for common exchange rate
#' datasets, stock returns and the like.
#' @param priorphi numeric vector of length 2, indicating the shape parameters
#' for the Beta prior distribution of the transformed parameter
#' \code{(phi + 1) / 2}, where \code{phi} denotes the persistence of the
#' log-volatility. The default value is \code{c(5, 1.5)}, which constitutes a
#' prior that puts some belief in a persistent log-volatility but also
#' encompasses the region where \code{phi} is around 0.
#' @param priorsigma single positive real number, which stands for the scaling
#' of the transformed parameter \code{sigma^2}, where \code{sigma} denotes the
#' volatility of log-volatility. More precisely, \code{sigma^2 ~ priorsigma *
#' chisq(df = 1)}. The default value is \code{1}, which constitutes a
#' reasonably vague prior for many common exchange rate datasets, stock returns
#' and the like.
#' @param priornu single non-negative number, indicating the rate parameter
#' for the exponential prior distribution of the parameter; can be \code{Inf}
#' \code{nu}, the degrees-of-freedom parameter of the conditional innovations
#' t-distribution. The default value is \code{0}, fixing the
#' degrees-of-freedom to infinity. This corresponds to conditional standard
#' normal innovations, the pre-1.1.0 behavior of \pkg{stochvol}.
#' @param priorbeta numeric vector of length 2, indicating the mean and
#' standard deviation of the Gaussian prior for the regression parameters. The
#' default value is \code{c(0, 10000)}, which constitutes a very vague prior
#' for many common datasets. Not used if \code{designmatrix} is \code{NA}.
#' @param priorlatent0 either a single non-negative number or the string
#' \code{'stationary'} (the default, also the behavior before version 1.3.0).
#' When \code{priorlatent0} is equal to \code{'stationary'}, the stationary
#' distribution of the latent AR(1)-process is used as the prior for the
#' initial log-volatility \code{h_0}. When \code{priorlatent0} is equal to a
#' number \eqn{B}, we have \eqn{h_0 \sim N(\mu, B\sigma^2)} a priori.
#' @param priorrho either \code{NA} for the no-leverage case or a numeric
#' vector of length 2 that specify the beta prior distribution for
#' \code{(rho+1)/2}
#' @param priorspec in case one needs different prior distributions than the
#' ones specified by \code{priormu}, \code{...}, \code{priorrho}, a \code{priorspec}
#' object can be supplied here. A smart constructor for this usecase is
#' \link{specify_priors}.
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
#' @param keeptau logical value indicating whether the 'variance inflation
#' factors' should be stored (used for the sampler with conditional t
#' innovations only). This may be useful to check at what point(s) in time the
#' normal disturbance had to be 'upscaled' by a mixture factor and when the
#' series behaved 'normally'.
#' @param quiet logical value indicating whether the progress bar and other
#' informative output during sampling should be omitted. The default value is
#' \code{FALSE}, implying verbose output.
#' @param startpara \emph{optional} named list, containing the starting values
#' for the parameter draws. If supplied, \code{startpara} must contain three
#' elements named \code{mu}, \code{phi}, and \code{sigma}, where \code{mu} is
#' an arbitrary numerical value, \code{phi} is a real number between \code{-1}
#' and \code{1}, and \code{sigma} is a positive real number. Moreover, if
#' \code{priornu} is not \code{0}, \code{startpara} must also contain an
#' element named \code{nu} (the degrees of freedom parameter for the
#' t-innovations). The default value is equal to the prior mean. 
#' @param startlatent \emph{optional} vector of length \code{length(y)},
#' containing the starting values for the latent log-volatility draws. The
#' default value is \code{rep(-10, length(y))}.
#' @param expert \emph{optional} named list of expert parameters. For most
#' applications, the default values probably work best. Interested users are
#' referred to the literature provided in the References section. If
#' \code{expert} is provided, it may contain the following named elements:
#' \code{interweave}: Logical value. If \code{TRUE} (the default),
#' then ancillarity-sufficiency interweaving strategy (ASIS) is applied
#' to improve on the sampling efficiency for the parameters.
#' Otherwise one parameterization is used.
#' \code{correct_model_misspecification}: Logical value. If \code{FALSE}
#' (the default), then auxiliary mixture sampling is used to sample the latent
#' states. If \code{TRUE}, extra computations are made to correct for model
#' misspecification either ex-post by reweighting or on-line using a
#' Metropolis-Hastings step.}
#' @param \dots Any extra arguments will be forwarded to
#' \code{\link{updatesummary}}, controlling the type of statistics calculated
#' for the posterior draws.
#' @return The value returned is a list object of class \code{svdraws} holding
#' \item{para}{\code{mcmc} object containing the \emph{parameter} draws from
#' the posterior distribution.}
#' \item{latent}{\code{mcmc} object containing the
#' \emph{latent instantaneous log-volatility} draws from the posterior
#' distribution.}
#' \item{latent0}{\code{mcmc} object containing the \emph{latent
#' initial log-volatility} draws from the posterior distribution.}
#' \item{tau}{\code{mcmc} object containing the \emph{latent variance inflation
#' factors} for the sampler with conditional t-innovations \emph{(optional)}.}
#' \item{beta}{\code{mcmc} object containing the \emph{regression coefficient}
#' draws from the posterior distribution \emph{(optional)}.}
#' \item{y}{the left hand side of the observation equation, usually
#' the argument \code{y}. In case of an AR(\code{k}) specification, the
#' first \code{k} elements are removed.}
#' \item{runtime}{\code{proc_time} object containing the
#' run time of the sampler.}
#' \item{priors}{a \code{priorspec} object containing the parameter
#' values of the prior distributions for \code{mu},
#' \code{phi}, \code{sigma}, \code{nu}, \code{rho}, and
#' \code{beta}s, and the variance of specification for \code{latent0}.}
#' \item{thinning}{\code{list} containing the thinning
#' parameters, i.e. the arguments \code{thinpara}, \code{thinlatent} and
#' \code{keeptime}.}
#' \item{summary}{\code{list} containing a collection of
#' summary statistics of the posterior draws for \code{para}, \code{latent},
#' and \code{latent0}.}
#' \item{meanmodel}{\code{character} containing information about how \code{designmatrix}
#' was employed.}
#' 
#' To display the output, use \code{print}, \code{summary} and \code{plot}. The
#' \code{print} method simply prints the posterior draws (which is very likely
#' a lot of output); the \code{summary} method displays the summary statistics
#' currently stored in the object; the \code{plot} method
#' \code{\link{plot.svdraws}} gives a graphical overview of the posterior
#' distribution by calling \code{\link{volplot}}, \code{\link{traceplot}} and
#' \code{\link{densplot}} and displaying the results on a single page.
#' @note If \code{y} contains zeros, you might want to consider de-meaning your
#' returns or use \code{designmatrix = "ar0"}.
#' @seealso \code{\link{svsim}}, \code{\link{specify_priors}}
#' @references Kastner, G. and Frühwirth-Schnatter, S. (2014).
#' Ancillarity-sufficiency interweaving strategy (ASIS) for boosting MCMC
#' estimation of stochastic volatility models. \emph{Computational Statistics &
#' Data Analysis}, \bold{76}, 408--423,
#' \url{http://dx.doi.org/10.1016/j.csda.2013.01.002}.
#' @keywords models ts
#' @examples
#' # Example 1
#' ## Simulate a short and highly persistent SV process 
#' sim <- svsim(100, mu = -10, phi = 0.99, sigma = 0.2)
#' 
#' ## Obtain 5000 draws from the sampler (that's not a lot)
#' draws <- svsample(sim$y, draws = 5000, burnin = 100,
#' 		  priormu = c(-10, 1), priorphi = c(20, 1.5), priorsigma = 0.2)
#' 
#' ## Check out the results
#' summary(draws)
#' plot(draws, simobj = sim)
#' 
#' 
#' # Example 2
#' ## Simulate a short and conditionally heavy-tailed SV process
#' sim <- svsim(100, mu = -10, phi = 0.96, sigma = 0.3, nu = 2.5)
#' 
#' ## Obtain 5000 draws from the sampler
#' tdraws <- svsample(sim$y, draws = 5000, burnin = 100,
#' 		  priormu = c(-10, 1), priorphi = c(20, 1.5), priorsigma = 0.5,
#'      priornu = 0.2)
#' 
#' ## Check out the results
#' summary(tdraws)
#' plot(tdraws, simobj = sim)
#' 
#' 
#' \dontrun{
#' # Example 3
#' ## AR(1) structure for the mean
#' data(exrates)
#' len <- 3000
#' ahead <- 100
#' y <- head(exrates$USD, len)
#' 
#' ## Fit AR(1)-SVL model to EUR-USD exchange rates
#' res <- svsample(y, designmatrix = "ar1")
#' 
#' ## Use predict.svdraws to obtain predictive distributions
#' preddraws <- predict(res, steps = ahead)
#' 
#' ## Calculate predictive quantiles
#' predquants <- apply(preddraws$y, 2, quantile, c(.1, .5, .9))
#' 
#' ## Visualize
#' expost <- tail(head(exrates$USD, len+ahead), ahead)
#' ts.plot(y, xlim = c(length(y)-4*ahead, length(y)+ahead),
#' 	       ylim = range(c(predquants, expost, tail(y, 4*ahead))))
#' for (i in 1:3) {
#'   lines((length(y)+1):(length(y)+ahead), predquants[i,],
#'         col = 3, lty = c(2, 1, 2)[i])
#' }
#' lines((length(y)+1):(length(y)+ahead), expost,
#'       col = 2)
#' 
#' 
#' # Example 4
#' ## Predicting USD based on JPY and GBP in the mean
#' data(exrates)
#' len <- 3000
#' ahead <- 30
#' ## Calculate log-returns
#' logreturns <- apply(exrates[, c("USD", "JPY", "GBP")], 2,
#'                     function (x) diff(log(x)))
#' logretUSD <- logreturns[2:(len+1), "USD"]
#' regressors <- cbind(1, as.matrix(logreturns[1:len, ]))  # lagged by 1 day
#' 
#' ## Fit SV model to EUR-USD exchange rates
#' res <- svsample(logretUSD, designmatrix = regressors)
#' 
#' ## Use predict.svdraws to obtain predictive distributions
#' predregressors <- cbind(1, as.matrix(logreturns[(len+1):(len+ahead), ]))
#' preddraws <- predict(res, steps = ahead,
#'                      newdata = predregressors)
#' predprice <- exrates[len+2, "USD"] * exp(t(apply(preddraws$y, 1, cumsum)))
#' 
#' ## Calculate predictive quantiles
#' predquants <- apply(predprice, 2, quantile, c(.1, .5, .9))
#' 
#' ## Visualize
#' priceUSD <- exrates[3:(len+2), "USD"]
#' expost <- exrates[(len+3):(len+ahead+2), "USD"]
#' ts.plot(priceUSD, xlim = c(len-4*ahead, len+ahead+1),
#' 	       ylim = range(c(expost, predquants, tail(priceUSD, 4*ahead))))
#' for (i in 1:3) {
#'   lines(len:(len+ahead), c(tail(priceUSD, 1), predquants[i,]),
#'         col = 3, lty = c(2, 1, 2)[i])
#' }
#' lines(len:(len+ahead), c(tail(priceUSD, 1), expost),
#'       col = 2)
#' }
#' @export
svsample <- function(y, draws = 10000, burnin = 1000, designmatrix = NA,
                     priormu = c(0, 100), priorphi = c(5, 1.5), priorsigma = 1,
                     priornu = 0, priorrho = NA,
                     priorbeta = c(0, 10000), priorlatent0 = "stationary",
                     priorspec = NULL,
                     thinpara = 1, thinlatent = 1, keeptime = "all", keeptau = FALSE,
                     quiet = FALSE, startpara = NULL, startlatent = NULL, expert = NULL, ...) {

  # Validation
  ## y
  if (inherits(y, "svsim")) {
    simobj <- y
    y <- simobj[["y"]]
    message("Extracted data vector from 'svsim'-object.")
  } else {
    simobj <- NULL
  }
  assert_numeric(y, "Argument 'y'")
  assert_gt(length(y), 1, "The length of the input time series 'y'")

  y_orig <- y
  y <- as.vector(y)

  myoffset <- if (any(is.na(designmatrix)) && any(y^2 == 0)) {
    warning("Argument 'y' (data vector) contains values very close to zero. I am applying an offset constant of size ", myoffset, " to do the auxiliary mixture sampling. If you want to avoid this, you might consider de-meaning the returns before calling this function.")
    sd(y)/10000
  } else {
    0
  }

  ## draws
  assert_numeric(draws, "Argument 'draws'")
  assert_single(draws, "Argument 'draws'")
  assert_positive(draws, "Argument 'draws'")
  draws <- as.integer(draws)

  ## burnin
  assert_numeric(burnin, "Argument 'burnin'")
  assert_single(burnin, "Argument 'burnin'")
  assert_ge(burnin, 0, "Argument 'burnin'")
  burnin <- as.integer(burnin)

  ## regression
  meanmodel <- "matrix"
  arorder <- 0L
  if (any(is.na(designmatrix))) {
    designmatrix <- matrix(NA_real_)
    meanmodel <- "none"
  } else {
    if (any(grep("ar[0-9]+$", as.character(designmatrix)[1]))) {
      arorder <- as.integer(gsub("ar", "", as.character(designmatrix)))
      if (length(y) <= arorder + 1) {
        stop("Time series 'y' is too short for this AR process.")
      }
      designmatrix <- matrix(rep.int(1, length(y) - arorder), ncol = 1)
      colnames(designmatrix) <- c("const")
      meanmodel <- "constant"
      if (arorder >= 1) {  # e.g. first row with "ar(3)": c(1, y_3, y_2, y_1)
        for (i in 1:arorder) {
          oldnames <- colnames(designmatrix)
          designmatrix <- cbind(designmatrix, y[(arorder-i+1):(length(y)-i)])
          colnames(designmatrix) <- c(oldnames, paste0("ar", i))
        }
        y <- y[-(1:arorder)]
        meanmodel <- paste0("ar", arorder)
      }
    } else if (is.character(designmatrix)) {
      stop("Argument 'designmatrix' must be a numeric matrix or an AR-specification.")
    }
    assert_numeric(designmatrix, "The processed argument 'designmatrix'")
    if (!is.matrix(designmatrix)) {
      designmatrix <- matrix(designmatrix, ncol = 1)
    }
    if (NROW(designmatrix) != length(y)) {
      stop("Number of columns of argument 'designmatrix' must be equal to length(y).")
    }
  }

  ## priors
  if (isTRUE(is.null(priorspec))) {
    validate_sv_priors(priormu, priorphi, priorsigma, priornu, priorrho, priorbeta, priorlatent0)
    priorspec <-
      specify_priors(mu = sv_normal(mean = priormu[1], sd = priormu[2]),
                     phi = sv_beta(alpha = priorphi[1], beta = priorphi[2]),
                     sigma2 = sv_gamma(shape = 0.5, rate = 0.5 / priorsigma),
                     nu = if (priornu == 0) sv_infinity() else sv_exponential(rate = priornu),
                     rho = if (isTRUE(is.na(priorrho))) sv_constant(value = 0) else sv_beta(alpha = priorrho[1], beta = priorrho[2]),
                     beta = sv_multinormal(mean = priorbeta[1], sd = priorbeta[2], dim = NCOL(designmatrix)),
                     latent0 = priorlatent0)
  } else if (!inherits(priorspec, "sv_priorspec")) {
    stop("Received argument 'priorspec' but it does not have the correct form. Please refer to the function called 'specify_priors'.")
  }

  ## thinning parameters
  validate_thinning(thinpara, thinlatent, keeptime)
  thinpara <- as.integer(thinpara)
  thinlatent <- as.integer(thinlatent)
  thintime <- switch(keeptime,
                     all = 1L,
                     last = length(y))

  ## expert
  expertdefault <-
    list(correct_model_misspecification = FALSE,  # online correction for general_sv and post-correction for fast_sv
         interweave = TRUE,
         fast_sv =  # UNDOCUMENTED; very expert settings of the fast_sv sampler
           default_fast_sv,
         general_sv =  # UNDOCUMENTED; very expert settings of the general_sv sampler
           default_general_sv)
  expert <- apply_default_list(expert, expertdefault, "Names in expert", "allowed names in expert")
  validate_expert(expert)
  correct_model_misspecification <- expert$correct_model_misspecification
  interweave <- expert$interweave
  fast_sv <- expert$fast_sv
  general_sv <- expert$general_sv

  # Initial values
  startparadefault <-
    list(mu = mean(priorspec$mu),  # init_mu
         phi = if (inherits(priorspec$phi, "sv_beta")) 2 * mean(priorspec$phi) - 1 else mean(priorspec$phi),
         sigma = sqrt(mean(priorspec$sigma2)),
         nu = 2 + mean(priorspec$nu),
         rho = if (inherits(priorspec$rho, "sv_beta")) 2 * mean(priorspec$rho) - 1 else mean(priorspec$rho),
         beta = rep.int(mean(priorspec$beta), NCOL(designmatrix)),
         latent0 = -10)
  startlatentdefault <- rep.int(-10, length(y))
  startpara <- apply_default_list(startpara, startparadefault)
  startlatent <- apply_default_list(startlatent, startlatentdefault)
  validate_initial_values(startpara, startlatent, y, designmatrix)

  # Decision about the sampler
  use_fast_sv <-
    # rho == 0
    (inherits(priorspec$rho, "sv_constant") && isTRUE(priorspec$rho$value == 0)) &&
    # mu is either 0 or normal
    ((inherits(priorspec$mu, "sv_constant") && isTRUE(priorspec$mu$value == 0)) ||
       inherits(priorspec$mu, "sv_normal")) &&
    # keeptime == "last" && correct_model_misspecification can't be done as of yet
    (keeptime == "all" || !correct_model_misspecification) &&
    # prior for phi is beta
    inherits(priorspec$phi, "sv_beta") &&
    # prior for sigma is gamma(0.5, _)
    (inherits(priorspec$sigma2, "sv_gamma") && priorspec$sigma2$shape == 0.5)

  # Call sampler
  myquiet <- (.Platform$OS.type != "unix") || quiet  # Hack to prevent console flushing problems with Windows
  if (use_fast_sv) {
    para <- 1 + (fast_sv$baseline_parameterization == "noncentered") + 2 * interweave
    parameterization <- c("centered", "noncentered", "GIS_C", "GIS_NC")[para]

    if (!quiet) {
      cat(paste("\nCalling ", parameterization, " MCMC sampler with ", draws+burnin, " iter. Series length is ", length(y), ".\n",sep=""), file=stderr())
      flush.console()
    }

    runtime <- system.time(res <-
      svsample_fast_cpp(y, draws, burnin, designmatrix, priorspec,
                        thinpara, thinlatent, keeptime,
                        startpara, startlatent, keeptau, myquiet,
                        correct_model_misspecification, interweave,
                        myoffset, fast_sv))
  } else {
    strategies <- if (interweave) c("centered", "noncentered") else expert$general_sv$starting_parameterization
    parameterization <- rep(strategies, general_sv$multi_asis)

    renameparam <- c("centered" = "C", "noncentered" = "NC")
    if (!quiet) {
      cat(paste("\nCalling ", asisprint(renameparam[parameterization], renameparam), " MCMC sampler with ", draws+burnin, " iter. Series length is ", length(y), ".\n",sep=""), file=stderr())
      flush.console()
    }

    runtime <- system.time(res <-
      svsample_general_cpp(y, draws, burnin, designmatrix, priorspec,
                           thinpara, thinlatent, keeptime,
                           startpara, startlatent, keeptau, quiet,
                           correct_model_misspecification, interweave,
                           myoffset, general_sv))
  }
  class(res) <- "svdraws"

  # Process results
  if (any(is.na(res))) {
    stop("Sampler returned NA. This is most likely due to bad input checks and shouldn't happen. Please report to package maintainer.")
  }

  if (!quiet) {
    cat("Timing (elapsed): ", file=stderr())
    cat(runtime["elapsed"], file=stderr())
    cat(" seconds.\n", file=stderr())
    cat(round((draws+burnin)/runtime[3]), "iterations per second.\n\n", file=stderr())
    cat("Converting results to coda objects... ", file=stderr())
  }

  # store results:
  res$runtime <- runtime
  res$y <- if (arorder == 0) y_orig else tail(y_orig, -arorder)
  res$simobj <- simobj
  res$para <- coda::mcmc(res$para, burnin+thinpara, burnin+draws, thinpara)
  res$latent <- coda::mcmc(res$latent, burnin+thinlatent, burnin+draws, thinlatent)
  res$latent0 <- coda::mcmc(res$latent0, burnin+thinlatent, burnin+draws, thinlatent)
  res$thinning <- list(para = thinpara, latent = thinlatent, time = keeptime)
  res$priors <- priorspec
  if (!any(is.na(designmatrix))) {
    res$beta <- coda::mcmc(res$beta, burnin+thinpara, burnin+draws, thinpara)
    res$designmatrix <- designmatrix
  } else {
    res$beta <- NULL
  }
  res$meanmodel <- meanmodel
  res$para_transform <- list(mu = function (x) {x},
                             phi = function (x) {(x+1)/2},
                             sigma = function (x) {x^2},
                             nu = function (x) {x-2},
                             rho = function (x) {(x+1)/2})
  res$para_inv_transform <- list(mu = function (x) {x},
                                 phi = function (x) {2*x-1},
                                 sigma = function (x) {sqrt(x)},
                                 nu = function (x) {x+2},
                                 rho = function (x) {2*x-1})
  res$para_transform_det <- list(mu = function (x) {1},
                                 phi = function (x) {.5},
                                 sigma = function (x) {2*x},
                                 nu = function (x) {1},
                                 rho = function (x) {.5})

  if (keeptau) {
    res$tau <- coda::mcmc(t(res$tau), burnin+thinlatent, burnin+draws, thinlatent)
  } else {
    res$tau <- NULL
  }

  if (!quiet) {
    cat("Done!\n", file=stderr())
    cat("Summarizing posterior draws... ", file=stderr())
  }
  res <- updatesummary(res, ...)

  if (!quiet) cat("Done!\n\n", file=stderr())
  res
}

#' @rdname svsample_cpp
#' @export
default_fast_sv <- 
  list(baseline_parameterization = "centered",  # "centered" or "noncentered"
       proposal_phi = "immediate acceptance-rejection",  # "immediate acceptance-rejection" or "repeated acceptance-rejection"
       proposal_sigma2 = "independence",  # "independence" or "log random walk"
       proposal_intercept_var = 1e12,  # positive number
       proposal_phi_var = 1e8,  # positive number
       proposal_sigma2_rw_scale = 0.1,  # positive number
       mh_blocking_steps = 2,  # 1/2/3
       store_indicators = FALSE,
       update = list(latent_vector = TRUE, parameters = TRUE, mixture_indicators = TRUE),
       init_indicators = 5)
#' @rdname svsample_cpp
#' @export
default_general_sv <-
  list(multi_asis = 5,  # positive integer
       starting_parameterization = "centered",  # "centered" or "noncentered"
       update = list(latent_vector = TRUE, parameters = TRUE),
       proposal_diffusion_ken = FALSE)  # FALSE turns on adaptation
