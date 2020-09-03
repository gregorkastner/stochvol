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
#' for the exponential prior distribution of the parameter
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
#' @param priorrho to be documented TODO
#' @param priorspec to be documented TODO
#' @param thinpara single number greater or equal to 1, coercible to integer.
#' Every \code{thinpara}th parameter draw is kept and returned. The default
#' value is 1, corresponding to no thinning of the parameter draws i.e. every
#' draw is stored.
#' @param thinlatent single number greater or equal to 1, coercible to integer.
#' Every \code{thinlatent}th latent variable draw is kept and returned. The
#' default value is 1, corresponding to no thinning of the latent variable
#' draws, i.e. every draw is kept.
#' @param keeptime Either 'all' (the default) or 'last'. Indicates which latent
#  volatility draws should be stored.
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
#' 
#' \code{interweaving}: Logical value. TODO
#' 
#' \code{correct_model_misspecification}: Logical value. TODO
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
#' \item{y}{the
#' argument \code{y}.}
#' \item{runtime}{\code{proc_time} object containing the
#' run time of the sampler.}
#' \item{priors}{\code{list} containing the parameter
#' values of the prior distribution, i.e. the arguments \code{priormu},
#' \code{priorphi}, \code{priorsigma}, and potentially \code{priornu} and
#' \code{priorbeta}.}
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
#' @author Gregor Kastner \email{gregor.kastner@@wu.ac.at}
#' @seealso \code{\link{svsim}}, \code{\link{updatesummary}},
#' \code{\link{predict.svdraws}}, \code{\link{plot.svdraws}}.
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
    y <- y[["y"]]
    message("Extracted data vector from 'svsim'-object.")
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
                     beta = if (meanmodel == "none") sv_constant(value = 0) else sv_normal(mean = priorbeta[1], sd = priorbeta[2]),
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
         latent0 = -10)  # init_beta
  startlatentdefault <- rep.int(-10, length(y))
  startpara <- apply_default_list(startpara, startparadefault)
  startlatent <- apply_default_list(startlatent, startlatentdefault)
  validate_initial_values(startpara, startlatent, y, designmatrix)

  # Decision about the sampler
  use_fast_sv <- TRUE
  if (!(inherits(priorspec$rho, "sv_constant") && isTRUE(priorspec$rho$value == 0))) {
    use_fast_sv <- FALSE
  }
  # TODO fully specify!!
  #  e.g. keeptime == "last" && correct_model_misspecification

  # Call sampler
  myquiet <- (.Platform$OS.type != "unix") || quiet  # Hack to prevent console flushing problems with Windows
  if (use_fast_sv) {
    para <- 1 + (fast_sv$baseline_parameterization == "noncentered") + 2 * interweave
    parameterization <- c("centered", "noncentered", "GIS_C", "GIS_NC")[para]
    mhcontrol <- if (fast_sv$proposal_sigma2 == "independence") -1 else fast_sv$proposal_sigma2_rw_scale
    gammaprior <- inherits(priorspec$sigma2, "sv_gamma")
    truncnormal <- fast_sv$proposal_phi == "repeated acceptance-rejection"
    mhsteps <- fast_sv$mh_blocking_steps
    B011 <- fast_sv$proposal_phi_var
    B022 <- fast_sv$proposal_intercept_var

    if (!quiet) {
      cat(paste("\nCalling ", parameterization, " MCMC sampler with ", draws+burnin, " iter. Series length is ", length(y), ".\n",sep=""), file=stderr())
      flush.console()
    }

    runtime <- system.time(res <-
      svsample_cpp(y, draws, burnin, designmatrix, priorspec,
                   thinpara, thinlatent, keeptime,
                   startpara, startlatent, keeptau, myquiet,
                   correct_model_misspecification, interweave,
                   myoffset, fast_sv))

    res$latent <- t(res$latent)
    class(res) <- "svdraws"
  } else {
    strategies <- if (interweave) c("centered", "noncentered") else "centered"
    parameterization <- rep(strategies, general_sv$multi_asis)
    use.mala <- general_sv$proposal_para == "metropolis-adjusted langevin algorithm"
    gammaprior <- inherits(priorspec$sigma2, "sv_gamma")
    correct.latent.draws <- correct_model_misspecification

    renameparam <- c("centered" = "C", "noncentered" = "NC")
    if (!quiet) {
      cat(paste("\nCalling ", asisprint(renameparam[parameterization], renameparam), " MCMC sampler with ", draws+burnin, " iter. Series length is ", length(y), ".\n",sep=""), file=stderr())
      flush.console()
    }

    runtime <- system.time(res <-
      svlsample_cpp(y, draws, burnin, designmatrix, priorspec,
                    thinpara, thinlatent, keeptime,
                    startpara, startlatent, keeptau, quiet,
                    correct_model_misspecification, interweave,
                    myoffset, general_sv))

    class(res) <- c("svldraws", "svdraws")
  }

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
  if (keeptime == "all") {
    colnames(res$latent) <- paste0('h_', arorder + seq_along(y))
  } else if (keeptime == "last") {
    colnames(res$latent) <- paste0('h_', arorder + length(y))
  }
  res$runtime <- runtime
  res$y <- y_orig
  res$para <- coda::mcmc(res$para, burnin+thinpara, burnin+draws, thinpara)
  res$latent <- coda::mcmc(res$latent, burnin+thinlatent, burnin+draws, thinlatent)
  res$latent0 <- coda::mcmc(res$latent0, burnin+thinlatent, burnin+draws, thinlatent)
  res$thinning <- list(para = thinpara, latent = thinlatent, time = keeptime)
  res$priors <- c(list(mu = priormu,
                       phi = priorphi,
                       sigma = priorsigma,
                       gammaprior = gammaprior),
                  if (priornu <= 0) list() else list(nu = priornu),
                  if (isTRUE(inherits(priorspec$rho, "sv_beta"))) list(rho = c(priorspec$rho$alpha, priorspec$rho$beta)) else list())
  if (!any(is.na(designmatrix))) {
    res$beta <- coda::mcmc(res$beta, burnin+thinpara, burnin+draws, thinpara)
    colnames(res$beta) <- paste0("b_", 0:(NCOL(designmatrix)-1))
    res$priors <- c(res$priors, "beta" = list(priorbeta), "designmatrix" = list(designmatrix))
  } else {
    res$beta <- NULL
  }
  res$meanmodel <- meanmodel

  if (keeptau) {
    res$tau <- coda::mcmc(t(res$tau), burnin+thinlatent, burnin+draws, thinlatent)
    if (keeptime == "all")
      colnames(res$tau) <- paste0('tau_', arorder + seq_along(y))
    else if (keeptime == "last")
      colnames(res$tau) <- paste0('tau_', arorder + length(y))
  }

  if (!quiet) {
    cat("Done!\n", file=stderr())
    cat("Summarizing posterior draws... ", file=stderr())
  }
  res <- updatesummary(res, ...)

  if (!quiet) cat("Done!\n\n", file=stderr())
  res
}

#' @rdname svsample
#' @export
svsample_fast <- function(y, draws = 1, burnin = 0, designmatrix = matrix(NA, 1, 1), priorspec,
                          thinpara = 1, thinlatent = 1, keeptime = "all", keeptau = FALSE,
                          quiet = TRUE, startpara, startlatent) {
  svsample_cpp(y, draws, burnin, designmatrix, priorspec,
               thinpara, thinlatent, keeptime,
               startpara, startlatent, keeptau, quiet,
               FALSE, TRUE,
               0, default_fast_sv)
}

#' @rdname svsample
#' @export
svsample_general <- function(y, draws = 1, burnin = 0, designmatrix = matrix(NA, 1, 1), priorspec,
                             thinpara = 1, thinlatent = 1, keeptime = "all", keeptau = FALSE,
                             quiet = TRUE, startpara, startlatent) {
  svlsample_cpp(y, draws, burnin, designmatrix, priorspec,
                thinpara, thinlatent, keeptime,
                startpara, startlatent, keeptau, quiet,
                FALSE, TRUE,
                0, default_general_sv)
}

default_fast_sv <- 
  list(baseline_parameterization = "centered",  # "centered" or "noncentered"
       proposal_phi = "immediate acceptance-rejection",  # "immediate acceptance-rejection" or "repeated acceptance-rejection"
       proposal_sigma2 = "independence",  # "independence" or "log random walk"
       proposal_intercept_var = 1e12,  # positive number
       proposal_phi_var = 1e8,  # positive number
       proposal_sigma2_rw_scale = 0.1,  # positive number
       mh_blocking_steps = 2,  # 1/2/3
       store_indicators = FALSE,
       init_indicators = 5)

default_general_sv <-
  list(multi_asis = 3,  # positive integer
       starting_parameterization = "centered",  # "centered" or "noncentered"
       proposal_para = "random walk",  # "random walk" or "metropolis-adjusted langevin algorithm
       proposal_diffusion_ken = NULL)

##' @rdname svsample
##' @export
#svtsample <- function(y, draws = 10000, burnin = 1000, designmatrix = NA,
#                      priormu = c(0, 100), priorphi = c(5, 1.5), priorsigma = 1,
#                      priornu = 0.1, priorrho = NA,
#                      priorbeta = c(0, 10000), priorlatent0 = "stationary",
#                      priorspec = NULL,
#                      thinpara = 1, thinlatent = 1, keeptime = "all", keeptau = FALSE,
#                      quiet = FALSE, startpara = NULL, startlatent = NULL, expert = NULL, ...) {
#  svsample(y, draws = draws, burnin = burnin, designmatrix = designmatrix,
#           priormu = priormu, priorphi = priorphi, priorsigma = priorsigma,
#           priornu = priornu, priorrho = priorrho,
#           priorbeta = priorbeta, priorlatent0 = priorlatent0,
#           priorspec = priorspec,
#           thinpara = thinpara, thinlatent = thinlatent, keeptime = keeptime,
#           keeptau = keeptau, quiet = quiet, startpara = startpara, startlatent = startlatent,
#           expert = expert, ...)
#}
#
#
##' Minimal overhead version of \code{\link{svsample}}.
##' 
##' \code{svsample2} is a minimal overhead version of \code{\link{svsample}}
##' with slightly different default arguments and a simplified return value
##' structure. It is intended to be used mainly for one-step updates where speed
##' is an issue, e.g., as a plug-in into other MCMC samplers. Note that
##' absolutely no input checking is performed, thus this function is to be used
##' with proper care!
##' 
##' As opposed to the ordinary \code{\link{svsample}}, the default values differ
##' for \code{draws}, \code{burnin}, and \code{quiet}. Note that currently
##' neither \code{expert} nor \code{\dots{}} arguments are provided.
##' 
##' @param y numeric vector containing the data (usually log-returns), which
##' must not contain zeroes.
##' @param draws single number greater or equal to 1, indicating the number of
##' draws after burn-in (see below). Will be automatically coerced to integer.
##' The defaults value is 1.
##' @param burnin single number greater or equal to 0, indicating the number of
##' draws discarded as burn-in. Will be automatically coerced to integer. The
##' default value is 0.
##' @param priormu numeric vector of length 2, indicating mean and standard
##' deviation for the Gaussian prior distribution of the parameter \code{mu},
##' the level of the log-volatility. The default value is \code{c(0, 100)},
##' which constitutes a practically uninformative prior for common exchange rate
##' datasets, stock returns and the like.
##' @param priorphi numeric vector of length 2, indicating the shape parameters
##' for the Beta prior distribution of the transformed parameter
##' \code{(phi + 1) / 2}, where \code{phi} denotes the persistence of the
##' log-volatility. The default value is \code{c(5, 1.5)}, which constitutes a
##' prior that puts some belief in a persistent log-volatility but also
##' encompasses the region where \code{phi} is around 0.
##' @param priorsigma single positive real number, which stands for the scaling
##' of the transformed parameter \code{sigma^2}, where \code{sigma} denotes the
##' volatility of log-volatility. More precisely, \code{sigma^2 ~ priorsigma *
##' chisq(df = 1)}. The default value is \code{1}, which constitutes a
##' reasonably vague prior for many common exchange rate datasets, stock returns
##' and the like.
##' @param priornu single non-negative number, indicating the rate parameter
##' for the exponential prior distribution of the parameter
##' \code{nu}, the degrees-of-freedom parameter of the conditional innovations
##' t-distribution. The default value is \code{0}, fixing the
##' degrees-of-freedom to infinity. This corresponds to conditional standard
##' normal innovations, the pre-1.1.0 behavior of \pkg{stochvol}.
##' @param priorlatent0 either a single non-negative number or the string
##' \code{'stationary'} (the default, also the behavior before version 1.3.0).
##' When \code{priorlatent0} is equal to \code{'stationary'}, the stationary
##' distribution of the latent AR(1)-process is used as the prior for the
##' initial log-volatility \code{h_0}. When \code{priorlatent0} is equal to a
##' number \eqn{B}, we have \eqn{h_0 \sim N(\mu, B\sigma^2)} a priori.
##' @param priorrho to be documented TODO
##' @param priorspec to be documented TODO
##' @param thinpara single number greater or equal to 1, coercible to integer.
##' Every \code{thinpara}th parameter draw is kept and returned. The default
##' value is 1, corresponding to no thinning of the parameter draws -- every
##' draw is stored.
##' @param thinlatent single number greater or equal to 1, coercible to integer.
##' Every \code{thinlatent}th latent variable draw is kept and returned. The
##' default value is 1, corresponding to no thinning of the latent variable
##' draws, i.e. every draw is kept.
##' @param keeptime Either 'all' (the default) or 'last'. Indicates which latent
##  volatility draws should be stored.
##' @param keeptau logical value indicating whether the 'variance inflation
##' factors' should be stored (used for the sampler with conditional t
##' innovations only). This may be useful to check at what point(s) in time the
##' normal disturbance had to be 'upscaled' by a mixture factor and when the
##' series behaved 'normally'.
##' @param quiet logical value indicating whether the progress bar and other
##' informative output during sampling should be omitted. The default value is
##' \code{TRUE}, implying non-verbose output.
##' @param startpara \emph{compulsory} named list, containing the starting
##' values for the parameter draws. \code{startpara} must contain three elements
##' named \code{mu}, \code{phi}, and \code{sigma}, where \code{mu} is an
##' arbitrary numerical value, \code{phi} is a real number between \code{-1} and
##' \code{1}, and \code{sigma} is a positive real number. Moreover, if
##' \code{priornu} is not \code{0}, \code{startpara} must also contain an
##' element named \code{nu} (the degrees of freedom parameter for the
##' t-innovations).
##' @param startlatent \emph{compulsory} vector of length \code{length(x$y)},
##' containing the starting values for the latent log-volatility draws.
##' @return A list with three components:
##' \item{para}{\code{3} times
##' \code{draws} matrix containing the parameter draws. If \code{priornu} is not
##' \code{0}, this is a \code{4} times \code{draws} matrix.}
##' \item{latent}{\code{length(y)} times \code{draws} matrix containing draws of
##' the latent variables \code{h_1, \dots{}, h_n}.}
##' \item{latent0}{Vector of
##' length \code{draws} containing the draw(s) of the initial latent variable
##' \code{h_0}.}
##' @note Please refer to the package vignette for an example.
##' @section Warning: Expert use only! For most applications, the use of
##' \code{\link{svsample}} is recommended.
##' @author Gregor Kastner \email{gregor.kastner@@wu.ac.at}
##' @seealso \code{\link{svsample}}
##' @keywords models ts
##' @examples
##' data(exrates)
##' aud.price <- subset(exrates,
##'   as.Date("2010-01-01") <= date & date < as.Date("2011-01-01"),
##'   "AUD")[,1]
##' draws <- svsample2(logret(aud.price),
##'                    draws = 10, burnin = 0,
##'                    startpara = list(phi = 0.95, mu = -10, sigma = 0.2),
##'                    startlatent = rep_len(-10, length(aud.price) - 1))
##' @export
#svsample2 <- function(y, draws = 1, burnin = 1, designmatrix = NA,
#                       priormu = c(0, 100), priorphi = c(5, 1.5), priorsigma = 1,
#                       priornu = 0, priorrho = c(4, 4),
#                       priorbeta = c(0, 10000), priorlatent0 = "stationary",
#                       priorspec = NULL,
#                       thinpara = 1, thinlatent = 1, keeptime = "all", keeptau = FALSE,
#                       quiet = FALSE, startpara = NULL, startlatent = NULL, expert = NULL, ...) {
#  # TODO make minimal overhead
#  svsample(y, draws = draws, burnin = burnin, designmatrix = designmatrix,
#           priormu = priormu, priorphi = priorphi, priorsigma = priorsigma,
#           priornu = priornu, priorrho = priorrho,
#           priorbeta = priorbeta, priorlatent0 = priorlatent0,
#           priorspec = priorspec,
#           thinpara = thinpara, thinlatent = thinlatent, keeptime = keeptime,
#           keeptau = keeptau, quiet = quiet, startpara = startpara, startlatent = startlatent,
#           expert = expert, ...)
#}
#
#svsample2.old <- function(y, draws = 1, burnin = 0, priormu = c(0, 100),
#                      priorphi = c(5, 1.5), priorsigma = 1, priornu = 0,
#                      priorlatent0 = "stationary", thinpara = 1, thinlatent = 1,
#                      keeptime = "all", keeptau = FALSE,
#                      quiet = TRUE, startpara, startlatent) {
#
#  if (priorlatent0 == "stationary")
#    priorlatent0 <- -1L
#
#  if (keeptime == "all")
#    thintime <- 1L
#  else if (keeptime == "last")
#    thintime <- length(y)
#  else
#    stop("Parameter 'keeptime' must be either 'all' or 'last'.")
#
#  res <- svsample_cpp(y, draws, burnin, matrix(NA), priormu[1], priormu[2]^2,
#                      priorphi[1], priorphi[2], priorsigma, thinpara, thinlatent,
#                      thintime, startpara, startlatent, keeptau, quiet, 3L, 2L, 10^8,
#                      10^12, -1, TRUE, FALSE, 0, FALSE, priornu, c(NA, NA), priorlatent0)
#
#  res$para <- t(res$para)
#  if (NROW(res$para) == 3) {
#    rownames(res$para) <- names(res$para) <- c("mu", "phi", "sigma")
#  } else {
#    rownames(res$para) <- names(res$para) <- c("mu", "phi", "sigma", "nu")
#  }
#  res.meanmodel <- "none"
#
#  res
#}
#
#
##' Markov Chain Monte Carlo (MCMC) Sampling for the Stochastic Volatility
##' Model with Leverage (SVL)
##' 
##' \code{svlsample} simulates from the joint posterior distribution of the SVL
##' parameters \code{mu}, \code{phi}, \code{sigma}, and \code{rho},
##' along with the latent log-volatilities \code{h_1,...,h_n} and returns the
##' MCMC draws. If a design matrix is provided, simple Bayesian regression can
##' also be conducted.
##' 
##' @param y numeric vector containing the data (usually log-returns), which
##' must not contain zeros. Alternatively, \code{y} can be an \code{svsim}
##' object. In this case, the returns will be extracted and a warning is thrown.
##' @param draws single number greater or equal to 1, indicating the number of
##' draws after burn-in (see below). Will be automatically coerced to integer.
##' The default value is 10000.
##' @param burnin single number greater or equal to 0, indicating the number of
##' draws discarded as burn-in. Will be automatically coerced to integer. The
##' default value is 1000.
##' @param designmatrix regression design matrix for modeling the mean. Must
##' have \code{length(y)} rows. Alternatively, \code{designmatrix} may be a
##' string of the form \code{"arX"}, where \code{X} is a nonnegative integer. To
##' fit a constant mean model, use \code{designmatrix = "ar0"} (which is
##' equivalent to \code{designmatrix = matrix(1, nrow = length(y))}). To fit an
##' AR(1) model, use \code{designmatrix = "ar1"}, and so on. If some elements of
##' \code{designmatrix} are \code{NA}, the mean is fixed to zero (pre-1.2.0
##' behavior of \pkg{stochvol}).
##' @param priormu numeric vector of length 2, indicating mean and standard
##' deviation for the Gaussian prior distribution of the parameter \code{mu},
##' the level of the log-volatility. The default value is \code{c(0, 100)},
##' which constitutes a practically uninformative prior for common exchange rate
##' datasets, stock returns and the like.
##' @param priorphi numeric vector of length 2, indicating the shape parameters
##' for the Beta prior distribution of the transformed parameter
##' \code{(phi + 1) / 2}, where \code{phi} denotes the persistence of the
##' log-volatility. The default value is \code{c(5, 1.5)}, which constitutes a
##' prior that puts some belief in a persistent log-volatility but also
##' encompasses the region where \code{phi} is around 0.
##' @param priorsigma single positive real number, which stands for the scaling
##' of the transformed parameter \code{sigma^2}, where \code{sigma} denotes the
##' volatility of log-volatility. More precisely, \code{sigma^2 ~ priorsigma *
##' chisq(df = 1)}. The default value is \code{1}, which constitutes a
##' reasonably vague prior for many common exchange rate datasets, stock returns
##' and the like.
##' @param priorrho numeric vector of length 2, indicating the shape parameters
##' for the Beta prior distribution of the transformed parameter
##' \code{(rho + 1) / 2}, where \code{rho} denotes the conditional correlation
##' between observation and the increment of the
##' log-volatility. The default value is \code{c(4, 4)}, which constitutes a
##' slightly informative prior around 0 (the no leverage case) to boost convergence.
##' @param priorbeta numeric vector of length 2, indicating the mean and
##' standard deviation of the Gaussian prior for the regression parameters. The
##' default value is \code{c(0, 10000)}, which constitutes a very vague prior
##' for many common datasets. Not used if \code{designmatrix} is \code{NA}.
##' @param priorspec to be documented TODO
##' @param thinpara single number greater or equal to 1, coercible to integer.
##' Every \code{thinpara}th parameter draw is kept and returned. The default
##' value is 1, corresponding to no thinning of the parameter draws i.e. every
##' draw is stored.
##' @param thinlatent single number greater or equal to 1, coercible to integer.
##' Every \code{thinlatent}th latent variable draw is kept and returned. The
##' default value is 1, corresponding to no thinning of the latent variable
##' draws, i.e. every draw is kept.
##' @param keeptime Either 'all' (the default) or 'last'. Indicates which latent
##  volatility draws should be stored.
##' @param quiet logical value indicating whether the progress bar and other
##' informative output during sampling should be omitted. The default value is
##' \code{FALSE}, implying verbose output.
##' @param startpara \emph{optional} named list, containing the starting values
##' for the parameter draws. If supplied, \code{startpara} must contain four
##' elements named \code{mu}, \code{phi}, \code{sigma}, and \code{rho}, where \code{mu} is
##' an arbitrary numerical value, \code{phi} is a real number between \code{-1}
##' and \code{1}, \code{sigma} is a positive real number, and \code{rho} is
##' a real number between \code{-1} and \code{1}. The default value is equal
##' to the prior mean.
##' @param startlatent \emph{optional} vector of length \code{length(y)},
##' containing the starting values for the latent log-volatility draws. The
##' default value is \code{rep(-10, length(y))}.
##' @param expert \emph{optional} named list of expert parameters. For most
##' applications, the default values probably work best. If
##' \code{expert} is provided, it may contain the following named elements:
##' 
##' \code{parameterization}: Character string containing values \code{"centered"},
##' and \code{"noncentered"}. Alternatively, it can be a single element character
##' vector of the form \code{"asisX"}, where \code{X} is an integer, which is
##' equivalent to \code{rep(c("centered", "noncentered"), X)}.
##' Defaults to \code{"asis5"}.
##' 
##' \code{gammaprior}: Single logical value indicating whether a Gamma prior for
##' \code{sigma^2} should be used. If set to \code{FALSE}, a moment-matched Inverse Gamma
##' prior is employed. Defaults to \code{TRUE}.
##' 
##' \code{init.with.svsample}: Single integer indicating the length of a ``pre-burnin'' run using
##' the computationally much more efficient \code{\link{svsample}}. This run helps
##' in finding good initial values for the latent states, giving \code{svlsample}
##' a considerable initial boost for convergence. Defaults to \code{1000L}.
##' 
##' \code{mhcontrol}: Either a single numeric value specifying the diagonal elements of
##' a diagonal covariance matrix, or a list with two elements, both single numeric values
##' (explained later), or a 4x4 covariance matrix. Argument \code{mhcontrol} controls the
##' proposal density of a Metropolis-Hastings (MH) update step when jointly sampling \code{mu},
##' \code{phi}, \code{sigma}, and \code{rho}. It specifies the covariance matrix of a
##' log-random-walk proposal. In case \code{mhcontrol} is a list of length two, its elements
##' have to be \code{scale} and \code{rho.var}. In this case, the covariance matrix is calculated
##' from the pre-burnin step with \code{\link{svsample}}, which gives an approximate
##' posterior structure of the second moment for \code{mu}, \code{phi}, and \code{sigma}.
##' This covariance matrix is then extended with \code{mhcontrol$rho.var}, specifying the
##' variance for \code{rho}. The off-diagonal elements belonging to \code{rho} are set
##' to 0. Finally, the whole covariance matrix is scaled by \code{mhcontrol$scale}. For
##' this case to work, \code{init.with.svsample} has to be positive.
##' Defaults to \code{list(scale=0.35, rho.var=0.02)}.
##' 
##' \code{correct.latent.draws}: Single logical value indicating whether to correct
##' the draws obtained from the auxiliary model of Omori, et al. (2007). Defaults
##' to \code{FALSE}.
##' @param \dots Any extra arguments will be forwarded to
##' \code{\link{updatesummary}}, controlling the type of statistics calculated
##' for the posterior draws.
##' @return The value returned is a list object of class \code{svldraws} holding
##' \item{para}{\code{mcmc} object containing the \emph{parameter} draws from
##' the posterior distribution.}
##' \item{latent}{\code{mcmc} object containing the
##' \emph{latent instantaneous log-volatility} draws from the posterior
##' distribution.}
##' \item{latent0}{\code{mcmc} object containing the \emph{latent
##' initial log-volatility} draws from the posterior distribution.}
##' \item{beta}{\code{mcmc} object containing the \emph{regression coefficient}
##' draws from the posterior distribution \emph{(optional)}.}
##' \item{y}{the argument \code{y}.}
##' \item{runtime}{\code{proc_time} object containing the
##' run time of the sampler.}
##' \item{priors}{\code{list} containing the parameter
##' values of the prior distribution, i.e. the arguments \code{priormu},
##' \code{priorphi}, \code{priorsigma}, and \code{priorrho}, and potentially
##' \code{priorbeta}.}
##' \item{thinning}{\code{list} containing the thinning
##' parameters, i.e. the arguments \code{thinpara}, \code{thinlatent} and
##' \code{keeptime}.}
##' \item{summary}{\code{list} containing a collection of
##' summary statistics of the posterior draws for \code{para}, \code{latent}, and \code{latent0}.}
##' \item{meanmodel}{\code{character} containing information about how \code{designmatrix}
##' was used.}
##' 
##' To display the output, use \code{print}, \code{summary} and \code{plot}. The
##' \code{print} method simply prints the posterior draws (which is very likely
##' a lot of output); the \code{summary} method displays the summary statistics
##' currently stored in the object; the \code{plot} method
##' \code{\link{plot.svdraws}} gives a graphical overview of the posterior
##' distribution by calling \code{\link{volplot}}, \code{\link{traceplot}} and
##' \code{\link{densplot}} and displaying the results on a single page.
##' @note If \code{y} contains zeros, you might want to consider de-meaning your
##' returns or use \code{designmatrix = "ar0"}. We use the Metropolis-Hastings
##' algorithm for sampling the latent vector \code{h}, where the proposal is a
##' draw from an auxiliary mixture approximation model [Omori, et al. (2007)].
##' We draw the parameters \code{mu}, \code{phi}, \code{sigma}, and \code{rho}
##' jointly by employing a Metropolis random walk step. By default, we boost the
##' random walk through the repeated application of the ancillarity-sufficiency
##' interweaving strategy (ASIS) [Yu, Meng (2011)]. A message in the beginning
##' of sampling indicates the interweaving strategy used, which can be modified
##' through parameter \code{expert}.
##' @author Darjus Hosszejni \email{darjus.hosszejni@@wu.ac.at}
##' @references
##' Yu, Y. and Meng, X.-L. (2011).
##' To Center or not to Center: That is not the Question---An Ancillarity-Sufficiency
##' Interweaving Strategy (ASIS) for Boosting MCMC Efficiency. \emph{Journal of
##' Computational and Graphical Statistics}, \bold{20}(3), 531--570,
##' \url{http://dx.doi.org/10.1198/jcgs.2011.203main}
##'
##' Omori, Y. and Chib, S. and Shephard, N. and Nakajima, J. (2007).
##' Stochastic Volatility with Leverage: Fast and Efficient Likelihood Inference.
##' \emph{Journal of Econometrics}, \bold{140}(2), 425--449,
##' \url{http://dx.doi.org/10.1016/j.jeconom.2006.07.008}
##' @seealso \code{\link{svsim}}, \code{\link{svsample}}, \code{\link{updatesummary}},
##' \code{\link{predict.svdraws}}, \code{\link{plot.svdraws}}.
##' @keywords models ts
##' @examples
##' \dontrun{
##' # Example 1
##' ## Simulate a short SVL process
##' sim <- svsim(200, mu = -10, phi = 0.95, sigma = 0.2, rho = -0.4)
##' 
##' ## Obtain 5000 draws from the sampler (that's not a lot)
##' draws <- svlsample(sim$y)
##' 
##' ## Check out the results
##' summary(draws)
##' plot(draws, simobj = sim)
##' 
##' 
##' # Example 2
##' ## AR(1) structure for the mean
##' data(exrates)
##' len <- 1200
##' ahead <- 100
##' y <- head(exrates$USD, len)
##' 
##' ## Fit AR(1)-SVL model to EUR-USD exchange rates
##' res <- svlsample(y, designmatrix = "ar1")
##' 
##' ## Use predict.svdraws to obtain predictive distributions
##' preddraws <- predict(res, steps = ahead)
##' 
##' ## Calculate predictive quantiles
##' predquants <- apply(preddraws$y, 2, quantile, c(.1, .5, .9))
##' 
##' ## Visualize
##' expost <- tail(head(exrates$USD, len+ahead), ahead)
##' ts.plot(y, xlim = c(length(y)-4*ahead, length(y)+ahead),
##' 	       ylim = range(c(predquants, expost, tail(y, 4*ahead))))
##' for (i in 1:3) {
##'   lines((length(y)+1):(length(y)+ahead), predquants[i,],
##'         col = 3, lty = c(2, 1, 2)[i])
##' }
##' lines((length(y)+1):(length(y)+ahead), expost,
##'       col = 2)
##' 
##' 
##' # Example 3
##' ## Predicting USD based on JPY and GBP in the mean
##' data(exrates)
##' len <- 1200
##' ahead <- 30
##' ## Calculate log-returns
##' logreturns <- apply(exrates[, c("USD", "JPY", "GBP")], 2,
##'                     function (x) diff(log(x)))
##' logretUSD <- logreturns[2:(len+1), "USD"]
##' regressors <- cbind(1, as.matrix(logreturns[1:len, ]))  # lagged by 1 day
##' 
##' ## Fit SV model to EUR-USD exchange rates
##' res <- svlsample(logretUSD, designmatrix = regressors)
##' 
##' ## Use predict.svdraws to obtain predictive distributions
##' predregressors <- cbind(1, as.matrix(logreturns[(len+1):(len+ahead), ]))
##' preddraws <- predict(res, steps = ahead,
##'                      newdata = predregressors)
##' predprice <- exrates[len+2, "USD"] * exp(t(apply(preddraws$y, 1, cumsum)))
##' 
##' ## Calculate predictive quantiles
##' predquants <- apply(predprice, 2, quantile, c(.1, .5, .9))
##' 
##' ## Visualize
##' priceUSD <- exrates[3:(len+2), "USD"]
##' expost <- exrates[(len+3):(len+ahead+2), "USD"]
##' ts.plot(priceUSD, xlim = c(len-4*ahead, len+ahead+1),
##' 	       ylim = range(c(expost, predquants, tail(priceUSD, 4*ahead))))
##' for (i in 1:3) {
##'   lines(len:(len+ahead), c(tail(priceUSD, 1), predquants[i,]),
##'         col = 3, lty = c(2, 1, 2)[i])
##' }
##' lines(len:(len+ahead), c(tail(priceUSD, 1), expost),
##'       col = 2)
##' }
##' @export
#svlsample <- function(y, draws = 10000, burnin = 1000, designmatrix = NA,
#                      priormu = c(0, 100), priorphi = c(5, 1.5), priorsigma = 1,
#                      priornu = 0, priorrho = c(4, 4),
#                      priorbeta = c(0, 10000), priorlatent0 = "stationary",
#                      priorspec = NULL,
#                      thinpara = 1, thinlatent = 1, keeptime = "all", keeptau = FALSE,
#                      quiet = FALSE, startpara = NULL, startlatent = NULL, expert = NULL, ...) {
#  svsample(y, draws = draws, burnin = burnin, designmatrix = designmatrix,
#           priormu = priormu, priorphi = priorphi, priorsigma = priorsigma,
#           priornu = priornu, priorrho = priorrho,
#           priorbeta = priorbeta, priorlatent0 = priorlatent0,
#           priorspec = priorspec,
#           thinpara = thinpara, thinlatent = thinlatent, keeptime = keeptime,
#           keeptau = keeptau, quiet = quiet, startpara = startpara, startlatent = startlatent,
#           expert = expert, ...)
#}
#
#svlsample.old <- function (y, draws = 10000, burnin = 1000, designmatrix = NA,
#                       priormu = c(0, 100), priorphi = c(5, 1.5), priorsigma = 1,
#                       priorrho = c(4, 4), priorbeta = c(0, 10000),
#                       thinpara = 1, thinlatent = 1, keeptime = "all",
#                       quiet = FALSE, startpara, startlatent, expert, ...) {
#  # Some error checking for y
#  if (inherits(y, "svsim")) {
#    y <- y[["y"]]
#    warning("Extracted data vector from 'svsim'-object.")
#  }
#  if (!is.numeric(y)) stop("Argument 'y' (data vector) must be numeric.")
#  if (length(y) < 2) stop("Argument 'y' (data vector) must contain at least two elements.")
#
#  y_orig <- y
#  y <- as.vector(y)
#  
#  myoffset <- if (any(is.na(designmatrix)) && any(y^2 == 0)) 1.0e-8 else 0
#  if (myoffset > 0) {
#    warning(paste("Argument 'y' (data vector) contains zeros. I am adding an offset constant of size ", myoffset, " to do the auxiliary mixture sampling. If you want to avoid this, you might consider de-meaning the returns before calling this function.", sep=""))
#  }
#
#  # Some error checking for draws
#  if (!is.numeric(draws) || length(draws) != 1 || draws < 1) {
#    stop("Argument 'draws' (number of MCMC iterations after burn-in) must be a single number >= 1.")
#  } else {
#    draws <- as.integer(draws)
#  }
#
#  # Some error checking for burnin
#  if (!is.numeric(burnin) || length(burnin) != 1 || burnin < 0) {
#    stop("Argument 'burnin' (burn-in period) must be a single number >= 0.")
#  } else {
#    burnin <- as.integer(burnin)
#  }
#
#  # Some error checking for designmatrix
#  meanmodel <- "matrix"
#  arorder <- 0L
#  if (any(is.na(designmatrix))) {
#    designmatrix <- matrix(NA)
#    meanmodel <- "none"
#  } else {
#    if (any(grep("ar[0-9]+$", as.character(designmatrix)[1]))) {
#      arorder <- as.integer(gsub("ar", "", as.character(designmatrix)))
#      if (length(y) <= (arorder + 1L)) stop("Time series 'y' is too short for this AR process.")
#      designmatrix <- matrix(1, nrow = length(y) - arorder, ncol = 1)
#      colnames(designmatrix) <- c("const")
#      meanmodel <- "constant"
#      if (arorder >= 1) {
#        for (i in seq_len(arorder)) {
#          oldnames <- colnames(designmatrix)
#          designmatrix <- cbind(designmatrix, y[(arorder-i+1):(length(y)-i)])
#          colnames(designmatrix) <- c(oldnames, paste0("ar", i))
#        }
#        y <- y[-(1:arorder)]
#        meanmodel <- paste0("ar", arorder)
#      }
#    }
#    if (!is.numeric(designmatrix)) stop("Argument 'designmatrix' must be a numeric matrix or an AR-specification.")
#    if (!is.matrix(designmatrix)) {
#      designmatrix <- matrix(designmatrix, ncol = 1)
#    }
#    if (nrow(designmatrix) != length(y)) stop("Number of columns of argument 'designmatrix' must be equal to length(y).")
#  }
#
#  # Some error checking for the prior parameters 
#  if (!is.numeric(priormu) || length(priormu) != 2) {
#    stop("Argument 'priormu' (mean and variance for the Gaussian prior for mu) must be numeric and of length 2.")
#  }
#
#  if (!is.numeric(priorphi) || length(priorphi) != 2) {
#    stop("Argument 'priorphi' (shape1 and shape2 parameters for the Beta prior for (phi + 1) / 2) must be numeric and of length 2.")
#  }
#
#  if (!is.numeric(priorsigma) || length(priorsigma) != 1 || priorsigma <= 0) {
#    stop("Argument 'priorsigma' (scaling of the chi-squared(df = 1) prior for sigma^2) must be a single number > 0.")
#  }
#
#  if (!is.numeric(priorrho) || length(priorrho) != 2) {
#    stop("Argument 'priorrho' (shape1 and shape2 parameters for the Beta prior for (rho + 1) / 2) must be numeric and of length 2.")
#  }
#
#  if (!is.numeric(priorbeta) || length(priorbeta) != 2) {
#    stop("Argument 'priorbeta' (means and sds for the independent Gaussian priors for beta) must be numeric and of length 2.")
#  }
#
#  # Some error checking for thinpara
#  if (!is.numeric(thinpara) || length(thinpara) != 1 || thinpara < 1) {
#    stop("Argument 'thinpara' (thinning parameter for mu, phi, sigma, and rho) must be a single number >= 1.")
#  } else {
#    thinpara <- as.integer(thinpara)
#  }
#
#  # Some error checking for thinlatent
#  if (!is.numeric(thinlatent) || length(thinlatent) != 1 || thinlatent < 1) {
#    stop("Argument 'thinlatent' (thinning parameter for the latent log-volatilities) must be a single number >= 1.")
#  } else {
#    thinlatent <- as.integer(thinlatent)
#  }
#
#  # Some error checking for keeptime
#  if (length(keeptime) != 1L || !is.character(keeptime) || !(keeptime %in% c("all", "last"))) {
#    stop("Parameter 'keeptime' must be either 'all' or 'last'.")
#  } else {
#    if (keeptime == "all") thintime <- 1L else if (keeptime == "last") thintime <- length(y)
#  }
#
#  # Some input checking for startpara
#  startparadefault <- list(mu = priormu[1],
#                           phi = 2 * (priorphi[1] / sum(priorphi)) - 1,
#                           sigma = priorsigma,
#                           rho = 2 * (priorrho[1] / sum(priorrho)) - 1)
#  if (missing(startpara)) {
#    startpara <- startparadefault
#  } else {
#    if (!is.list(startpara))
#      stop("Argument 'startpara' must be a list. Its elements must be named 'mu', 'phi', 'sigma', and 'rho'.")
#
#    if (!is.numeric(startpara[["mu"]]))
#      stop('Argument \'startpara[["mu"]]\' must exist and be numeric.')
#
#    if (!is.numeric(startpara[["phi"]]))
#      stop('Argument \'startpara[["phi"]]\' must exist and be numeric.')
#
#    if (abs(startpara[["phi"]]) >= 1)
#      stop('Argument \'startpara[["phi"]]\' must be between -1 and 1.')
#
#    if (!is.numeric(startpara[["sigma"]]))
#      stop('Argument \'startpara[["sigma"]]\' must exist and be numeric.')
#
#    if (startpara[["sigma"]] <= 0)
#      stop('Argument \'startpara[["sigma"]]\' must be positive.')
#
#    if (!is.numeric(startpara[["rho"]]))
#      stop('Argument \'startpara[["rho"]]\' must exist and be numeric.')
#
#    if (abs(startpara[["rho"]]) >= 1)
#      stop('Argument \'startpara[["rho"]]\' must be between -1 and 1.')
#  }
#
#  # Some input checking for startlatent
#  if (missing(startlatent)) {
#    startlatent <- rep(-10, length(y))
#  } else {
#    if (!is.numeric(startlatent) | length(startlatent) != length(y))
#      stop("Argument 'startlatent' must be numeric and of the same length as the data 'y'.")
#  }
#
#  # Some error checking for expert
#  strategies <- c("centered", "noncentered")
#  expertdefault <- list(parameterization = rep(strategies, 5),  # default: ASISx5
#                        mhcontrol = list(use.mala = FALSE),
#                        gammaprior = TRUE,
#                        correct.latent.draws = FALSE)
#  if (missing(expert)) {
#    parameterization <- expertdefault$parameterization
#    mhcontrol <- expertdefault$mhcontrol
#    use.mala <- mhcontrol$use.mala
#    gammaprior <- expertdefault$gammaprior
#    correct.latent.draws <- expertdefault$correct.latent.draws
#  } else {
#    expertnames <- names(expert)
#    if (!is.list(expert) || is.null(expertnames) || any(expertnames == ""))
#      stop("Argument 'expert' must be a named list with nonempty names.")
#    if (length(unique(expertnames)) != length(expertnames))
#      stop("No duplicate elements allowed in argument 'expert'.")
#    allowednames <- c("parameterization", "mhcontrol", "gammaprior", "correct.latent.draws")
#    exist <- pmatch(expertnames, allowednames)
#    if (any(is.na(exist)))
#      stop(paste("Illegal element '", paste(expertnames[is.na(exist)], collapse="' and '"), "' in argument 'expert'.", sep=''))
#
#    expertenv <- list2env(expert) 
#
#    if (exists("parameterization", expertenv)) {
#      if (!is.character(expert[["parameterization"]]))
#        stop("Argument 'parameterization' must be either a vector of 'centered', 'noncentered' values or a character string of form 'asis#' with # a positive integer.")
#      nmatches <- grep("^asis[1-9][0-9]*$", expert[["parameterization"]])
#      if (length(nmatches) == 0) {
#        parameterization <- match.arg(expert[["parameterization"]], strategies, several.ok = TRUE)
#      } else if (length(nmatches) == 1) {
#        parameterization <- rep(strategies, nmatches)
#      } else {
#        parameterization <- NA
#      }
#      if (!all(parameterization %in% strategies)) {
#        stop("Argument 'parameterization' must be either a vector of 'centered', 'noncentered' values or a character string of form 'asis#' with # a positive integer.")
#      }
#    } else {
#      parameterization <- expertdefault$parameterization
#    }
#
#    if (exists("mhcontrol", expertenv)) {
#      mhcontrol <- expert[["mhcontrol"]]
#      if (!is.list(mhcontrol))  # TODO write proper validation
#        stop("Argument 'mhcontrol' must be a list.")
#    } else {
#      mhcontrol <- expertdefault$mhcontrol
#    }
#    use.mala <- mhcontrol$use.mala
#    if (!isTRUE(use.mala)) {
#      use.mala <- FALSE
#    }
#
#    if (exists("gammaprior", expertenv)) {
#      gammaprior <- expert[["gammaprior"]]
#      if (!is.logical(gammaprior)) stop("Argument 'gammaprior' must be TRUE or FALSE.")
#    } else {
#      gammaprior <- expertdefault$gammaprior
#    }
#
#    if (exists("correct.latent.draws", expertenv)) {
#      correct.latent.draws <- expert[["correct.latent.draws"]]
#      if (!is.logical(correct.latent.draws)) stop("Argument 'correct.latent.draws' must be TRUE or FALSE.")
#    } else {
#      correct.latent.draws <- expertdefault$correct.latent.draws
#    }
#  }
#
#  renameparam <- c("centered" = "C", "noncentered" = "NC")
#  if (!quiet) {
#    cat(paste("\nCalling ", asisprint(renameparam[parameterization], renameparam), " MCMC sampler with ", draws+burnin, " iter. Series length is ", length(y), ".\n",sep=""), file=stderr())
#    flush.console()
#  }
#
#  if (.Platform$OS.type != "unix") myquiet <- TRUE else myquiet <- quiet  # Hack to prevent console flushing problems with Windows
#
#  runtime <- system.time({
#    res <- svlsample_cpp(y, draws, burnin, designmatrix, thinpara, thinlatent, thintime,
#                                 startpara, startlatent,
#                                 priorphi[1], priorphi[2], priorrho[1], priorrho[2],
#                                 0.5, 0.5/priorsigma, priormu[1], priormu[2],
#                                 priorbeta[1], priorbeta[2], !myquiet,
#                                 myoffset, use.mala,
#                                 gammaprior, correct.latent.draws,
#                                 parameterization, FALSE)
#  })
#
#  if (any(is.na(res))) stop("Sampler returned NA. This is most likely due to bad input checks and shouldn't happen. Please report to package maintainer.")
#
#  if (!quiet) {
#    cat("Timing (elapsed): ", file=stderr())
#    cat(runtime["elapsed"], file=stderr())
#    cat(" seconds.\n", file=stderr())
#    cat(round((draws+burnin)/runtime[3]), "iterations per second.\n\n", file=stderr())
#    cat("Converting results to coda objects... ", file=stderr())
#  }
#  
#  # create svldraws class
#  colnames(res$para) <- c("mu", "phi", "sigma", "rho")
#  if (keeptime == "all")
#    colnames(res$latent) <- paste0('h_', arorder + seq_along(y))
#  else if (keeptime == "last")
#    colnames(res$latent) <- paste0('h_', arorder + length(y))
#  res$runtime <- runtime
#  res$y <- y_orig
#  res$para <- coda::mcmc(res$para, burnin+thinpara, burnin+draws, thinpara)
#  res$latent <- coda::mcmc(res$latent, burnin+thinlatent, burnin+draws, thinlatent)
#  res$latent0 <- coda::mcmc(res$latent0, burnin+thinlatent, burnin+draws, thinlatent)
#  res$thinning <- list(para = thinpara, latent = thinlatent, time = keeptime)
#  res$priors <- list(mu = priormu, phi = priorphi, sigma = priorsigma, rho = priorrho, gammaprior = gammaprior)
#  if (!any(is.na(designmatrix))) {
#    res$beta <- coda::mcmc(res$beta, burnin+thinpara, burnin+draws, thinpara)
#    colnames(res$beta) <- paste0("b_", 0:(NCOL(designmatrix)-1))
#    res$priors <- c(res$priors, "beta" = list(priorbeta), "designmatrix" = list(designmatrix))
#  } else {
#    res$beta <- NULL
#  }
#  res$meanmodel <- meanmodel
#  class(res) <- c("svldraws", "svdraws")
#
#  if (!quiet) {
#    cat("Done!\n", file=stderr())
#    cat("Summarizing posterior draws... ", file=stderr())
#  }
#  res <- updatesummary(res, ...)
#
#  if (!quiet) cat("Done!\n\n", file=stderr())
#  res
#}
#
##' Minimal overhead version of \code{\link{svlsample}}.
##' 
##' \code{svlsample2} is a minimal overhead version of \code{\link{svlsample}}
##' with slightly different default arguments and a simplified return value
##' structure. It is intended to be used mainly for one-step updates where speed
##' is an issue, e.g., as a plug-in into other MCMC samplers. Note that
##' absolutely no input checking is performed, thus this function is to be used
##' with proper care!
##' 
##' As opposed to the ordinary \code{\link{svlsample}}, the default values differ
##' for \code{draws}, \code{burnin}, and \code{quiet}. Note that currently
##' neither \code{expert} nor \code{\dots{}} arguments are provided.
##' 
##' @param y numeric vector containing the data (usually log-returns), which
##' must not contain zeros. Alternatively, \code{y} can be an \code{svsim}
##' object. In this case, the returns will be extracted and a warning is thrown.
##' @param draws single number greater or equal to 1, indicating the number of
##' draws after burn-in (see below). Will be automatically coerced to integer.
##' The default value is 1.
##' @param burnin single number greater or equal to 0, indicating the number of
##' draws discarded as burn-in. Will be automatically coerced to integer. The
##' default value is 0.
##' @param priormu numeric vector of length 2, indicating mean and standard
##' deviation for the Gaussian prior distribution of the parameter \code{mu},
##' the level of the log-volatility. The default value is \code{c(0, 100)},
##' which constitutes a practically uninformative prior for common exchange rate
##' datasets, stock returns and the like.
##' @param priorphi numeric vector of length 2, indicating the shape parameters
##' for the Beta prior distribution of the transformed parameter
##' \code{(phi + 1) / 2}, where \code{phi} denotes the persistence of the
##' log-volatility. The default value is \code{c(5, 1.5)}, which constitutes a
##' prior that puts some belief in a persistent log-volatility but also
##' encompasses the region where \code{phi} is around 0.
##' @param priorsigma single positive real number, which stands for the scaling
##' of the transformed parameter \code{sigma^2}, where \code{sigma} denotes the
##' volatility of log-volatility. More precisely, \code{sigma^2 ~ priorsigma *
##' chisq(df = 1)}. The default value is \code{1}, which constitutes a
##' reasonably vague prior for many common exchange rate datasets, stock returns
##' and the like.
##' @param priorrho numeric vector of length 2, indicating the shape parameters
##' for the Beta prior distribution of the transformed parameter
##' \code{(rho + 1) / 2}, where \code{rho} denotes the conditional correlation
##' between observation and the increment of the
##' log-volatility. The default value is \code{c(4, 4)}, which constitutes a
##' slightly informative prior around 0 (the no leverage case) to boost convergence.
##' @param priorspec to be documented TODO
##' @param thinpara single number greater or equal to 1, coercible to integer.
##' Every \code{thinpara}th parameter draw is kept and returned. The default
##' value is 1, corresponding to no thinning of the parameter draws i.e. every
##' draw is stored.
##' @param thinlatent single number greater or equal to 1, coercible to integer.
##' Every \code{thinlatent}th latent variable draw is kept and returned. The
##' default value is 1, corresponding to no thinning of the latent variable
##' draws, i.e. every draw is kept.
##' @param keeptime Either 'all' (the default) or 'last'. Indicates which latent
##  volatility draws should be stored.
##' @param quiet logical value indicating whether the progress bar and other
##' informative output during sampling should be omitted. The default value is
##' \code{TRUE}.
##' @param startpara \emph{compulsory} named list, containing the starting values
##' for the parameter draws. It must contain four
##' elements named \code{mu}, \code{phi}, \code{sigma}, and \code{rho}, where \code{mu} is
##' an arbitrary numerical value, \code{phi} is a real number between \code{-1}
##' and \code{1}, \code{sigma} is a positive real number, and \code{rho} is
##' a real number between \code{-1} and \code{1}.
##' @param startlatent \emph{compulsory} vector of length \code{length(y)},
##' containing the starting values for the latent log-volatility draws.
##' @return The value returned is a list object holding
##' \item{para}{matrix of dimension \code{4 x draws} containing
##' the \emph{parameter} draws from the posterior distribution.}
##' \item{latent}{matrix of dimension \code{length(y) x draws} containing the
##' \emph{latent instantaneous log-volatility} draws from the posterior
##' distribution.}
##' \item{latent0}{Vector of
##' length \code{draws} containing the draw(s) of the initial latent variable
##' \code{h_0}.}
##' \item{meanmodel}{always equals \code{"none"}}
##' @author Darjus Hosszejni \email{darjus.hosszejni@@wu.ac.at}
##' @seealso \code{\link{svlsample}}
##' @keywords models ts
##' @examples
##' data(exrates)
##' aud.price <- subset(exrates,
##'   as.Date("2010-01-01") <= date & date < as.Date("2011-01-01"),
##'   "AUD")[,1]
##' draws <- svlsample2(logret(aud.price),
##'                     draws = 10, burnin = 0,
##'                     startpara = list(phi=0.95, mu=-10, sigma=0.2, rho=-0.1),
##'                     startlatent = rep_len(-10, length(aud.price)-1))
##' @export
#svlsample2 <- function(y, draws = 10000, burnin = 1000, designmatrix = NA,
#                       priormu = c(0, 100), priorphi = c(5, 1.5), priorsigma = 1,
#                       priornu = 0, priorrho = c(4, 4),
#                       priorbeta = c(0, 10000), priorlatent0 = "stationary",
#                       priorspec = NULL,
#                       thinpara = 1, thinlatent = 1, keeptime = "all", keeptau = FALSE,
#                       quiet = FALSE, startpara = NULL, startlatent = NULL, expert = NULL, ...) {
#  # TODO make minimal overhead
#  svsample(y, draws = draws, burnin = burnin, designmatrix = designmatrix,
#           priormu = priormu, priorphi = priorphi, priorsigma = priorsigma,
#           priornu = priornu, priorrho = priorrho,
#           priorbeta = priorbeta, priorlatent0 = priorlatent0,
#           priorspec = priorspec,
#           thinpara = thinpara, thinlatent = thinlatent, keeptime = keeptime,
#           keeptau = keeptau, quiet = quiet, startpara = startpara, startlatent = startlatent,
#           expert = expert, ...)
#}
#
#svlsample2.old <- function(y, draws = 1, burnin = 0,
#                       priormu = c(0, 100), priorphi = c(5, 1.5), priorsigma = 1, priorrho = c(4, 4),
#                       thinpara = 1, thinlatent = 1, keeptime = "all",
#                       quiet = TRUE, startpara, startlatent) {
#
#  if (keeptime == "all") thintime <- 1L else if (keeptime == "last") thintime <- length(y) else stop("Parameter 'keeptime' must be either 'all' or 'last'.")
#
#  res <- svlsample_cpp(y, draws, burnin, matrix(NA), thinpara, thinlatent, thintime,
#                               startpara, startlatent,
#                               priorphi[1], priorphi[2], priorrho[1], priorrho[2],
#                               0.5, 0.5/priorsigma, priormu[1], priormu[2],
#                               0, 1, !quiet,
#                               0, FALSE, TRUE,
#                               FALSE, rep(c("centered", "noncentered"), 5), FALSE)
#
#  res$para <- t(res$para)
#  res$latent <- t(res$latent)
#  res$meanmodel <- "none"
#  res
#}

