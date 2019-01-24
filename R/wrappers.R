# R wrapper function for the main MCMC loop



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
#' object. In this case, the returns will be extracted and a warning is thrown.
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
#' \code{(phi+1)/2}, where \code{phi} denotes the persistence of the
#' log-volatility. The default value is \code{c(5, 1.5)}, which constitutes a
#' prior that puts some belief in a persistent log-volatility but also
#' encompasses the region where \code{phi} is around 0.
#' @param priorsigma single positive real number, which stands for the scaling
#' of the transformed parameter \code{sigma^2}, where \code{sigma} denotes the
#' volatility of log-volatility. More precisely, \code{sigma^2 ~ priorsigma *
#' chisq(df = 1)}. The default value is \code{1}, which constitutes a
#' reasonably vague prior for many common exchange rate datasets, stock returns
#' and the like.
#' @param priornu numeric vector of length 2 (or \code{NA}), indicating the
#' lower and upper bounds for the uniform prior distribution of the parameter
#' \code{nu}, the degrees-of-freedom parameter of the conditional innovations
#' t-distribution. The default value is \code{NA}, fixing the
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
#' @param thinpara single number greater or equal to 1, coercible to integer.
#' Every \code{thinpara}th parameter draw is kept and returned. The default
#' value is 1, corresponding to no thinning of the parameter draws i.e. every
#' draw is stored.
#' @param thinlatent single number greater or equal to 1, coercible to integer.
#' Every \code{thinlatent}th latent variable draw is kept and returned. The
#' default value is 1, corresponding to no thinning of the latent variable
#' draws, i.e. every draw is kept.
#' @param thintime \emph{deprecated}
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
#' and \code{1}, and \code{sigma} is a positive real number. The default value
#' is \code{list(mu = -10, phi = 0.9, sigma = 0.3)}. Moreover, if
#' \code{priornu} is not \code{NA}, \code{startpara} must also contain an
#' element named \code{nu} (the degrees of freedom parameter for the
#' t-innovations).
#' @param startlatent \emph{optional} vector of length \code{length(y)},
#' containing the starting values for the latent log-volatility draws. The
#' default value is \code{rep(-10, length(y))}.
#' @param expert \emph{optional} named list of expert parameters. For most
#' applications, the default values probably work best. Interested users are
#' referred to the literature provided in the References section. If
#' \code{expert} is provided, it may contain the following named elements:
#' 
#' \code{parameterization}: Character string equal to \code{"centered"},
#' \code{"noncentered"}, \code{"GIS_C"}, or \code{"GIS_NC"}. Defaults to
#' \code{"GIS_C"}.
#' 
#' \code{mhcontrol}: Single numeric value controlling the proposal density of a
#' Metropolis-Hastings (MH) update step when sampling \code{sigma}. If
#' \code{mhcontrol} is smaller than 0, an independence proposal will be used,
#' while values greater than zero control the stepsize of a log-random-walk
#' proposal. Defaults to \code{-1}.
#' 
#' \code{gammaprior}: Single logical value indicating whether a Gamma prior for
#' \code{sigma^2} should be used. If set to \code{FALSE}, an Inverse Gamma
#' prior is employed. Defaults to \code{TRUE}.
#' 
#' \code{truncnormal}: Single logical value indicating whether a truncated
#' Gaussian distribution should be used as proposal for draws of \code{phi}. If
#' set to \code{FALSE}, a regular Gaussian prior is employed and the draw is
#' immediately discarded when values outside the unit ball happen to be drawn.
#' Defaults to \code{FALSE}.
#' 
#' \code{mhsteps}: Either \code{1}, \code{2}, or \code{3}. Indicates the number
#' of blocks used for drawing from the posterior of the parameters. Defaults to
#' \code{2}.
#' 
#' \code{proposalvar4sigmaphi}: Single positive number indicating the
#' conditional prior variance of \code{sigma*phi} in the ridge \emph{proposal}
#' density for sampling \code{(mu, phi)}. Defaults to \code{10^8}.
#' 
#' \code{proposalvar4sigmatheta}: Single positive number indicating the
#' conditional prior variance of \code{sigma*theta} in the ridge
#' \emph{proposal} density for sampling \code{(mu, phi)}. Defaults to
#' \code{10^12}.
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
#' \code{thintime}.}
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
#' @seealso \code{\link{svsim}}, \code{\link{svlsample}}, \code{\link{updatesummary}},
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
#' plot(draws)
#' 
#' 
#' #\dontrun{
#' # Example 2
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
#' # Example 3
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
#' #}
#' @export
svsample <- function(y, draws = 10000, burnin = 1000, designmatrix = NA,
                     priormu = c(0, 100), priorphi = c(5, 1.5), priorsigma = 1,
                     priornu = NA, priorbeta = c(0, 10000), priorlatent0 = "stationary",
                     thinpara = 1, thinlatent = 1, thintime = 1,
                     keeptau = FALSE, quiet = FALSE, startpara, startlatent, expert, ...) {

  # Some error checking for y
  if (inherits(y, "svsim")) {
    y <- y[["y"]]
    warning("Extracted data vector from 'svsim'-object.")
  }
  if (!is.numeric(y)) stop("Argument 'y' (data vector) must be numeric.")

  if (length(y) < 2) stop("Argument 'y' (data vector) must contain at least two elements.")

  if (any(is.na(designmatrix)) && any(y^2 == 0)) {
    myoffset <- sd(y)/10000
    warning(paste("Argument 'y' (data vector) contains zeros. I am adding an offset constant of size ", myoffset, " to do the auxiliary mixture sampling. If you want to avoid this, you might consider de-meaning the returns before calling this function.", sep=""))
  } else myoffset <- 0

  # Some error checking for draws
  if (!is.numeric(draws) || length(draws) != 1 || draws < 1) {
    stop("Argument 'draws' (number of MCMC iterations after burn-in) must be a single number >= 1.")
  } else {
    draws <- as.integer(draws)
  }

  # Some error checking for burnin
  if (!is.numeric(burnin) || length(burnin) != 1 || burnin < 0) {
    stop("Argument 'burnin' (burn-in period) must be a single number >= 0.")
  } else {
    burnin <- as.integer(burnin)
  }

  # Some error checking for designmatrix and save the mode of mean modelling
  meanmodel <- "matrix"
  if (any(is.na(designmatrix))) {
    designmatrix <- matrix(NA)
    meanmodel <- "none"
  } else {
    if (any(grep("ar[0-9]+$", as.character(designmatrix)[1]))) {
      arorder <- as.integer(gsub("ar", "", as.character(designmatrix)))
      if (length(y) <= arorder + 1) stop("Time series 'y' is too short for this AR process.")
      designmatrix <- matrix(rep(1, length(y) - arorder), ncol = 1)
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
    }
    if (!is.numeric(designmatrix)) stop("Argument 'designmatrix' must be a numeric matrix or an AR-specification.")
    if (!is.matrix(designmatrix)) {
      designmatrix <- matrix(designmatrix, ncol = 1)
    }
    if (!nrow(designmatrix) == length(y)) stop("Number of columns of argument 'designmatrix' must be equal to length(y).")
  }

  # Some error checking for the prior parameters 
  if (!is.numeric(priormu) || length(priormu) != 2) {
    stop("Argument 'priormu' (mean and sd for the Gaussian prior for mu) must be numeric and of length 2.")
  }

  if (!is.numeric(priorphi) | length(priorphi) != 2) {
    stop("Argument 'priorphi' (shape1 and shape2 parameters for the Beta prior for (phi+1)/2) must be numeric and of length 2.")
  }

  if (!is.numeric(priorsigma) | length(priorsigma) != 1 | priorsigma <= 0) {
    stop("Argument 'priorsigma' (scaling of the chi-squared(df = 1) prior for sigma^2) must be a single number > 0.")
  }

  if (any(is.na(priornu))) {
    priornu <- NA
  } else {
    if (!is.numeric(priornu) || length(priornu) != 2 || priornu[1] >= priornu[2] || priornu[1] < 0) {
      stop("If not NA, argument 'priornu' (lower and upper bounds for the uniform prior for the df) must be numeric and of length 2. Moreover, 0 <= priornu[1] < priornu[2].")
    }
  }

  if (!is.numeric(priorbeta) || length(priorbeta) != 2) {
    stop("Argument 'priorbeta' (means and sds for the independent Gaussian priors for beta) must be numeric and of length 2.")
  }

  if (!is.numeric(priorlatent0) || length(priorlatent0) != 1 || priorlatent0 < 0) {
    if (priorlatent0 == "stationary") priorlatent0 <- -1L else
      stop("Argument 'priorlatent0' must be 'stationary' or a single non-negative number.")
  }

  # Some error checking for thinpara
  if (!is.numeric(thinpara) || length(thinpara) != 1 || thinpara < 1) {
    stop("Argument 'thinpara' (thinning parameter for mu, phi, and sigma) must be a single number >= 1.")
  } else {
    thinpara <- as.integer(thinpara)
  }

  # Some error checking for thinlatent
  if (!is.numeric(thinlatent) || length(thinlatent) != 1 || thinlatent < 1) {
    stop("Argument 'thinlatent' (thinning parameter for the latent log-volatilities) must be a single number >= 1.")
  } else {
    thinlatent <- as.integer(thinlatent)
  }

  # thintime deprecated
  if (!(identical(thintime, 1L) || identical(thintime, 1.0))) {
    thintime <- 1L
    warning("Parameter 'thintime' is deprecated. Setting 'thintime' to 1.")
  }

  # Some error checking for expert
  if (missing(expert)) {
    para <- 3L ; parameterization <- 'GIS_C'
    mhcontrol <- -1
    gammaprior <- TRUE
    truncnormal <- FALSE
    mhsteps <- 2L
    B011 <- 10^8
    B022 <- 10^12
  } else {
    expertnames <- names(expert)
    if (!is.list(expert) | is.null(expertnames) | any(expertnames == ""))
      stop("Argument 'expert' must be a named list with nonempty names.")
    if (length(unique(expertnames)) != length(expertnames))
      stop("No duplicate elements allowed in argument 'expert'.")
    allowednames <- c("parameterization", "mhcontrol", "gammaprior", "truncnormal", "mhsteps", "proposalvar4sigmaphi", "proposalvar4sigmatheta")
    exist <- pmatch(expertnames, allowednames)
    if (any(is.na(exist)))
      stop(paste("Illegal element '", paste(expertnames[is.na(exist)], collapse="' and '"), "' in argument 'expert'.", sep=''))

    expertenv <- list2env(expert) 

    if (exists("parameterization", expertenv)) {
      parameterization <- expert[["parameterization"]]
      if (!is.character(parameterization) | is.na(parameterization)) {
        stop("Argument 'parameterization' must be either 'centered', 'noncentered', 'GIS_C', or 'GIS_NC'.")
      }
      switch(parameterization,
             centered = para <- 1L,
             noncentered = para <- 2L,
             GIS_C = para <- 3L,
             GIS_NC = para <- 4L,
             stop("Unknown parameterization. Currently you can only use 'centered', 'noncentered', 'GIS_C', and 'GIS_NC'.")
             )
    } else {
      para <- 3L ; parameterization <- 'GIS_C'
    }

    # Remark: mhcontrol < 0 means independence proposal,
    #         mhcontrol > 0 controls stepsize of log-random-walk proposal
    if (exists("mhcontrol", expertenv)) {
      mhcontrol <- expert[["mhcontrol"]]
      if (!is.numeric(mhcontrol) | length(mhcontrol) != 1)
        stop("Argument 'mhcontrol' must be a single number.")
    } else {
      mhcontrol <- -1
    }

    # use a Gamma prior for sigma^2 in C?
    if (exists("gammaprior", expertenv)) {
      gammaprior <- expert[["gammaprior"]]
      if (!is.logical(gammaprior)) stop("Argument 'gammaprior' must be TRUE or FALSE.")
    } else {
      gammaprior <- TRUE
    }

    # use a truncated normal as proposal? (or normal with rejection step)
    if (exists("truncnormal", expertenv)) {
      truncnormal <- expert[["truncnormal"]]
      if (!is.logical(truncnormal)) stop("Argument 'truncnormal' must be TRUE or FALSE.")
    } else {
      truncnormal <- FALSE
    }

    if (exists("mhsteps", expertenv)) {
      mhsteps <- as.integer(expert[["mhsteps"]])
      if (mhsteps != 2L & mhsteps != 1L & mhsteps != 3L) stop("mhsteps must be 1, 2, or 3")
      if (mhsteps != 2L & mhcontrol >= 0)
        stop("Log normal random walk proposal currently only implemented for mhsteps==2.")
      if (mhsteps != 2L & !isTRUE(gammaprior))
        stop("Inverse Gamma prior currently only implemented for mhsteps==2.")
    } else {
      mhsteps <- 2L
    }

    # prior for ridge _proposal_ (variance of sigma*phi)
    if (exists("proposalvar4sigmaphi", expertenv)) {
      B011 <- expert[["proposalvar4sigmaphi"]]
      if (!is.numeric(B011) | length(B011) != 1 | B011 <= 0)
        stop("Argument 'proposalvar4sigmaphi' must be a positive number.")
    } else {
      B011 <- 10^8
    }

    # prior for ridge _proposal_ (variance of sigma*theta)
    if (exists("proposalvar4sigmatheta", expertenv)) {
      B022 <- expert[["proposalvar4sigmatheta"]]
      if (!is.numeric(B022) | length(B022) != 1 | B022 <= 0)
        stop("Argument 'proposalvar4sigmatheta' must be a positive number.")
    } else {
      B022 <- 10^12
    }
  }

  # Some input checking for startpara
  if (missing(startpara)) {
    if (any(is.na(priornu))) {
      startpara <- list(mu = priormu[1],
                        phi = 2 * (priorphi[1] / sum(priorphi)) - 1,
                        sigma = priorsigma)
    } else {
      startpara <- list(mu = priormu[1],
                        phi = 2 * (priorphi[1] / sum(priorphi)) - 1,
                        sigma = priorsigma,
                        nu = mean(priornu))
    }
  } else {
    if (!is.list(startpara))
      stop("Argument 'startpara' must be a list. Its elements must be named 'mu', 'phi', 'sigma'. Moreover, if !is.na(priornu), an element named 'nu' must exist.")

    if (!is.numeric(startpara[["mu"]]))
      stop('Argument \'startpara[["mu"]]\' must exist and be numeric.')

    if (!is.numeric(startpara[["phi"]]))
      stop('Argument \'startpara[["phi"]]\' must exist and be numeric.')

    if (abs(startpara[["phi"]]) >= 1)
      stop('Argument \'startpara[["phi"]]\' must be between -1 and 1.')

    if (!is.numeric(startpara[["sigma"]]))
      stop('Argument \'startpara[["sigma"]]\' must exist and be numeric.')

    if (startpara[["sigma"]] <= 0)
      stop('Argument \'startpara[["sigma"]]\' must be positive.')

    if (!is.na(priornu) && !is.numeric(startpara[["nu"]]))
      stop('Argument \'startpara[["nu"]]\' must exist and be numeric.')

    if (!is.na(priornu) && (startpara[["nu"]] > priornu[2] || startpara[["nu"]] < priornu[1]))
      stop('Argument \'startpara[["nu"]]\' must be within range(priornu).')
  }

  # Some input checking for startlatent
  if (missing(startlatent)) {
    startlatent <- rep(-10, length(y))
  } else {
    if (!is.numeric(startlatent) | length(startlatent) != length(y))
      stop("Argument 'startlatent' must be numeric and of the same length as the data 'y'.")
  }

  if (length(keeptau) != 1 || !is.logical(keeptau)) {
    stop("Argument 'keeptau' must be TRUE or FALSE.")
  }

  if (is.na(priornu) && keeptau) {
    warning("Setting argument 'keeptau' to FALSE, as 'priornu' is NA.")
    keeptau <- FALSE
  }

  if (!quiet) {
    cat(paste("\nCalling ", parameterization, " MCMC sampler with ", draws+burnin, " iter. Series length is ", length(y), ".\n",sep=""), file=stderr())
    flush.console()
  }

  if (.Platform$OS.type != "unix") myquiet <- TRUE else myquiet <- quiet  # Hack to prevent console flushing problems with Windows

  runtime <- system.time(res <-
    svsample_cpp(y, draws, burnin, designmatrix,
          priormu[1], priormu[2]^2, priorphi[1], priorphi[2], priorsigma, 
          thinlatent, thintime, startpara, startlatent, keeptau, myquiet, para,
          mhsteps, B011, B022, mhcontrol, gammaprior, truncnormal,
          myoffset, FALSE, priornu, priorbeta, priorlatent0))

  if (any(is.na(res))) stop("Sampler returned NA. This is most likely due to bad input checks and shouldn't happen. Please report to package maintainer.")

  if (!quiet) {
    cat("Timing (elapsed): ", file=stderr())
    cat(runtime["elapsed"], file=stderr())
    cat(" seconds.\n", file=stderr())
    cat(round((draws+burnin)/runtime[3]), "iterations per second.\n\n", file=stderr())
    cat("Converting results to coda objects... ", file=stderr())
  }

  # store results:
  # remark: +1, because C-sampler also returns the first value
  res$y <- y
  res$para <- mcmc(res$para[seq(burnin+thinpara+1, burnin+draws+1, thinpara),,drop=FALSE], burnin+thinpara, burnin+draws, thinpara)
  res$latent <- mcmc(t(res$latent), burnin+thinlatent, burnin+draws, thinlatent)
  attr(res$latent, "dimnames") <- list(NULL, paste('h_', seq(1, length(y), by=thintime), sep=''))
  res$latent0 <- mcmc(res$latent0, burnin+thinlatent, burnin+draws, thinlatent)
  if (!any(is.na(designmatrix))) {
    res$beta <- mcmc(res$beta[seq(burnin+thinpara+1, burnin+draws+1, thinpara),,drop=FALSE], burnin+thinpara, burnin+draws, thinpara)
    attr(res$beta, "dimnames") <- list(NULL, paste("b", 0:(ncol(designmatrix)-1), sep = "_"))
  } else res$beta <- NULL
  res$meanmodel <- meanmodel

  if (ncol(res$para) == 3) {
    attr(res$para, "dimnames") <- list(NULL, c("mu", "phi", "sigma"))
    res$priors <- list(mu = priormu, phi = priorphi, sigma = priorsigma, gammaprior = gammaprior)
  } else {
    attr(res$para, "dimnames") <- list(NULL, c("mu", "phi", "sigma", "nu"))
    res$priors <- list(mu = priormu, phi = priorphi, sigma = priorsigma, nu = priornu, gammaprior = gammaprior)
    #res$tau <- mcmc(t(res$tau), burnin+thinlatent, burnin+draws, thinlatent)
  }

  if (!any(is.na(designmatrix))) {
    res$priors <- c(res$priors, "beta" = list(priorbeta), "designmatrix" = list(designmatrix))
  }

  if (keeptau) {
    res$tau <- mcmc(t(res$tau), burnin+thinlatent, burnin+draws, thinlatent)
    attr(res$tau, "dimnames") <- list(NULL, paste('tau_', seq(1, length(y), by=thintime), sep=''))
  }

  res$runtime <- runtime
  res$thinning <- list(para = thinpara, latent = thinlatent, time = thintime)
  class(res) <- "svdraws"

  if (!quiet) {
    cat("Done!\n", file=stderr())
    cat("Summarizing posterior draws... ", file=stderr())
  }
  res <- updatesummary(res, ...)

  if (!quiet) cat("Done!\n\n", file=stderr())
  res
}

# This function does not check input nor converts the result to coda objects!



#' Minimal overhead version of \code{\link{svsample}}.
#' 
#' \code{svsample2} is a minimal overhead version of \code{\link{svsample}}
#' with slightly different default arguments and a simplified return value
#' structure. It is intended to be used mainly for one-step updates where speed
#' is an issue, e.g., as a plug-in into other MCMC samplers. Note that
#' absolutely no input checking is performed, thus this function is to be used
#' with proper care!
#' 
#' As opposed to the ordinary \code{\link{svsample}}, the default values differ
#' for \code{draws}, \code{burnin}, and \code{quiet}. Note that currently
#' neither \code{expert} nor \code{\dots{}} arguments are provided.
#' 
#' @param y numeric vector containing the data (usually log-returns), which
#' must not contain zeroes.
#' @param draws single number greater or equal to 1, indicating the number of
#' draws after burn-in (see below). Will be automatically coerced to integer.
#' The defaults value is 1.
#' @param burnin single number greater or equal to 0, indicating the number of
#' draws discarded as burn-in. Will be automatically coerced to integer. The
#' default value is 0.
#' @param priormu numeric vector of length 2, indicating mean and standard
#' deviation for the Gaussian prior distribution of the parameter \code{mu},
#' the level of the log-volatility. The default value is \code{c(0, 100)},
#' which constitutes a practically uninformative prior for common exchange rate
#' datasets, stock returns and the like.
#' @param priorphi numeric vector of length 2, indicating the shape parameters
#' for the Beta prior distribution of the transformed parameter
#' \code{(phi+1)/2}, where \code{phi} denotes the persistence of the
#' log-volatility. The default value is \code{c(5, 1.5)}, which constitutes a
#' prior that puts some belief in a persistent log-volatility but also
#' encompasses the region where \code{phi} is around 0.
#' @param priorsigma single positive real number, which stands for the scaling
#' of the transformed parameter \code{sigma^2}, where \code{sigma} denotes the
#' volatility of log-volatility. More precisely, \code{sigma^2 ~ priorsigma *
#' chisq(df = 1)}. The default value is \code{1}, which constitutes a
#' reasonably vague prior for many common exchange rate datasets, stock returns
#' and the like.
#' @param priornu numeric vector of length 2 (or \code{NA}), indicating the
#' lower and upper bounds for the uniform prior distribution of the parameter
#' \code{nu}, the degrees-of-freedom parameter of the conditional innovations
#' t-distribution. The default value is \code{NA}, fixing the
#' degrees-of-freedom to infinity. This corresponds to conditional standard
#' normal innovations, the pre-1.1.0 behavior of \pkg{stochvol}.
#' @param priorlatent0 either a single non-negative number or the string
#' \code{'stationary'} (the default, also the behavior before version 1.3.0).
#' When \code{priorlatent0} is equal to \code{'stationary'}, the stationary
#' distribution of the latent AR(1)-process is used as the prior for the
#' initial log-volatility \code{h_0}. When \code{priorlatent0} is equal to a
#' number \eqn{B}, we have \eqn{h_0 \sim N(\mu, B\sigma^2)} a priori.
#' @param thinpara single number greater or equal to 1, coercible to integer.
#' Every \code{thinpara}th parameter draw is kept and returned. The default
#' value is 1, corresponding to no thinning of the parameter draws -- every
#' draw is stored.
#' @param thinlatent single number greater or equal to 1, coercible to integer.
#' Every \code{thinlatent}th latent variable draw is kept and returned. The
#' default value is 1, corresponding to no thinning of the latent variable
#' draws, i.e. every draw is kept.
#' @param thintime single number greater or equal to 1, coercible to integer.
#' If \code{thintime} is different from 1, only every \code{thintime}th latent
#' log-volatility is being monitored. If, e.g., \code{thintime = 3}, the latent
#' log-volatilities \code{h_1,h_4,h_7,...} will be kept. The default value is
#' 1, meaning that all latent variables \code{h_1,h_2,h_3,...} are stored.
#' @param keeptau logical value indicating whether the 'variance inflation
#' factors' should be stored (used for the sampler with conditional t
#' innovations only). This may be useful to check at what point(s) in time the
#' normal disturbance had to be 'upscaled' by a mixture factor and when the
#' series behaved 'normally'.
#' @param quiet logical value indicating whether the progress bar and other
#' informative output during sampling should be omitted. The default value is
#' \code{TRUE}, implying non-verbose output.
#' @param startpara \emph{compulsory} named list, containing the starting
#' values for the parameter draws. \code{startpara} must contain three elements
#' named \code{mu}, \code{phi}, and \code{sigma}, where \code{mu} is an
#' arbitrary numerical value, \code{phi} is a real number between \code{-1} and
#' \code{1}, and \code{sigma} is a positive real number. Moreover, if
#' \code{priornu} is not \code{NA}, \code{startpara} must also contain an
#' element named \code{nu} (the degrees of freedom parameter for the
#' t-innovations).
#' @param startlatent \emph{compulsory} vector of length \code{length(x$y)},
#' containing the starting values for the latent log-volatility draws.
#' @return A list with three components:
#' \item{para}{\code{3} times
#' \code{draws} matrix containing the parameter draws. If \code{priornu} is not
#' \code{NA}, this is a \code{4} times \code{draws} matrix.}
#' \item{latent}{\code{length(y)} times \code{draws} matrix containing draws of
#' the latent variables \code{h_1, \dots{}, h_n}.}
#' \item{latent0}{Vector of
#' length \code{draws} containing the draw(s) of the initial latent variable
#' \code{h_0}.}
#' @note Please refer to the package vignette for an example.
#' @section Warning: Expert use only! For most applications, the use of
#' \code{\link{svsample}} is recommended.
#' @author Gregor Kastner \email{gregor.kastner@@wu.ac.at}
#' @seealso \code{\link{svsample}}
#' @keywords models ts
#' @export
svsample2 <- function(y, draws = 1, burnin = 0, priormu = c(0, 100), priorphi = c(5, 1.5), priorsigma = 1, priornu = NA, priorlatent0 = "stationary", thinpara = 1, thinlatent = 1, thintime = 1, keeptau = FALSE, quiet = TRUE, startpara, startlatent) {

 if (priorlatent0 == "stationary") priorlatent0 <- -1L

  # thintime deprecated
  if (!(identical(thintime, 1L) || identical(thintime, 1.0))) {
    thintime <- 1L
    warning("Parameter 'thintime' is deprecated. Setting 'thintime' to 1.")
  }

 res <- svsample_cpp(y, draws, burnin, matrix(NA), priormu[1], priormu[2]^2,
	      priorphi[1], priorphi[2], priorsigma, thinlatent,
	      thintime, startpara, startlatent, keeptau, quiet, 3L, 2L, 10^8,
	      10^12, -1, TRUE, FALSE, 0, FALSE, priornu, c(NA, NA), priorlatent0)

 res$para <- t(res$para[seq(burnin+thinpara+1, burnin+draws+1, thinpara),,drop=FALSE])
 if (nrow(res$para) == 3) {
  rownames(res$para) <- names(res$para) <- c("mu", "phi", "sigma")
 } else {
  rownames(res$para) <- names(res$para) <- c("mu", "phi", "sigma", "nu")
 }
 res.meanmodel <- "none"

 res
}

#' Markov Chain Monte Carlo (MCMC) Sampling for the Stochastic Volatility
#' Model with Leverage (SVL)
#' 
#' \code{svlsample} simulates from the joint posterior distribution of the SVL
#' parameters \code{mu}, \code{phi}, \code{sigma}, and \code{rho},
#' along with the latent log-volatilities \code{h_1,...,h_n} and returns the
#' MCMC draws. If a design matrix is provided, simple Bayesian regression can
#' also be conducted.
#' 
#' @param y numeric vector containing the data (usually log-returns), which
#' must not contain zeros. Alternatively, \code{y} can be an \code{svsim}
#' object. In this case, the returns will be extracted and a warning is thrown.
#' @param draws single number greater or equal to 1, indicating the number of
#' draws after burn-in (see below). Will be automatically coerced to integer.
#' The default value is 10000.
#' @param burnin single number greater or equal to 0, indicating the number of
#' draws discarded as burn-in. Will be automatically coerced to integer. The
#' default value is 3000.
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
#' \code{(phi+1)/2}, where \code{phi} denotes the persistence of the
#' log-volatility. The default value is \code{c(5, 1.5)}, which constitutes a
#' prior that puts some belief in a persistent log-volatility but also
#' encompasses the region where \code{phi} is around 0.
#' @param priorsigma single positive real number, which stands for the scaling
#' of the transformed parameter \code{sigma^2}, where \code{sigma} denotes the
#' volatility of log-volatility. More precisely, \code{sigma^2 ~ priorsigma *
#' chisq(df = 1)}. The default value is \code{1}, which constitutes a
#' reasonably vague prior for many common exchange rate datasets, stock returns
#' and the like.
#' @param priorrho numeric vector of length 2, indicating the shape parameters
#' for the Beta prior distribution of the transformed parameter
#' \code{(rho+1)/2}, where \code{rho} denotes the conditional correlation
#' between observation and the increment of the
#' log-volatility. The default value is \code{c(4, 4)}, which constitutes a
#' slightly informative prior around 0 (the no leverage case) to boost convergence.
#' @param priorbeta numeric vector of length 2, indicating the mean and
#' standard deviation of the Gaussian prior for the regression parameters. The
#' default value is \code{c(0, 10000)}, which constitutes a very vague prior
#' for many common datasets. Not used if \code{designmatrix} is \code{NA}.
#' @param thinpara single number greater or equal to 1, coercible to integer.
#' Every \code{thinpara}th parameter draw is kept and returned. The default
#' value is 1, corresponding to no thinning of the parameter draws i.e. every
#' draw is stored.
#' @param thinlatent single number greater or equal to 1, coercible to integer.
#' Every \code{thinlatent}th latent variable draw is kept and returned. The
#' default value is 1, corresponding to no thinning of the latent variable
#' draws, i.e. every draw is kept.
#' @param thintime \emph{deprecated}
#' @param quiet logical value indicating whether the progress bar and other
#' informative output during sampling should be omitted. The default value is
#' \code{FALSE}, implying verbose output.
#' @param startpara \emph{optional} named list, containing the starting values
#' for the parameter draws. If supplied, \code{startpara} must contain four
#' elements named \code{mu}, \code{phi}, \code{sigma}, and \code{rho}, where \code{mu} is
#' an arbitrary numerical value, \code{phi} is a real number between \code{-1}
#' and \code{1}, \code{sigma} is a positive real number, and \code{rho} is
#' a real number between \code{-1} and \code{1}. The default value is
#' is \code{list(mu = -10, phi = 0.9, sigma = 0.3, rho = 0)}.
#' @param startlatent \emph{optional} vector of length \code{length(y)},
#' containing the starting values for the latent log-volatility draws. The
#' default value is \code{rep(-10, length(y))}.
#' @param expert \emph{optional} named list of expert parameters. For most
#' applications, the default values probably work best. If
#' \code{expert} is provided, it may contain the following named elements:
#' 
#' \code{parameterization}: Character string containing values \code{"centered"},
#' and \code{"noncentered"}. Alternatively, it can be a single element character
#' vector of the form \code{"asisX"}, where \code{X} is an integer, which is
#' equivalent to \code{rep_len(c("centered", "noncentered"), X)}.
#' Defaults to \code{"asis5"}.
#' 
#' \code{mhcontrol}: Single numeric value controlling the proposal density of a
#' Metropolis-Hastings (MH) update step when jointly sampling \code{mu}, \code{phi},
#' \code{sigma}, and \code{rho}. It controls the stepsize of a log-random-walk
#' proposal. Defaults to \code{0.1}.
#' 
#' \code{gammaprior}: Single logical value indicating whether a Gamma prior for
#' \code{sigma^2} should be used. If set to \code{FALSE}, a moment-matched Inverse Gamma
#' prior is employed. Defaults to \code{TRUE}.
#' 
#' \code{init.with.svsample}: Single integer indicating the length of a ``pre-burnin'' run using
#' the computationally much more efficient \code{\link{svsample}}. This run helps
#' in finding good initial values for the latent states, giving \code{svlsample}
#' a considerable initial boost for convergence. Defaults to \code{100L}.
#' @param \dots Any extra arguments will be forwarded to
#' \code{\link{updatesummary}}, controlling the type of statistics calculated
#' for the posterior draws.
#' @return The value returned is a list object of class \code{svldraws} holding
#' \item{para}{\code{mcmc} object containing the \emph{parameter} draws from
#' the posterior distribution.}
#' \item{latent}{\code{mcmc} object containing the
#' \emph{latent instantaneous log-volatility} draws from the posterior
#' distribution.}
#' \item{beta}{\code{mcmc} object containing the \emph{regression coefficient}
#' draws from the posterior distribution \emph{(optional)}.}
#' \item{y}{the argument \code{y}.}
#' \item{runtime}{\code{proc_time} object containing the
#' run time of the sampler.}
#' \item{priors}{\code{list} containing the parameter
#' values of the prior distribution, i.e. the arguments \code{priormu},
#' \code{priorphi}, \code{priorsigma}, and \code{priorrho}, and potentially
#' \code{priorbeta}.}
#' \item{thinning}{\code{list} containing the thinning
#' parameters, i.e. the arguments \code{thinpara}, \code{thinlatent} and
#' \code{thintime}.}
#' \item{summary}{\code{list} containing a collection of
#' summary statistics of the posterior draws for \code{para}, and \code{latent}.}
#' \item{meanmodel}{\code{character} containing information about how \code{designmatrix}
#' was used.}
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
#' @author Darjus Hosszejni \email{darjus.hosszejni@@wu.ac.at}
#' @seealso \code{\link{svsim}}, \code{\link{svsample}}, \code{\link{updatesummary}},
#' \code{\link{predict.svdraws}}, \code{\link{plot.svdraws}}.
#' @keywords models ts
#' @examples
#' #\dontrun{
#' # Example 1
#' ## Simulate a short SVL process
#' sim <- svsim(200, mu = -10, phi = 0.95, sigma = 0.2, rho = -0.4)
#' 
#' ## Obtain 5000 draws from the sampler (that's not a lot)
#' draws <- svlsample(sim$y, draws = 5000, burnin = 2000,
#' 		  priormu = c(-10, 1), priorphi = c(20, 1.5), priorsigma = 0.05,
#'      priorrho = c(3, 9))
#' 
#' ## Check out the results
#' summary(draws)
#' plot(draws, simobj = sim)
#' 
#' 
#' # Example 2
#' ## AR(1) structure for the mean
#' data(exrates)
#' len <- 1200
#' ahead <- 100
#' y <- head(exrates$USD, len)
#' 
#' ## Fit AR(1)-SVL model to EUR-USD exchange rates
#' res <- svlsample(y, designmatrix = "ar1")
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
#' # Example 3
#' ## Predicting USD based on JPY and GBP in the mean
#' data(exrates)
#' len <- 1200
#' ahead <- 30
#' ## Calculate log-returns
#' logreturns <- apply(exrates[, c("USD", "JPY", "GBP")], 2,
#'                     function (x) diff(log(x)))
#' logretUSD <- logreturns[2:(len+1), "USD"]
#' regressors <- cbind(1, as.matrix(logreturns[1:len, ]))  # lagged by 1 day
#' 
#' ## Fit SV model to EUR-USD exchange rates
#' res <- svlsample(logretUSD, designmatrix = regressors)
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
#' #}
#' @export
svlsample <- function (y, draws = 10000, burnin = 3000, designmatrix = NA,
                       priormu = c(0, 100), priorphi = c(5, 1.5), priorsigma = 1,
                       priorrho = c(4, 4), priorbeta = c(0, 10000),
                       thinpara = 1, thinlatent = 1, thintime = 1,
                       quiet = FALSE, startpara, startlatent, expert, ...) {
  # Some error checking for y
  if (inherits(y, "svsim")) {
    y <- y[["y"]]
    warning("Extracted data vector from 'svsim'-object.")
  }
  if (!is.numeric(y)) stop("Argument 'y' (data vector) must be numeric.")
  if (length(y) < 2) stop("Argument 'y' (data vector) must contain at least two elements.")
  
  myoffset <- if (any(is.na(designmatrix)) && any(y^2 == 0)) 1.0e-8 else 0
  if (myoffset > 0) {
    warning(paste("Argument 'y' (data vector) contains zeros. I am adding an offset constant of size ", myoffset, " to do the auxiliary mixture sampling. If you want to avoid this, you might consider de-meaning the returns before calling this function.", sep=""))
  }

  # Some error checking for draws
  if (!is.numeric(draws) || length(draws) != 1 || draws < 1) {
    stop("Argument 'draws' (number of MCMC iterations after burn-in) must be a single number >= 1.")
  } else {
    draws <- as.integer(draws)
  }

  # Some error checking for burnin
  if (!is.numeric(burnin) || length(burnin) != 1 || burnin < 0) {
    stop("Argument 'burnin' (burn-in period) must be a single number >= 0.")
  } else {
    burnin <- as.integer(burnin)
  }

  # Some error checking for designmatrix
  meanmodel <- "matrix"
  if (any(is.na(designmatrix))) {
    designmatrix <- matrix(NA)
    meanmodel <- "none"
  } else {
    if (any(grep("ar[0-9]+$", as.character(designmatrix)[1]))) {
      arorder <- as.integer(gsub("ar", "", as.character(designmatrix)))
      meanmodel <- paste0("ar", arorder)
      if (length(y) <= arorder + 1) stop("Time series 'y' is to short for this AR process.")
      designmatrix <- matrix(rep(1, length(y) - arorder), ncol = 1)
      colnames(designmatrix) <- c("const")
      if (arorder >= 1) {
        for (i in 1:arorder) {
          oldnames <- colnames(designmatrix)
          designmatrix <- cbind(designmatrix, y[(arorder-i+1):(length(y)-i)])
          colnames(designmatrix) <- c(oldnames, paste0("ar", i))
        }
        y <- y[-(1:arorder)]
      }
    }
    if (!is.numeric(designmatrix)) stop("Argument 'designmatrix' must be a numeric matrix or an AR-specification.")
    if (!is.matrix(designmatrix)) {
      designmatrix <- matrix(designmatrix, ncol = 1)
    }
    if (!nrow(designmatrix) == length(y)) stop("Number of columns of argument 'designmatrix' must be equal to length(y).")
  }

  # Some error checking for the prior parameters 
  if (!is.numeric(priormu) || length(priormu) != 2) {
    stop("Argument 'priormu' (mean and variance for the Gaussian prior for mu) must be numeric and of length 2.")
  }

  if (!is.numeric(priorphi) || length(priorphi) != 2) {
    stop("Argument 'priorphi' (shape1 and shape2 parameters for the Beta prior for (phi+1)/2) must be numeric and of length 2.")
  }

  if (!is.numeric(priorsigma) || length(priorsigma) != 1 || priorsigma <= 0) {
    stop("Argument 'priorsigma' (scaling of the chi-squared(df = 1) prior for sigma^2) must be a single number > 0.")
  }

  if (!is.numeric(priorrho) || length(priorrho) != 2) {
    stop("Argument 'priorrho' (shape1 and shape2 parameters for the Beta prior for (rho+1)/2) must be numeric and of length 2.")
  }

  if (!is.numeric(priorbeta) || length(priorbeta) != 2) {
    stop("Argument 'priorbeta' (means and sds for the independent Gaussian priors for beta) must be numeric and of length 2.")
  }

  # Some error checking for thinpara
  if (!is.numeric(thinpara) || length(thinpara) != 1 || thinpara < 1) {
    stop("Argument 'thinpara' (thinning parameter for mu, phi, sigma, and rho) must be a single number >= 1.")
  } else {
    thinpara <- as.integer(thinpara)
  }

  # Some error checking for thinlatent
  if (!is.numeric(thinlatent) || length(thinlatent) != 1 || thinlatent < 1) {
    stop("Argument 'thinlatent' (thinning parameter for the latent log-volatilities) must be a single number >= 1.")
  } else {
    thinlatent <- as.integer(thinlatent)
  }

  # thintime deprecated
  if (!(identical(thintime, 1L) || identical(thintime, 1.0))) {
    thintime <- 1L
    warning("Parameter 'thintime' is deprecated. Setting 'thintime' to 1.")
  }

  # Some input checking for startpara
  startparadefault <- list(mu = priormu[1],
                           phi = 2 * (priorphi[1] / sum(priorphi)) - 1,
                           sigma = priorsigma,
                           rho = 2 * (priorrho[1] / sum(priorrho)) - 1)
  if (missing(startpara)) {
    startpara <- startparadefault
  } else {
    if (!is.list(startpara))
      stop("Argument 'startpara' must be a list. Its elements must be named 'mu', 'phi', 'sigma', and 'rho'.")

    if (!is.numeric(startpara[["mu"]]))
      stop('Argument \'startpara[["mu"]]\' must exist and be numeric.')

    if (!is.numeric(startpara[["phi"]]))
      stop('Argument \'startpara[["phi"]]\' must exist and be numeric.')

    if (abs(startpara[["phi"]]) >= 1)
      stop('Argument \'startpara[["phi"]]\' must be between -1 and 1.')

    if (!is.numeric(startpara[["sigma"]]))
      stop('Argument \'startpara[["sigma"]]\' must exist and be numeric.')

    if (startpara[["sigma"]] <= 0)
      stop('Argument \'startpara[["sigma"]]\' must be positive.')

    if (!is.numeric(startpara[["rho"]]))
      stop('Argument \'startpara[["rho"]]\' must exist and be numeric.')

    if (abs(startpara[["rho"]]) >= 1)
      stop('Argument \'startpara[["rho"]]\' must be between -1 and 1.')
  }

  # Some input checking for startlatent
  if (missing(startlatent)) {
    startlatent <- rep(-10, length(y))
  } else {
    if (!is.numeric(startlatent) | length(startlatent) != length(y))
      stop("Argument 'startlatent' must be numeric and of the same length as the data 'y'.")
  }

  # Some error checking for expert
  strategies <- c("centered", "noncentered")
  expertdefault <- list(parameterization = rep(strategies, 5),  # default: ASISx5
                        mhcontrol = 0.1,
                        gammaprior = TRUE,
                        init.with.svsample = 100L)
  if (missing(expert)) {
    parameterization <- expertdefault$parameterization
    mhcontrol <- expertdefault$mhcontrol
    gammaprior <- expertdefault$gammaprior
    init.with.svsample <- expertdefault$init.with.svsample
  } else {
    expertnames <- names(expert)
    if (!is.list(expert) || is.null(expertnames) || any(expertnames == ""))
      stop("Argument 'expert' must be a named list with nonempty names.")
    if (length(unique(expertnames)) != length(expertnames))
      stop("No duplicate elements allowed in argument 'expert'.")
    allowednames <- c("parameterization", "mhcontrol", "gammaprior", "init.with.svsample")
    exist <- pmatch(expertnames, allowednames)
    if (any(is.na(exist)))
      stop(paste("Illegal element '", paste(expertnames[is.na(exist)], collapse="' and '"), "' in argument 'expert'.", sep=''))

    expertenv <- list2env(expert) 

    if (exists("parameterization", expertenv)) {
      if (!is.character(expert[["parameterization"]]))
        stop("Argument 'parameterization' must be either a vector of 'centered', 'noncentered' values or a character string of form 'asis#' with # a positive integer.")
      nmatches <- grep("^asis[1-9][0-9]*$", expert[["parameterization"]])
      if (length(nmatches) == 0) {
        parameterization <- match.arg(expert[["parameterization"]], strategies, several.ok = TRUE)
      } else if (length(nmatches) == 1) {
        parameterization <- rep(strategies, nmatches)
      } else {
        parameterization <- NA
      }
      if (!all(parameterization %in% strategies)) {
        stop("Argument 'parameterization' must be either a vector of 'centered', 'noncentered' values or a character string of form 'asis#' with # a positive integer.")
      }
    } else {
      parameterization <- expertdefault$parameterization
    }

    # Remark: mhcontrol > 0 controls stepsize of log-random-walk proposal
    if (exists("mhcontrol", expertenv)) {
      mhcontrol <- expert[["mhcontrol"]]
      if (!is.numeric(mhcontrol) || length(mhcontrol) != 1 || mhcontrol <= 0)
        stop("Argument 'mhcontrol' must be a single positive number.")
    } else {
      mhcontrol <- expertdefault$mhcontrol
    }

    # use a Gamma prior for sigma^2 in C?
    if (exists("gammaprior", expertenv)) {
      gammaprior <- expert[["gammaprior"]]
      if (!is.logical(gammaprior)) stop("Argument 'gammaprior' must be TRUE or FALSE.")
    } else {
      gammaprior <- expertdefault$gammaprior
    }

    if (exists("init.with.svsample", expertenv)) {
      init.with.svsample <- expert[["init.with.svsample"]]
      if (!is.numeric(init.with.svsample) || length(init.with.svsample) != 1 || init.with.svsample < 0)
        stop("Argument 'init.with.svsample' must be a single non-negative number.")
      init.with.svsample <- as.integer(init.with.svsample)
    } else {
      init.with.svsample <- expertdefault$init.with.svsample
    }
  }
  
  if (init.with.svsample > 0L) {
    if (!quiet) {
      cat(paste0("\nInitial values: calling function 'svsample' with ", init.with.svsample, " iter. Series length is ", length(y), ".\n"), file=stderr())
      flush.console()
    }

    init.runtime <- system.time({
      init.res <- svsample(y,  # not svsample2 because of 'designmatrix'
                           draws = as.integer(init.with.svsample/2),
                           burnin = as.integer(init.with.svsample/2),
                           designmatrix = designmatrix, priormu = priormu,
                           priorphi = priorphi, priorsigma = priorsigma,
                           priornu = NA, priorbeta = priorbeta,
                           priorlatent0 = "stationary",
                           thinpara = 1L, thinlatent = 1L, thintime = 1L,
                           keeptau = FALSE, quiet = TRUE,
                           startpara = startpara[c("phi", "mu", "sigma")],
                           startlatent = startlatent,
                           expert = list(gammaprior = gammaprior,
                                         parameterization = "GIS_C"),
                           quantiles = .5, esspara = FALSE,  # updatesummary
                           esslatent = FALSE)
    })

    if (!quiet) {
      cat(paste0("\nInitial values: time taken by 'svsample': ", round(init.runtime["elapsed"], 3), " seconds.\n"), file=stderr())
    }

    startpara[c("mu", "phi", "sigma")] <-
      as.list(apply(init.res$para, 2, median))
    startlatent <- as.numeric(apply(init.res$latent, 2, median))
  }

  renameparam <- c("centered" = "C", "noncentered" = "NC")
  if (!quiet) {
    cat(paste("\nCalling (", paste(renameparam[parameterization], collapse=", "), ") MCMC sampler with ", draws+burnin, " iter. Series length is ", length(y), ".\n",sep=""), file=stderr())
    flush.console()
  }

  if (.Platform$OS.type != "unix") myquiet <- TRUE else myquiet <- quiet  # Hack to prevent console flushing problems with Windows

  runtime <- system.time({
    res <- svlsample_cpp(y, draws, burnin, designmatrix, thinpara, thinlatent, thintime,
                                 startpara, startlatent,
                                 priorphi[1], priorphi[2], priorrho[1], priorrho[2],
                                 0.5, 0.5/priorsigma, priormu[1], priormu[2],
                                 priorbeta[1], priorbeta[2], !myquiet,
                                 myoffset, mhcontrol, gammaprior, parameterization)
  })

  if (any(is.na(res))) stop("Sampler returned NA. This is most likely due to bad input checks and shouldn't happen. Please report to package maintainer.")

  if (!quiet) {
    cat("Timing (elapsed): ", file=stderr())
    cat(runtime["elapsed"], file=stderr())
    cat(" seconds.\n", file=stderr())
    cat(round((draws+burnin)/runtime[3]), "iterations per second.\n\n", file=stderr())
    cat("Converting results to coda objects... ", file=stderr())
  }
  
  colnames(res$para) <- c("mu", "phi", "sigma", "rho")
  colnames(res$latent) <- paste0('h_', seq(1, length(y), by=thintime))
  # create svldraws class
  res$runtime <- runtime
  res$y <- y
  res$para <- coda::mcmc(res$para, burnin+thinpara, burnin+draws, thinpara)
  res$latent <- coda::mcmc(res$latent, burnin+thinlatent, burnin+draws, thinlatent)
  res$thinning <- list(para = thinpara, latent = thinlatent, time = thintime)
  res$priors <- list(mu = priormu, phi = priorphi, sigma = priorsigma, rho = priorrho, gammaprior = gammaprior)
  if (!any(is.na(designmatrix))) {
    res$beta <- coda::mcmc(res$beta, burnin+thinpara, burnin+draws, thinpara)
    colnames(res$beta) <- paste0("b_", 0:(NCOL(designmatrix)-1))
    res$priors <- c(res$priors, "beta" = list(priorbeta), "designmatrix" = list(designmatrix))
  } else {
    res$beta <- NULL
  }
  res$meanmodel <- meanmodel
  class(res) <- c("svldraws", "svdraws")

  if (!quiet) {
    cat("Done!\n", file=stderr())
    cat("Summarizing posterior draws... ", file=stderr())
  }
  res <- updatesummary(res, ...)

  if (!quiet) cat("Done!\n\n", file=stderr())
  res
}

#' Minimal overhead version of \code{\link{svlsample}}.
#'
#' \code{svsample2} is a minimal overhead version of \code{\link{svsample}}
#' with slightly different default arguments and a simplified return value
#' structure. It is intended to be used mainly for one-step updates where speed
#' is an issue, e.g., as a plug-in into other MCMC samplers. Note that
#' absolutely no input checking is performed, thus this function is to be used
#' with proper care!
#'
#' @export
svlsample2 <- function (y, draws = 1, burnin = 0,
                       priormu = c(0, 100), priorphi = c(5, 1.5), priorsigma = 1, priorrho = c(3, 5),
                       thinpara = 1, thinlatent = 1, thintime = 1,
                       quiet = FALSE, startpara, startlatent) {

  # thintime deprecated
  if (!(identical(thintime, 1L) || identical(thintime, 1.0))) {
    thintime <- 1L
    warning("Parameter 'thintime' is deprecated. Setting 'thintime' to 1.")
  }

  res <- svlsample_cpp(y, draws, burnin, matrix(NA), thinpara, thinlatent, thintime,
                               startpara, startlatent,
                               priorphi[1], priorphi[2], priorrho[1], priorrho[2],
                               0.5, 0.5/priorsigma, priormu[1], priormu[2],
                               0, 1, !quiet,
                               0, 0.1, TRUE, rep(c("centered", "noncentered"), 5))

  res$para <- t(res$para)
  res$latent <- t(res$latent)
  res$meanmodel <- "none"
  rownames(res$para) <- c("mu", "phi", "sigma", "rho")
  res
}


.svsample <- function (...) {
  .Deprecated("svsample2")
  svsample(...)
}
