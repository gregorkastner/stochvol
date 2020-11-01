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

#' Markov Chain Monte Carlo (MCMC) Sampling for the Stochastic Volatility (SV)
#' Model
#' 
#' \code{svsample} simulates from the joint posterior distribution of the SV
#' parameters \code{mu}, \code{phi}, \code{sigma} (and potentially \code{nu} and \code{rho}),
#' along with the latent log-volatilities \code{h_0,...,h_n} and returns the
#' MCMC draws. If a design matrix is provided, simple Bayesian regression can
#' also be conducted.
#' 
#' Functions \code{svtsample}, \code{svlsample}, and \code{svtlsample} are
#' wrappers around \code{svsample} with convenient default values for the SV
#' model with t-errors, leverage, and both t-errors and leverage, respectively.
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
#' @param thin single number greater or equal to 1, coercible to integer.
#' Every \code{thinpara}th parameter and latent draw is kept and returned. The default
#' value is 1, corresponding to no thinning of the parameter draws i.e. every
#' draw is stored.
#' @param thinpara single number greater or equal to 1, coercible to integer.
#' Every \code{thinpara}th parameter draw is kept and returned. The default
#' value is \code{thin}.
#' @param thinlatent single number greater or equal to 1, coercible to integer.
#' Every \code{thinlatent}th latent variable draw is kept and returned. The
#' default value is \code{thin}
#' @param keeptime Either 'all' (the default) or 'last'. Indicates which latent
#' volatility draws should be stored.
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
#' In case of parallel execution with \code{cl} provided, \code{startpara} can be a list of
#' named lists that initialize the parallel chains.
#' @param startlatent \emph{optional} vector of length \code{length(y)},
#' containing the starting values for the latent log-volatility draws. The
#' default value is \code{rep(-10, length(y))}.
#' In case of parallel execution with \code{cl} provided, \code{startlatent} can be a list of
#' named lists that initialize the parallel chains.
#' @param parallel \emph{optional} one of \code{"no"} (default), \code{"multicore"}, or \code{"snow"},
#' indicating what type of parallellism is to be applied. Option
#' \code{"multicore"} is not available on Windows.
#' @param n_cpus \emph{optional} positive integer, the number of CPUs to be used in case of
#' parallel computations. Defaults to \code{1L}. Ignored if parameter
#' \code{cl} is supplied and \code{parallel != "snow"}.
#' @param cl \emph{optional} so-called SNOW cluster object as implemented in package
#' \code{parallel}. Ignored unless \code{parallel == "snow"}.
#' @param n_chains \emph{optional} positive integer specifying the number of independent MCMC chains
#' @param print_progress \emph{optional} one of \code{"automatic"}, \code{"progressbar"},
#' or \code{"iteration"}, controls the output. Ignored if \code{quiet} is \code{TRUE}.
#' @param expert \emph{optional} named list of expert parameters. For most
#' applications, the default values probably work best. Interested users are
#' referred to the literature provided in the References section. If
#' \code{expert} is provided, it may contain the following named elements:
#' \itemize{
#' \item{interweave}{Logical value. If \code{TRUE} (the default),
#' then ancillarity-sufficiency interweaving strategy (ASIS) is applied
#' to improve on the sampling efficiency for the parameters.
#' Otherwise one parameterization is used.}
#' \item{correct_model_misspecification}{Logical value. If \code{FALSE}
#' (the default), then auxiliary mixture sampling is used to sample the latent
#' states. If \code{TRUE}, extra computations are made to correct for model
#' misspecification either ex-post by reweighting or on-line using a
#' Metropolis-Hastings step.}
#' }
#' @param \dots Any extra arguments will be forwarded to
#' \code{\link{updatesummary}}, controlling the type of statistics calculated
#' for the posterior draws.
#' @return The value returned is a list object of class \code{svdraws} holding
#' \item{para}{\code{mcmc.list} object containing the \emph{parameter} draws from
#' the posterior distribution.}
#' \item{latent}{\code{mcmc.list} object containing the
#' \emph{latent instantaneous log-volatility} draws from the posterior
#' distribution.}
#' \item{latent0}{\code{mcmc.list} object containing the \emph{latent
#' initial log-volatility} draws from the posterior distribution.}
#' \item{tau}{\code{mcmc.list} object containing the \emph{latent variance inflation
#' factors} for the sampler with conditional t-innovations \emph{(optional)}.}
#' \item{beta}{\code{mcmc.list} object containing the \emph{regression coefficient}
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
#' @example inst/examples/svsample.R
#' @export
svsample <- function(y, draws = 10000, burnin = 1000, designmatrix = NA,
                     priormu = c(0, 100), priorphi = c(5, 1.5), priorsigma = 1,
                     priornu = 0, priorrho = NA,
                     priorbeta = c(0, 10000), priorlatent0 = "stationary",
                     priorspec = NULL, thin = 1,
                     thinpara = thin, thinlatent = thin, keeptime = "all",
                     quiet = FALSE, startpara = NULL, startlatent = NULL,
                     parallel = c("no", "multicore", "snow"),
                     n_cpus = 1L, cl = NULL, n_chains = 1L,
                     print_progress = "automatic",
                     expert = NULL, ...) {

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
                     phi = sv_beta(shape1 = priorphi[1], shape2 = priorphi[2]),
                     sigma2 = sv_gamma(shape = 0.5, rate = 0.5 / priorsigma),
                     nu = if (priornu == 0) sv_infinity() else sv_exponential(rate = priornu),
                     rho = if (isTRUE(is.na(priorrho))) sv_constant(value = 0) else sv_beta(shape1 = priorrho[1], shape2 = priorrho[2]),
                     beta = sv_multinormal(mean = priorbeta[1], sd = priorbeta[2], dim = NCOL(designmatrix)),
                     latent0_variance = priorlatent0)
  } else if (!inherits(priorspec, "sv_priorspec")) {
    stop("Received argument 'priorspec' but it does not have the correct form. Please refer to the function called 'specify_priors'.")
  }
  keeptau <- !isTRUE(inherits(priorspec$nu, "sv_infinity"))

  ## thinning parameters
  validate_thinning(thinpara, thinlatent, keeptime)
  thinpara <- as.integer(thinpara)
  thinlatent <- as.integer(thinlatent)
  thintime <- switch(keeptime,
                     all = 1L,
                     last = length(y))

  ## parallel (strongly influenced by package 'boot')
  parallel <- match.arg(parallel)
  have_mc <- FALSE
  have_snow <- FALSE
  if (parallel != "no") {
    if (parallel == "multicore") {
      have_mc <- .Platform$OS.type != "windows"
    } else if (parallel == "snow") {
      have_snow <- TRUE
    }
    if (!have_mc && !have_snow) {
      n_cpus <- 1L
    }
    requireNamespace("parallel")
  }
  create_cluster <- isTRUE(is.null(cl))
  have_cluster <- isTRUE(inherits(cl, "SOCKcluster"))
  if (have_snow && !create_cluster && !have_cluster) {
    warning("Unknown type of input in parameter 'cl'. Should be either NULL or of class 'cluster' created by the parallel package. Turning off parallelism")
    cl <- NULL
    have_mc <- FALSE
    have_snow <- FALSE
    n_cpus <- 1L
  }
  if (have_snow && n_cpus == 1L && !have_cluster) {
    warning("Inconsistent settings for parallel execution: snow parallelism combined with n_cpus == 1 and no supplied cluster object. Turning off parallelism")
    have_snow <- FALSE
  }
  if (have_mc && n_cpus == 1L) {
    warning("Inconsistent settings for parallel execution: multicore parallelism combined with n_cpus == 1. Turning off parallelism")
    have_mc <- FALSE
  }
  n_workers <- if (have_snow) {
    if (have_cluster) {
      length(cl)
    } else {
      n_cpus
    }
  } else {
    1L
  }
  assert_single(n_chains, "Parameter for number of chains")
  assert_numeric(n_chains, "Parameter for number of chains")
  n_chains <- as.integer(n_chains)
  assert_positive(n_chains, "Parameter for number of chains")

  ## expert
  expert <- validate_and_process_expert(expert)
  correct_model_misspecification <- expert$correct_model_misspecification
  interweave <- expert$interweave
  fast_sv <- expert$fast_sv
  general_sv <- expert$general_sv

  # Initial values
  have_parallel_startpara <- !is.null(startpara) && isTRUE(all(sapply(startpara, is.list)))
  have_parallel_startlatent <- !is.null(startlatent) && isTRUE(is.list(startlatent))
  if (have_parallel_startpara && length(startpara) != n_chains) {
    stop("Got list of lists for parameter 'startpara' of length ", length(startpara), " but the number of chains is ", n_chains, ". The two numbers should match.")
  }
  if (have_parallel_startlatent && length(startlatent) != n_chains) {
    stop("Got list for parameter 'startlatent' of length ", length(startlatent), " but the number of chains is ", n_chains, ". The two numbers should match.")
  }
  startbetadefault <- if (meanmodel == "none") {
    mean(priorspec$beta)
  } else {
    init_beta(y, designmatrix)
  }
  startmudefault <- if (meanmodel == "none") {
    init_mu(y, priorspec)
  } else {
    init_mu(y, priorspec, X = designmatrix, beta_hat = startbetadefault)
  }
  startparadefault <-
    list(mu = startmudefault,
         phi = if (inherits(priorspec$phi, "sv_beta")) 2 * mean(priorspec$phi) - 1 else mean(priorspec$phi),
         sigma = sqrt(mean(priorspec$sigma2)),
         nu = 2 + mean(priorspec$nu),
         rho = if (inherits(priorspec$rho, "sv_beta")) 2 * mean(priorspec$rho) - 1 else mean(priorspec$rho),
         beta = startbetadefault,
         latent0 = -10)
  startpara <- if (have_parallel_startpara) {
    lapply(startpara, function (x, def) apply_default_list(x, def), def = startparadefault)
  } else {
    replicate(n_chains, apply_default_list(startpara, startparadefault), simplify = FALSE)
  }
  startlatent <- if (have_parallel_startlatent) {
    lapply(startlatent, function (x, def) apply_default_list(x, def), def = rep.int(NA, length(y)))
  } else if (is.null(startlatent)) {
    lapply(startpara, function (x, len) rep.int(x$mu, len), len = length(y))  # default startlatent is constant mu
  } else {
    replicate(n_chains, apply_default_list(startlatent, rep.int(NA, length(y))), simplify = FALSE)
  }
  validate_initial_values(startpara, startlatent, y, designmatrix)

  # Decision about the sampler
  use_fast_sv <-
    # rho == 0
    (inherits(priorspec$rho, "sv_constant") && isTRUE(priorspec$rho$value == 0)) &&
    # mu is either 0 or normal
    ((inherits(priorspec$mu, "sv_constant") && isTRUE(priorspec$mu$value == 0)) ||
       inherits(priorspec$mu, "sv_normal")) &&
    # fast SV can't correct for misspecification and also do t-errors and/or regression
    (!correct_model_misspecification || (inherits(priorspec$nu, "sv_infinity") && meanmodel == "none")) &&
    # prior for phi is beta
    inherits(priorspec$phi, "sv_beta") &&
    # prior for sigma is gamma(0.5, _)
    (inherits(priorspec$sigma2, "sv_gamma") && priorspec$sigma2$shape == 0.5)

  # Pick sampling function
  myquiet <- (.Platform$OS.type != "unix") || quiet  # Hack to prevent console flushing problems with Windows
  ## print progress
  print_settings <- if (is.character(print_progress)) {
    print_progress <- match.arg(print_progress, c("automatic", "progressbar", "iteration"))
    if (print_progress == "automatic") {
      n_chains
    } else if (print_progress == "progressbar") {
      1
    } else {
      max(c(2, n_chains))
    }
  } else {  # hidden feature
    print_progress
  }

  sampling_function <- if (use_fast_sv) {
    function (chain) {
      res <- stochvol::svsample_fast_cpp(y, draws, burnin, designmatrix, priorspec,
                                         thinpara, thinlatent, keeptime,
                                         startpara[[chain]], startlatent[[chain]], keeptau,
                                         list(quiet = myquiet, chain = chain, n_chains = print_settings),
                                         correct_model_misspecification, interweave, myoffset, fast_sv)
      if (correct_model_misspecification) {
        para_indices <- sample.int(n = NROW(res$para), size = NROW(res$para), replace = TRUE, prob = res$correction_weight_para, useHash = FALSE)
        res$para <- res$para[para_indices, , drop = FALSE]
        latent_indices <- if (thinpara == thinlatent) {  # same re-sampling if thinning is the same
          para_indices
        } else {  # separate re-sampling if thinning is different => WARNING! joint distribution of para and latent is gone!
          sample.int(n = NROW(res$latent), size = NROW(res$latent), replace = TRUE, prob = res$correction_weight_latent, useHash = FALSE)
        }
        res$latent <- res$latent[latent_indices, , drop = FALSE]
      }
      res
    }
  } else {
    function (chain) {
      stochvol::svsample_general_cpp(y, draws, burnin, designmatrix, priorspec,
                                     thinpara, thinlatent, keeptime,
                                     startpara[[chain]], startlatent[[chain]], keeptau,
                                     list(quiet = myquiet, chain = chain, n_chains = print_settings),
                                     correct_model_misspecification, interweave, myoffset, general_sv)
    }
  }

  # Print sampling info
  if (use_fast_sv) {
    para <- 1 + (fast_sv$baseline_parameterization == "noncentered") + 2 * interweave
    parameterization <- c("centered", "noncentered", "GIS_C", "GIS_NC")[para]

    if (!quiet) {
      cat(paste("\nCalling ", parameterization, " MCMC sampler with ", draws+burnin, " iter. Series length is ", length(y), ".\n",sep=""), file=stderr())
      flush.console()
    }
  } else {
    strategies <- if (interweave) c("centered", "noncentered") else expert$general_sv$starting_parameterization
    parameterization <- rep(strategies, general_sv$multi_asis)

    renameparam <- c("centered" = "C", "noncentered" = "NC")
    if (!quiet) {
      cat(paste("\nCalling ", asisprint(renameparam[parameterization], renameparam), " MCMC sampler with ", draws+burnin, " iter. Series length is ", length(y), ".\n",sep=""), file=stderr())
      flush.console()
    }
  }

  # Call sampler
  runtime <- system.time(reslist <-
    if ((n_cpus > 1L || have_cluster) && (have_mc || have_snow)) {
      if (have_mc) {
        parallel::mclapply(seq_len(n_chains), sampling_function, mc.cores = n_cpus)
      } else if (have_snow) {
        list(...) # evaluate any promises
        if (create_cluster) {
          cl <- parallel::makePSOCKcluster(rep("localhost", n_workers), outfile = NULL)
        }
        RNGkind(kind = "L'Ecuyer-CMRG")
        parallel::clusterSetRNGStream(cl)
        parallel::clusterEvalQ(cl, library(stochvol))
        parallel::clusterExport(cl, c("y", "draws", "burnin", "designmatrix", "priorspec",
                                      "thinpara", "thinlatent", "keeptime",
                                      "startpara", "startlatent", "keeptau",
                                      "myquiet", "n_chains", "print_settings",
                                      "correct_model_misspecification", "interweave", "myoffset",
                                      "fast_sv", "general_sv"),
                                envir = environment())
        if (create_cluster) {
          reslist <- tryCatch(parallel::parLapply(cl, seq_len(n_chains), sampling_function),
                              finally = parallel::stopCluster(cl))
          reslist
        } else {
          parallel::parLapply(cl, seq_len(n_chains), sampling_function)
        }
      }
    } else {
      lapply(seq_len(n_chains), sampling_function)
    }
  )
  res <- list()
  class(res) <- "svdraws"

  # Process results
  if (any(is.na(res))) {
    stop("Sampler returned NA. This is most likely due to bad input checks and shouldn't happen. Please report to package maintainer.")
  }

  if (!quiet) {
    cat("Timing (elapsed): ", file=stderr())
    cat(runtime["elapsed"], file=stderr())
    cat(" seconds.\n", file=stderr())
    cat(round((draws+burnin)*n_chains/runtime[3]), "iterations per second.\n\n", file=stderr())
    cat("Converting results to coda objects... ", file=stderr())
  }

  # store results:
  res$runtime <- runtime
  res$y <- y
  res$y_orig <- y_orig
  res$simobj <- simobj
  res$para <- coda::mcmc.list(lapply(reslist, function (x, d, b, th) coda::mcmc(x$para, b+th, b+d, th), d=draws, b=burnin, th=thinpara))
  res$latent <- coda::mcmc.list(lapply(reslist, function (x, d, b, th) coda::mcmc(x$latent, b+th, b+d, th), d=draws, b=burnin, th=thinlatent))
  res$latent0 <- coda::mcmc.list(lapply(reslist, function (x, d, b, th) coda::mcmc(x$latent0, b+th, b+d, th), d=draws, b=burnin, th=thinlatent))
  res$thinning <- list(para = thinpara, latent = thinlatent, time = keeptime)
  res$priors <- priorspec
  if (!any(is.na(designmatrix))) {
    res$beta <- coda::mcmc.list(lapply(reslist, function (x, d, b, th) coda::mcmc(x$beta, b+th, b+d, th), d=draws, b=burnin, th=thinpara))
    res$designmatrix <- designmatrix
  } else {
    res$beta <- NULL
  }
  if (keeptau) {
    res$tau <- coda::mcmc.list(lapply(reslist, function (x, d, b, th) coda::mcmc(x$tau, b+th, b+d, th), d=draws, b=burnin, th=thinlatent))
  } else {
    res$tau <- NULL
  }
  res$meanmodel <- meanmodel
  res$para_transform <- list(mu = function (x) {x},
                             phi = if (inherits(priorspec$phi, "sv_beta")) function (x) {(x+1)/2} else function (x) {x},
                             sigma = function (x) {x^2},
                             nu = if (inherits(priorspec$nu, "sv_exponential")) function (x) {x-2} else function (x) {x},
                             rho = if (inherits(priorspec$rho, "sv_beta")) function (x) {(x+1)/2} else function (x) {x})
  res$para_inv_transform <- list(mu = function (x) {x},
                                 phi = if (inherits(priorspec$phi, "sv_beta")) function (x) {2*x-1} else function (x) {x},
                                 sigma = function (x) {sqrt(x)},
                                 nu = if (inherits(priorspec$nu, "sv_exponential")) function (x) {x+2} else function (x) {x},
                                 rho = if (inherits(priorspec$rho, "sv_beta")) function (x) {2*x-1} else function (x) {x})
  res$para_transform_det <- list(mu = function (x) {1},
                                 phi = if (inherits(priorspec$phi, "sv_beta")) function (x) {.5} else function (x) {1},
                                 sigma = function (x) {2*x},
                                 nu = function (x) {1},
                                 rho = if (inherits(priorspec$rho, "sv_beta")) function (x) {.5} else function (x) {1})
  res$resampled <- use_fast_sv && correct_model_misspecification
  if (res$resampled) {
    res$correction_weight_para <- lapply(reslist, function (x) x$correction_weight_para)
    res$correction_weight_latent <- lapply(reslist, function (x) x$correction_weight_latent)
  }

  if (!quiet) {
    message("Done!")
    message("Summarizing posterior draws...")
  }
  res <- if (res$resampled) {
    message("No computation of effective sample size after re-sampling")
    updatesummary(res, esspara = FALSE, esslatent = FALSE, ...)
  } else {
    updatesummary(res, ...)
  }

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
       init_indicators = 5,
       init_tau = 1)
#' @rdname svsample_cpp
#' @export
default_general_sv <-
  list(multi_asis = 5,  # positive integer
       starting_parameterization = "centered",  # "centered" or "noncentered"
       update = list(latent_vector = TRUE, parameters = TRUE),
       init_tau = 1,
       proposal_diffusion_ken = FALSE)  # FALSE turns on adaptation

#' @rdname svsample
#' @export
svtsample <- function(y, draws = 10000, burnin = 1000, designmatrix = NA,
                      priormu = c(0, 100), priorphi = c(5, 1.5), priorsigma = 1,
                      priornu = 0.1, priorrho = NA,
                      priorbeta = c(0, 10000), priorlatent0 = "stationary",
                      priorspec = NULL, thin = 1,
                      thinpara = thin, thinlatent = thin, keeptime = "all",
                      quiet = FALSE, startpara = NULL, startlatent = NULL,
                      parallel = c("no", "multicore", "snow"),
                      n_cpus = 1L, cl = NULL, n_chains = 1L,
                      print_progress = "automatic",
                      expert = NULL, ...) {
  svsample(y = y, draws = draws, burnin = burnin, designmatrix = designmatrix,
           priormu = priormu, priorphi = priorphi, priorsigma = priorsigma,
           priornu = priornu, priorrho = priorrho,
           priorbeta = priorbeta, priorlatent0 = priorlatent0,
           priorspec = priorspec,
           thinpara = thinpara, thinlatent = thinlatent, keeptime = keeptime,
           quiet = quiet, startpara = startpara, startlatent = startlatent,
           parallel = parallel,
           n_cpus = n_cpus, cl = cl, n_chains = n_chains,
           print_progress = print_progress,
           expert = expert, ...)
}

#' @rdname svsample
#' @export
svlsample <- function(y, draws = 20000, burnin = 2000, designmatrix = NA,
                      priormu = c(0, 100), priorphi = c(5, 1.5), priorsigma = 1,
                      priornu = 0, priorrho = c(4, 4),
                      priorbeta = c(0, 10000), priorlatent0 = "stationary",
                      priorspec = NULL, thin = 1,
                      thinpara = thin, thinlatent = thin, keeptime = "all",
                      quiet = FALSE, startpara = NULL, startlatent = NULL,
                      parallel = c("no", "multicore", "snow"),
                      n_cpus = 1L, cl = NULL, n_chains = 1L,
                      print_progress = "automatic",
                      expert = NULL, ...) {
  svsample(y = y, draws = draws, burnin = burnin, designmatrix = designmatrix,
           priormu = priormu, priorphi = priorphi, priorsigma = priorsigma,
           priornu = priornu, priorrho = priorrho,
           priorbeta = priorbeta, priorlatent0 = priorlatent0,
           priorspec = priorspec,
           thinpara = thinpara, thinlatent = thinlatent, keeptime = keeptime,
           quiet = quiet, startpara = startpara, startlatent = startlatent,
           parallel = parallel,
           n_cpus = n_cpus, cl = cl, n_chains = n_chains,
           print_progress = print_progress,
           expert = expert, ...)
}

#' @rdname svsample
#' @export
svtlsample <- function(y, draws = 20000, burnin = 2000, designmatrix = NA,
                       priormu = c(0, 100), priorphi = c(5, 1.5), priorsigma = 1,
                       priornu = 0.1, priorrho = c(4, 4),
                       priorbeta = c(0, 10000), priorlatent0 = "stationary",
                       priorspec = NULL, thin = 1,
                       thinpara = thin, thinlatent = thin, keeptime = "all",
                       quiet = FALSE, startpara = NULL, startlatent = NULL,
                       parallel = c("no", "multicore", "snow"),
                       n_cpus = 1L, cl = NULL, n_chains = 1L,
                       print_progress = "automatic",
                       expert = NULL, ...) {
  svsample(y = y, draws = draws, burnin = burnin, designmatrix = designmatrix,
           priormu = priormu, priorphi = priorphi, priorsigma = priorsigma,
           priornu = priornu, priorrho = priorrho,
           priorbeta = priorbeta, priorlatent0 = priorlatent0,
           priorspec = priorspec,
           thinpara = thinpara, thinlatent = thinlatent, keeptime = keeptime,
           quiet = quiet, startpara = startpara, startlatent = startlatent,
           parallel = parallel,
           n_cpus = n_cpus, cl = cl, n_chains = n_chains,
           print_progress = print_progress,
           expert = expert, ...)
}

#' @rdname svsample
#' @export
svsample2 <- function(y, draws = 10000, burnin = 1000, designmatrix = NA,
                      priormu = c(0, 100), priorphi = c(5, 1.5), priorsigma = 1,
                      priornu = 0, priorrho = NA,
                      priorbeta = c(0, 10000), priorlatent0 = "stationary",
                      thinpara = 1, thinlatent = 1, keeptime = "all",
                      quiet = FALSE, startpara = NULL, startlatent = NULL) {
  .Deprecated("svsample_fast_cpp")
  svsample(y = y, draws = draws, burnin = burnin, designmatrix = designmatrix,
           priormu = priormu, priorphi = priorphi, priorsigma = priorsigma,
           priornu = priornu, priorrho = priorrho,
           priorbeta = priorbeta, priorlatent0 = priorlatent0,
           thinpara = thinpara, thinlatent = thinlatent, keeptime = keeptime,
           quiet = quiet, startpara = startpara, startlatent = startlatent)
}

#' Rolling Estimation of Stochastic Volatility Models
#' 
#' \code{svsample_roll} performs rolling window estimation based on \link{svsample}.
#' A convenience function for backtesting purposes.
#' 
#' Functions \code{svtsample_roll}, \code{svlsample_roll}, and \code{svtlsample_roll} are
#' wrappers around \code{svsample_roll} with convenient default values for the SV
#' model with t-errors, leverage, and both t-errors and leverage, respectively.
#' 
#' @param y numeric vector containing the data (usually log-returns), which
#' must not contain zeros. Alternatively, \code{y} can be an \code{svsim}
#' object. In this case, the returns will be extracted and a message is signalled.
#' @param designmatrix regression design matrix for modeling the mean. Must
#' have \code{length(y)} rows. Alternatively, \code{designmatrix} may be a
#' string of the form \code{"arX"}, where \code{X} is a nonnegative integer. To
#' fit a constant mean model, use \code{designmatrix = "ar0"} (which is
#' equivalent to \code{designmatrix = matrix(1, nrow = length(y))}). To fit an
#' AR(1) model, use \code{designmatrix = "ar1"}, and so on. If some elements of
#' \code{designmatrix} are \code{NA}, the mean is fixed to zero (pre-1.2.0
#' behavior of \pkg{stochvol}).
#' @param n_ahead number of time steps to predict from each time window.
#' @param forecast_length the time horizon at the end of the data set
#' that is used for backtesting.
#' @param n_start \emph{optional} the starting time point for backtesting.
#' Computed from \code{forecast_length} if omitted.
#' @param refit_every the SV model is refit every \code{refit_every} time steps.
#' Only the value \code{1} is allowed.
#' @param refit_window one of \code{"moving"} or \code{"expanding"}. If
#' \code{"expanding"}, then the start of the time window stays
#' at the beginning of the data set. If \code{"moving"}, then the
#' length of the time window is constant throughout backtesting.
#' @param calculate_quantile vector of numbers between 0 and 1.
#' These quantiles are predicted using \code{\link{predict.svdraws}}
#' for each time window.
#' @param calculate_predictive_likelihood boolean. If \code{TRUE},
#' the \code{n_ahead} predictive density is evaluated at the 
#' \code{n_ahead} time observation after each time window.
#' @param keep_draws boolean. If \code{TRUE}, the \code{svdraws} and
#' the \code{svpredict} objects are kept from each time window.
#' @param parallel one of \code{"no"} (default), \code{"multicore"}, or \code{"snow"},
#' indicating what type of parallellism is to be applied. Option
#' \code{"multicore"} is not available on Windows.
#' @param n_cpus \emph{optional} positive integer, the number of CPUs to be used in case of
#' parallel computations. Defaults to \code{1L}. Ignored if parameter
#' \code{cl} is supplied and \code{parallel != "snow"}.
#' @param cl \emph{optional} so-called SNOW cluster object as implemented in package
#' \code{parallel}. Ignored unless \code{parallel == "snow"}.
#' @param \dots Any extra arguments will be forwarded to
#' \code{\link{svsample}}, controlling the prior setup, the starting values for the
#' MCMC chains, the number of independent MCMC chains, thinning and other expert
#' settings.
#' @return The value returned is a list object of class \code{svdraws_roll}
#' holding a list item for every time window. The elements of these list items are
#' \item{indices}{a list object containing two elements: \code{train} is the vector
#' of indices used for fitting the model, and \code{test} is the vector of indices
#' used for prediction. The latter is mainly useful if a \code{designmatrix} is provided.}
#' \item{quantiles}{the input parameter \code{calculate_quantiles}.}
#' \item{refit_every}{the input parameter \code{refit_every}.}
#' \item{predictive_likelihood}{present only if \code{calculate_predictive_likelihood}
#' is \code{TRUE}. Then it is a number, the expected predictive density
#' of the observation. The expecation is taken over the joint \code{n_ahead} predictive
#' distribution of all model parameters.}
#' \item{predictive_quantile}{present only if \code{calculate_quantile} is a non-empty
#' vector. Then it is a vector of quantiles from the \code{n_ahead} predictive
#' distribution of \code{y}. It is based on MCMC simulation by using \code{\link{predict}}.}
#' \item{fit}{present only if \code{keep_draws} is \code{TRUE}. Then it is an
#' \code{svdraws} object as returned by \code{\link{svsample}}.}
#' \item{prediction}{present only if \code{keep_draws} is \code{TRUE}. Then it is an
#' \code{svpredict} object as returned by \code{\link{predict.svdraws}}.}
#' 
#' To display the output, use \code{print} and \code{summary}. The
#' \code{print} method simply prints a short summary of the setup;
#' the \code{summary} method displays the summary statistics
#' of the backtesting.
#' @note
#' The function executes \code{\link{svsample}} \code{(length(y) - arorder - n_ahead - n_start + 1) \%/\% refit_every} times.
#' @seealso \code{\link{svsim}}, \code{\link{specify_priors}}, \code{\link{svsample}}
#' @keywords models ts
#' @example inst/examples/svsample_roll.R
#' @export
svsample_roll <- function (y, designmatrix = NA,
                           n_ahead = 1, forecast_length = 500,
                           n_start = NULL, refit_every = 1,
                           refit_window = c("moving", "expanding"),
                           calculate_quantile = c(0.01),
                           calculate_predictive_likelihood = TRUE,
                           keep_draws = FALSE,
                           parallel = c("no", "multicore", "snow"),
                           n_cpus = 1L, cl = NULL,
                           ...) {

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

  ## regression
  arorder <- 0L
  if (any(is.na(designmatrix))) {
    designmatrix <- matrix(NA_real_)
  } else {
    if (any(grep("ar[0-9]+$", as.character(designmatrix)[1]))) {
      arorder <- as.integer(gsub("ar", "", as.character(designmatrix)))
      if (length(y) <= arorder + 1) {
        stop("Time series 'y' is too short for this AR process.")
      }
    } else if (is.character(designmatrix)) {
      stop("Argument 'designmatrix' must be a numeric matrix or an AR-specification.")
    } else {
      assert_numeric(designmatrix, "The processed argument 'designmatrix'")
      if (!is.matrix(designmatrix)) {
        designmatrix <- matrix(designmatrix, ncol = 1)
      }
      if (NROW(designmatrix) != length(y)) {
        stop("Number of columns of argument 'designmatrix' must be equal to length(y).")
      }
    }
  }

  refit_window <- match.arg(refit_window)
  if (!identical(refit_every, 1)) {
    stop("Parameter 'refit_every' has to be 1")
  }

  assert_numeric(n_ahead, "Parameter 'n_ahead'")
  n_ahead <- as.integer(n_ahead)
  assert_single(n_ahead, "Parameter 'n_ahead'")
  assert_positive(n_ahead, "Parameter 'n_ahead'")

  assert_numeric(forecast_length, "Parameter 'forecast_length'")
  forecast_length <- as.integer(forecast_length)
  assert_single(forecast_length, "Parameter 'forecast_length'")
  assert_positive(forecast_length, "Parameter 'forecast_length'")

  if (!is.null(n_start)) {
    assert_numeric(n_start, "Parameter 'n_start'")
    n_start <- as.integer(n_start)
    assert_single(n_start, "Parameter 'n_start'")
    assert_gt(n_start, 2, "Parameter 'n_start'")
  } else {
    n_start <- length(y) - forecast_length + 1
  }

  assert_numeric(calculate_quantile, "Parameter 'calculate_quantile'")
  assert_positive(calculate_quantile, "Parameter 'calculate_quantile'")
  assert_lt(calculate_quantile, 1, "Parameter 'calculate_quantile'")

  assert_logical(calculate_predictive_likelihood, "Parameter 'calculate_predictive_likelihood'")
  assert_single(calculate_predictive_likelihood, "Parameter 'calculate_predictive_likelihood'")

  assert_logical(keep_draws, "Parameter 'keep_draws'")
  assert_single(keep_draws, "Parameter 'keep_draws'")

  ## parallel (strongly influenced by package 'boot')
  parallel <- match.arg(parallel)
  have_mc <- FALSE
  have_snow <- FALSE
  if (parallel != "no") {
    if (parallel == "multicore") {
      have_mc <- .Platform$OS.type != "windows"
    } else if (parallel == "snow") {
      have_snow <- TRUE
    }
    if (!have_mc && !have_snow) {
      n_cpus <- 1L
    }
    requireNamespace("parallel")
  }
  create_cluster <- isTRUE(is.null(cl))
  have_cluster <- isTRUE(inherits(cl, "SOCKcluster"))
  if (have_snow && !create_cluster && !have_cluster) {
    warning("Unknown type of input in parameter 'cl'. Should be either NULL or of class 'cluster' created by the parallel package. Turning off parallelism")
    cl <- NULL
    have_mc <- FALSE
    have_snow <- FALSE
    n_cpus <- 1L
  }
  if (have_snow && n_cpus == 1L && !have_cluster) {
    warning("Inconsistent settings for parallel execution: snow parallelism combined with n_cpus == 1 and no supplied cluster object. Turning off parallelism")
    have_snow <- FALSE
  }
  if (have_mc && n_cpus == 1L) {
    warning("Inconsistent settings for parallel execution: multicore parallelism combined with n_cpus == 1. Turning off parallelism")
    have_mc <- FALSE
  }
  n_workers <- if (have_snow) {
    if (have_cluster) {
      length(cl)
    } else {
      n_cpus
    }
  } else {
    1L
  }
  n_chains <- n_cpus

  # Define accessors for the windows
  choose_indices <- function (i) {
    n <- length(y)
    train_begin <- if (refit_window == "expanding") {
      1
    } else {
      1 + (i-1) * refit_every
    }
    train_end <- train_begin + n_start - 2 + arorder
    list(train = seq(train_begin, train_end),
         test = seq(train_end+1, train_end+n_ahead))
  }

  # Define parallel function
  sampling_function <- function (i) {
    ind <- choose_indices(i)
    arguments$y <- y[ind$train]
    arguments$designmatrix <- if (any(is.na(designmatrix)) || is.character(designmatrix)) designmatrix else designmatrix[ind$train, ]
    arguments$parallel <- "no"
    arguments$print_progress <- n_windows
    svdraws <- do.call(stochvol::svsample, arguments)
    new_data <- if (is.character(designmatrix) || any(is.na(designmatrix))) {
      NULL
    } else {
      designmatrix[ind$test, ]
    }
    pred <- predict(svdraws, steps = n_ahead, newdata = new_data)
    predobs <- predy(pred, "concat")
    predvol <- predvola(pred, "concat")
    predquant <- if (length(calculate_quantile) > 0) {
      quantile(predobs[, n_ahead, drop=TRUE], prob = calculate_quantile)
    } else {
      NULL
    }
    predlik <- if (calculate_predictive_likelihood) {
      last_test_ind <- tail(ind$test, 1)
      regression_mean <- if (is.character(designmatrix)) {
        if (arorder == 0L) {
          svbeta(svdraws, "concat")
        } else {
          c(1, y[last_test_ind - seq(arorder, 1, by = -1)]) %*% t(svbeta(svdraws, "concat"))
        }
      } else if (any(is.na(designmatrix))) {
        0
      } else {
        designmatrix[last_test_ind, ] %*% t(svbeta(svdraws, "concat"))
      }
      mean(dnorm(y[last_test_ind], mean = regression_mean, sd = predvol[, n_ahead]))
    } else {
      NULL
    }
    ret <- if (keep_draws) {
      attr(svdraws, "args") <- arguments
      attr(pred, "args") <- list(steps = n_ahead, newdata = new_data)
      list(fit = svdraws,
           prediction = pred)
    } else {
      list()
    }
    ret$predicted_quantile <- predquant
    ret$predictive_likelihood <- predlik
    ret$indices <- ind
    ret$quantiles <- calculate_quantile
    ret$refit_window <- refit_window
    ret
  }

  # Call sampler
  n_windows <- (length(y) - arorder - n_ahead - n_start + 1) %/% refit_every
  arguments <- list(...)
  runtime <- system.time(reslist <-
    if ((n_cpus > 1L || have_cluster) && (have_mc || have_snow)) {
      if (have_mc) {
        parallel::mclapply(seq_len(n_windows), sampling_function, mc.cores = n_cpus)
      } else if (have_snow) {
        if (create_cluster) {
          cat("\nStarting cluster...\n")
          cl <- parallel::makePSOCKcluster(rep("localhost", n_workers), outfile = NULL)
        }
        RNGkind(kind = "L'Ecuyer-CMRG")
        parallel::clusterSetRNGStream(cl)
        parallel::clusterEvalQ(cl, library(stochvol))
        parallel::clusterExport(cl, c("y", "designmatrix",
                                      "arorder", "refit_window", "refit_every",
                                      "n_start", "n_ahead", "choose_indices",
                                      "n_windows", "calculate_quantile",
                                      "calculate_predictive_likelihood",
                                      "keep_draws", "arguments"),
                                envir = environment())
        if (create_cluster) {
          reslist <- tryCatch(parallel::parLapply(cl, seq_len(n_windows), sampling_function),
                              finally = {
                                parallel::stopCluster(cl)
                                cat("\nCluster stopped.\n")
                                })
          reslist
        } else {
          parallel::parLapply(cl, seq_len(n_windows), sampling_function)
        }
      }
    } else {
      lapply(seq_len(n_windows), sampling_function)
    }
  )
  class(reslist) <- "svdraws_roll"
  reslist
}

#' @rdname svsample_roll
#' @export
svtsample_roll <- function (y, designmatrix = NA,
                            n_ahead = 1, forecast_length = 500,
                            n_start = NULL, refit_every = 1,
                            refit_window = c("moving", "expanding"),
                            calculate_quantile = c(0.01),
                            calculate_predictive_likelihood = TRUE,
                            keep_draws = FALSE,
                            parallel = c("no", "multicore", "snow"),
                            n_cpus = 1L, cl = NULL,
                            ...) {
  local_svsample_roll <- function (y, designmatrix,
                                   n_ahead, forecast_length,
                                   n_start, refit_every,
                                   refit_window,
                                   calculate_quantile,
                                   calculate_predictive_likelihood,
                                   keep_draws,
                                   parallel,
                                   n_cpus, cl,
                                   priornu = 0.1,  # set prior
                                   ...) {
    svsample_roll(y = y, designmatrix = designmatrix,
                  n_ahead = n_ahead, forecast_length = forecast_length,
                  n_start = n_start, refit_every = refit_every,
                  refit_window = refit_window,
                  calculate_quantile = calculate_quantile,
                  calculate_predictive_likelihood = calculate_predictive_likelihood,
                  keep_draws = keep_draws,
                  parallel = parallel,
                  n_cpus = n_cpus, cl = cl,
                  priornu = priornu,
                  ...)
  }
  local_svsample_roll(y = y, designmatrix = designmatrix,
                      n_ahead = n_ahead, forecast_length = forecast_length,
                      n_start = n_start, refit_every = refit_every,
                      refit_window = refit_window,
                      calculate_quantile = calculate_quantile,
                      calculate_predictive_likelihood = calculate_predictive_likelihood,
                      keep_draws = keep_draws,
                      parallel = parallel,
                      n_cpus = n_cpus, cl = cl,
                      ...)
}

#' @rdname svsample_roll
#' @export
svlsample_roll <- function (y, designmatrix = NA,
                            n_ahead = 1, forecast_length = 500,
                            n_start = NULL, refit_every = 1,
                            refit_window = c("moving", "expanding"),
                            calculate_quantile = c(0.01),
                            calculate_predictive_likelihood = TRUE,
                            keep_draws = FALSE,
                            parallel = c("no", "multicore", "snow"),
                            n_cpus = 1L, cl = NULL,
                            ...) {
  local_svsample_roll <- function (y, designmatrix,
                                   n_ahead, forecast_length,
                                   n_start, refit_every,
                                   refit_window,
                                   calculate_quantile,
                                   calculate_predictive_likelihood,
                                   keep_draws,
                                   parallel,
                                   n_cpus, cl,
                                   priorrho = c(4, 4),  # set prior
                                   ...) {
    svsample_roll(y = y, designmatrix = designmatrix,
                  n_ahead = n_ahead, forecast_length = forecast_length,
                  n_start = n_start, refit_every = refit_every,
                  refit_window = refit_window,
                  calculate_quantile = calculate_quantile,
                  calculate_predictive_likelihood = calculate_predictive_likelihood,
                  keep_draws = keep_draws,
                  parallel = parallel,
                  n_cpus = n_cpus, cl = cl,
                  priorrho = priorrho,
                  ...)
  }
  local_svsample_roll(y = y, designmatrix = designmatrix,
                      n_ahead = n_ahead, forecast_length = forecast_length,
                      n_start = n_start, refit_every = refit_every,
                      refit_window = refit_window,
                      calculate_quantile = calculate_quantile,
                      calculate_predictive_likelihood = calculate_predictive_likelihood,
                      keep_draws = keep_draws,
                      parallel = parallel,
                      n_cpus = n_cpus, cl = cl,
                      ...)
}

#' @rdname svsample_roll
#' @export
svtlsample_roll <- function (y, designmatrix = NA,
                             n_ahead = 1, forecast_length = 500,
                             n_start = NULL, refit_every = 1,
                             refit_window = c("moving", "expanding"),
                             calculate_quantile = c(0.01),
                             calculate_predictive_likelihood = TRUE,
                             keep_draws = FALSE,
                             parallel = c("no", "multicore", "snow"),
                             n_cpus = 1L, cl = NULL,
                             ...) {
  local_svsample_roll <- function (y, designmatrix,
                                   n_ahead, forecast_length,
                                   n_start, refit_every,
                                   refit_window,
                                   calculate_quantile,
                                   calculate_predictive_likelihood,
                                   keep_draws,
                                   parallel,
                                   n_cpus, cl,
                                   priornu = 0.1,  # set prior
                                   priorrho = c(4, 4),  # set prior
                                   ...) {
    svsample_roll(y = y, designmatrix = designmatrix,
                  n_ahead = n_ahead, forecast_length = forecast_length,
                  n_start = n_start, refit_every = refit_every,
                  refit_window = refit_window,
                  calculate_quantile = calculate_quantile,
                  calculate_predictive_likelihood = calculate_predictive_likelihood,
                  keep_draws = keep_draws,
                  parallel = parallel,
                  n_cpus = n_cpus, cl = cl,
                  priornu = priornu,
                  priorrho = priorrho,
                  ...)
  }
  local_svsample_roll(y = y, designmatrix = designmatrix,
                      n_ahead = n_ahead, forecast_length = forecast_length,
                      n_start = n_start, refit_every = refit_every,
                      refit_window = refit_window,
                      calculate_quantile = calculate_quantile,
                      calculate_predictive_likelihood = calculate_predictive_likelihood,
                      keep_draws = keep_draws,
                      parallel = parallel,
                      n_cpus = n_cpus, cl = cl,
                      ...)
}

