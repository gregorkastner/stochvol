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

#' @rdname extractors
#' @export
para <- function(x) {
 x$para
}

#' @rdname extractors
#' @export
latent <- function(x) {
 x$latent
}

#' @rdname extractors
#' @export
latent0 <- function(x) {
 x$latent0
}

#' @rdname extractors
#' @export
priors <- function(x) {
 x$priors
}

#' @rdname extractors
#' @export
thinning <- function(x) {
 x$thinning
}

#' @rdname extractors
#' @export
runtime <- function(x) {
 x$runtime
}

#' @rdname extractors
#' @export
sampled_parameters <- function(x) {
  res <- c()
  pr <- priors(x)
  if (!inherits(pr$mu, "sv_constant")) {
    res <- c(res, "mu")
  }
  if (!inherits(pr$phi, "sv_constant")) {
    res <- c(res, "phi")
  }
  if (!inherits(pr$sigma2, "sv_constant")) {
    res <- c(res, "sigma")
  }
  if (!inherits(pr$nu, "sv_constant") && !inherits(pr$nu, "sv_infinity")) {
    res <- c(res, "nu")
  }
  if (!inherits(pr$rho, "sv_constant")) {
    res <- c(res, "rho")
  }
  res
}

#' @export
"[.svdraws" <- function (x, n, drop=FALSE) {
  if (drop) {
    stop("Parameter 'drop' is only allowed to be FALSE")
  }
  res <- x
  mcmcnames <- c("para", "latent", "latent0", "beta", "tau")
  for (mcmcname in mcmcnames) {
    if (!is.null(res[[mcmcname]])) {
      if (length(n) > 1) {
        res[[mcmcname]] <- mcmc.list(res[[mcmcname]][n])
      } else {
        res[[mcmcname]] <- mcmc.list(list(res[[mcmcname]][[n]]))
      }
    }
  }
  res
}

#' @export
"[.svpredict" <- function (x, n, drop=FALSE) {
  if (drop) {
    stop("Parameter 'drop' is only allowed to be FALSE")
  }
  res <- x
  mcmcnames <- c("y", "h")
  for (mcmcname in mcmcnames) {
    if (!is.null(res[[mcmcname]])) {
      if (length(n) > 1) {
        res[[mcmcname]] <- mcmc.list(res[[mcmcname]][n])
      } else {
        res[[mcmcname]] <- mcmc.list(list(res[[mcmcname]][[n]]))
      }
    }
  }
  res
}

#' @rdname extractors
#' @export
vola <- function (x) {
  vol <- list()
  for (chain in seq_len(coda::nchain(x$latent))) {
    tau_chain <- if (!is.null(x$tau)) x$tau[[chain]] else 1
    pars <- mcpar(x$latent[[chain]])
    vol[[chain]] <- coda::mcmc(exp(x$latent[[chain]]/2) * sqrt(tau_chain), pars[1], pars[2], pars[3])
  }
  coda::mcmc.list(vol)
}

#' @method as.array svdraws
#' @export
as.array.svdraws <- function (x, ...) {
  aperm(simplify2array(lapply(seq_along(para(x)), function (i) cbind(x$para[[i]], x$latent0[[i]], x$latent[[i]], x$tau[[i]], x$beta[[i]])), higher = TRUE), c(1, 3, 2))
}

contains_beta <- function(x) {
  NROW(x$beta) > 0
}

contains_tau <- function(x) {
  NROW(x$tau) > 0
}

contains_indicators <- function(x) {
  NROW(x$indicators) > 0
}

join_mcmclist <- function (x) {
  do.call(rbind, x)
}

# x is svdraws with a single chain
# output: svdraws-like object but with mcmc objects instead of the mcmc.list objects
flatten <- function (x) {
  res <- x
  mcmcnames <- c("para", "latent", "latent0", "beta", "tau")
  for (mcmcname in mcmcnames) {
    if (!is.null(res[[mcmcname]])) {
      res[[mcmcname]] <- res[[mcmcname]][[1]]
    }
  }
  class(res) <- NULL
  res
}

#' Updating the Summary of MCMC Draws
#' 
#' Creates or updates a summary of an \code{svdraws} object.
#' 
#' \code{updatesummary} will always calculate the posterior mean and the
#' posterior standard deviation of the raw draws and some common
#' transformations thereof. Moroever, the posterior quantiles, specified by the
#' argument \code{quantiles}, are computed. If \code{esspara} and/or
#' \code{esslatent} are \code{TRUE}, the corresponding effective sample size
#' (ESS) will also be included.
#' 
#' @param x \code{svdraws} object.
#' @param quantiles numeric vector of posterior quantiles to be computed. The
#' default is \code{c(0.05, 0.5, 0.95)}.
#' @param esspara logical value which indicates whether the effective sample
#' size (ESS) should be calculated for the \emph{parameter draws}. This is
#' achieved by calling \code{\link[coda]{effectiveSize}} from the \code{coda}
#' package. The default is \code{TRUE}.
#' @param esslatent logical value which indicates whether the effective sample
#' size (ESS) should be calculated for the \emph{latent log-volatility} draws.
#' This is achieved by calling \code{\link[coda]{effectiveSize}} from the
#' \code{coda} package. The default is \code{FALSE}, because this can be quite
#' time-consuming when many latent variables are present.
#' @return The value returned is an updated list object of class \code{svdraws}
#' holding \item{para}{\code{mcmc} object containing the \emph{parameter} draws
#' from the posterior distribution.} \item{latent}{\code{mcmc} object
#' containing the \emph{latent instantaneous log-volatility} draws from the
#' posterior distribution.} \item{latent0}{\code{mcmc} object containing the
#' \emph{latent initial log-volatility} draws from the posterior distribution.}
#' \item{y}{argument \code{y}.} \item{runtime}{\code{"proc_time"} object
#' containing the run time of the sampler.} \item{priors}{\code{list}
#' containing the parameter values of the prior distribution, i.e. the
#' arguments \code{priormu}, \code{priorphi}, \code{priorsigma} (and
#' potentially \code{nu}).} \item{thinning}{\code{list} containing the thinning
#' parameters, i.e. the arguments \code{thinpara}, \code{thinlatent} and
#' \code{keeptime}.} \item{summary}{\code{list} containing a collection of
#' summary statistics of the posterior draws for \code{para}, \code{latent},
#' and \code{latent0}.}
#' 
#' To display the output, use \code{print}, \code{summary} and \code{plot}. The
#' \code{print} method simply prints the posterior draws (which is very likely
#' a lot of output); the \code{summary} method displays the summary statistics
#' currently stored in the object; the \code{plot} method gives a graphical
#' overview of the posterior distribution by calling \code{\link{volplot}},
#' \code{\link{traceplot}} and \code{\link{densplot}} and displaying the
#' results on a single page.
#' @note \code{updatesummary} does not actually overwrite the object's current
#' summary, but in fact creates a new object with an updated summary. Thus,
#' don't forget to overwrite the old object if this is want you intend to do.
#' See the examples below for more details.
#' @seealso \code{\link{svsample}}
#' @keywords utilities
#' @example inst/examples/updatesummary.R
#' @export
updatesummary <- function(x, quantiles = c(.05, .5, .95), esspara = TRUE, esslatent = FALSE) {
  sampled_para <- sampled_parameters(x)

  matrixsummary <- function (x, quants = quantiles) {
    cbind(mean = colMeans(x), sd = apply(x, 2, sd), t(apply(x, 2, quantile, prob = quants)))
  }

  # x is an mcmc.list object containing matrices
  # transformation is a singleton named list with a vectorized function of interest for the variable(s)
  summaryfunction <- function(x, quants = quantiles, ess = TRUE, transformation = NULL) {
    stopifnot(inherits(x, "mcmc.list"))
    x_join <- join_mcmclist(x)
    res <- matrixsummary(x_join, quants = quants)
    if (ess) {
      res <- cbind(res, ESS = as.numeric(coda::effectiveSize(x)))  # computes the sum of effective sample sizes over the chains
    }
    if (!is.null(transformation)) {
      stopifnot(length(transformation) == 1)
      name <- names(transformation)
      stopifnot(length(name) == 1)
      x_trans <- transformation[[name]](x_join)
      toadd <- colMeans(x_trans)
      toadd <- cbind(toadd, apply(x_trans, 2, sd))
      colnames(toadd) <- paste0(c("mean", "sd"), "(", name, ")")
      res <- cbind(res, toadd)
    }
    res
  }

  res <- list()

  res$para <- summaryfunction(x$para[, sampled_para, drop=FALSE], ess = esspara)
  if ("mu" %in% sampled_para) {
    toadd <- matrixsummary(exp(join_mcmclist(x$para[, "mu", drop=FALSE])/2))
    if (esspara) {
      toadd <- cbind(toadd, ESS = res$para["mu", "ESS"])
    }
    rownames(toadd) <- "exp(mu/2)"
    res$para <- rbind(res$para, toadd)
  }
  if ("sigma" %in% sampled_para) {
    toadd <- matrixsummary(join_mcmclist(x$para[, "sigma", drop=FALSE])^2)
    if (esspara) {
      toadd <- cbind(toadd, ESS = res$para["sigma", "ESS"])
    }
    rownames(toadd) <- "sigma^2"
    res$para <- rbind(res$para, toadd)
  }

  res$latent0 <- summaryfunction(x$latent0, ess = esslatent, transformation = list("exp(h_t/2)" = function (x) exp(x/2)))
  res$latent <- summaryfunction(x$latent, ess = esslatent, transformation = list("exp(h_t/2)" = function (x) exp(x/2)))

  vol <- vola(x)
  res$sd <- summaryfunction(vol, ess = esslatent)
  rownames(res$sd) <- gsub("h", "sd", rownames(vol))

  if (contains_beta(x)) {
    res$beta <- summaryfunction(x$beta, ess = esspara)
  }

  x$summary <- res
  x
}

#' @export
summary.svdraws <- function (object, showpara = TRUE, showlatent = TRUE, ...) {
  ret <- vector("list")
  class(ret) <- "summary.svdraws"
  ret$mcp <- mcpar(para(object)[[1]])
  ret$mcl <- mcpar(latent(object)[[1]])
  ret$priors <- priors(object)
  if (isTRUE(showpara)) {
    ret$para <- para(object$summary)
  }
  if (isTRUE(showlatent)) {
    ret$latent <- latent(object$summary)
    ret$latent <- rbind("h_0" = latent0(object$summary), ret$latent)
  }
  if (isTRUE(object$resampled)) {
    ret$resampled <- list(para = list(max_entropy = log(length(object$correction_weight_para[[1]])),
                                      entropy = min(sapply(object$correction_weight_para, function (w) -sum(w * log(w))))),
                          latent = list(max_entropy = log(length(object$correction_weight_latent[[1]])),
                                        entropy = min(sapply(object$correction_weight_latent, function (w) -sum(w * log(w))))),
                          same_resampling = object$thinning$para == object$thinning$latent)
  }
  ret
}

#' @export
print.svdraws <- function (x, showpara = TRUE, showlatent = FALSE, ...) {
  print(summary(x, showpara = showpara, showlatent = showlatent), ...)
  invisible(x)
}

#' @method print summary.svdraws
#' @export
print.summary.svdraws <- function (x, ...) { 
  cat("\nSummary of ", x$mcp[2]-x$mcp[1]+x$mcp[3], " MCMC draws after a burn-in of ", x$mcp[1]-x$mcp[3], ".\n", sep="")
  print(x$priors)
  
  if (exists("para", x)) {
    cat("\nPosterior draws of parameters (thinning = ", x$mcp[3], "):\n", sep='')
    print(x$para, digits=2, ...)
  }

  if (exists("latent", x)) {
    cat("\nPosterior draws of initial and contemporaneous latents (thinning = ", x$mcl[3], "):\n", sep='')
    print(x$latent, digits=2, ...)
  }

  if (exists("resampled", x)) {
    cat("\nPosterior draws from the auxiliary SV model were re-sampled after the MCMC\nprocedure ended to correct for model mis-specification.\n")
    print_resampling <- function (entropy_list, text) {
      if (entropy_list$entropy >= entropy_list$max_entropy*0.99) {
        cat("  Re-sampling of ", text, " was of little practical importance:\n", sep = "")
      } else {
        cat("  Re-sampling of ", text, " might have made practical a difference:\n", sep = "")
      }
      cat("    - max reachable entropy for this sample size: ", entropy_list$max_entropy, ",\n    - entropy of the re-sampling weight distribution: ", entropy_list$entropy, "\n", sep = "")
    }
    if (x$resampled$same_resampling) {
      print_resampling(x$resampled$para, "parameters and latents")
    } else {
      print_resampling(x$resampled$para, "parameters")
      print_resampling(x$resampled$latent, "latents")
    }
  }
  invisible(x)
}

#' @export
residuals.svdraws <- function(object, type = "mean", ...) {
  if (!inherits(object, "svdraws")) stop("This function expects an 'svdraws' object.")
  if (coda::nchain(para(object)) > 1) stop("Multiple chains in the svdraws object: not yet implemented")
  if (!type %in% c("mean", "median")) stop("Argument 'type' must currently be either 'mean' or 'median'.")

  if (object$thinning$time != 'all') stop("Not every point in time has been stored ('keeptime' was not set to 'all' during sampling), thus residuals cannot be extracted.")

  y <- as.vector(object$y)
  if (contains_beta(object)) {
    y <- y - object$designmatrix %*% t(object$beta)
  }

  if (type == "mean") {
    res <- rowMeans(y[seq(1, length(y))] / exp(t(object$latent)/2))
  }

  if (type == "median") {
    res <- apply(y[seq(1, length(y))] / exp(t(object$latent)/2), 1, median)
  }

  names(res) <- sub("h", "r", colnames(object$latent))
  class(res) <- "svresid"
  attr(res, "type") <- type

  # Also return posterior mean/median of df parameter if terr = TRUE
  if ("nu" %in% sampled_parameters(object)) {
    attr(res, "nu") <- get(type)(object$para[,"nu"])
  }

  res
}


#' Prediction of Future Returns and Log-Volatilities
#' 
#' Simulates draws from the predictive density of the returns and the latent log-volatility
#' process. The same mean model is used for prediction as was used for fitting, which is
#' either a) no mean parameter, b) constant mean, c) AR(k) structure, or d) general
#' Bayesian regression. In the last case, new regressors need to be provided for prediction.
#' 
#' 
#' @param object \code{svdraws} or \code{svldraws} object.
#' @param steps \emph{optional} single number, coercible to integer. Denotes the number of
#' steps to forecast.
#' @param newdata \emph{only in case d) of the description} corresponds to input
#' parameter \code{designmatrix} in \code{\link{svsample}}.
#' A matrix of regressors with number of rows equal to parameter \code{steps}.
#' @param \dots currently ignored.
#' @return Returns an object of class \code{svpredict}, a list containing
#' three elements:
#' \item{vol}{\code{mcmc.list} object of simulations from the predictive density of the standard deviations \code{sd_(n+1),...,sd_(n+steps)}}
#' \item{h}{\code{mcmc.list} object of simulations from the predictive density of \code{h_(n+1),...,h_(n+steps)}}
#' \item{y}{\code{mcmc.list} object of simulations from the predictive density of \code{y_(n+1),...,y_(n+steps)}}
#' @note You can use the resulting object within \code{\link{plot.svdraws}} (see example below), or use
#' the list items in the usual \code{coda} methods for \code{mcmc} objects to
#' print, plot, or summarize the predictions.
#' @seealso \code{\link{plot.svdraws}}, \code{\link{volplot}}.
#' @keywords ts
#' @example inst/examples/predict.R
#' @export
predict.svdraws <- function(object, steps = 1L, newdata = NULL, ...) {
  if (!(inherits(object, "svdraws"))) stop("Argument 'object' must be of class 'svdraws'.")

  steps <- as.integer(steps)
  if (steps < 1) stop("Argument 'steps' must be greater or equal to 1.")
  # Error checking for the mean model
  arorder <- 0  # AR(0) means either a constant mean or a no mean model here
  if (object$meanmodel == "none") {
    if (!is.null(newdata)) warning("No regression coefficients were used when estimating the model. Omitting 'newdata'.")
    regressors <- function (y, newdata, stepind) 0
  } else if (object$meanmodel == "matrix") {
    if (is.null(newdata)) stop("Regressors needed for prediction. Please provide regressors through parameter 'newdata'.")
    newdata <- as.matrix(newdata)
    if (is.null(object$beta) || coda::nvar(object$beta[[1]]) != NCOL(newdata)) stop(paste0("The number of fitted regression coefficients (", coda::nvar(object$beta[[1]]), ") does not equal the number of given regressors (", NCOL(newdata), ")."))
    if (NROW(newdata) != steps) stop("The size of the design matrix (", NROW(newdata), " rows) does not match the number of steps to predict (", steps, ").")
    regressors <- function (y, newdata, stepind) matrix(newdata[stepind, ], nrow = NROW(y), ncol = NCOL(newdata), byrow = TRUE)  # matches the format in the ar* case
  } else if (object$meanmodel == "constant") {
    if (!is.null(newdata)) warning("Constant mean was assumed when estimating the model. Omitting 'newdata'.")
    regressors <- function (y, newdata, stepind) 1
  } else if (any(grep("ar[1-9][0-9]*", object$meanmodel))) {
    arorder <- as.integer(gsub("ar", "", object$meanmodel))
    if (!is.null(newdata)) warning(paste0("An AR(", arorder, ") mean was assumed when estimating the model. Omitting 'newdata'."))
    regressors <- function (y, newdata, stepind) cbind(1, y[,seq.int(stepind, stepind-1+arorder),drop=FALSE])
  } else {
    stop("Unknown mean model. Please contact the developer.")
  }

  object_1 <- flatten(object[1])
  thinlatent <- object_1$thinning$latent
  thinpara <- object_1$thinning$para
  if (thinpara != thinlatent) {
    warning("Thinning of parameters is different from thinning of latent variables. Trying to sort this out.")  # TODO use lowest common multiple
    if (thinpara %% thinlatent == 0) {
      usepara <- seq_len(NROW(para(object_1)))
      uselatent <- seq(thinpara %/% thinlatent, NROW(latent(object_1)), by = thinpara %/% thinlatent)
    } else if (thinlatent %% thinpara == 0) {
      uselatent <- seq_len(NROW(latent(object_1)))
      usepara <- seq(thinlatent %/% thinpara, NROW(para(object_1)), by = thinlatent %/% thinpara)
    } else stop("Incompatible thinning parameters. Prediction currently not implemented.")
  } else {
    usepara <- uselatent <- seq.int(NROW(para(object_1)))
  }

  ret <- list(y = list(), h = list(), vol = list())
  for (chain in seq_len(coda::nchain(para(object)))) {
    object_i <- flatten(object[chain])
    mu <- object_i$para[usepara,"mu"]
    phi <- object_i$para[usepara,"phi"]
    sigma <- object_i$para[usepara,"sigma"]
    rho <- object_i$para[usepara,"rho"]
    nu <- object_i$para[usepara,"nu"]
    heavy_tailed <- is.finite(tail(nu, 1))
    hlast <- object_i$latent[uselatent,NCOL(object_i$latent)]
    taulast <- if (heavy_tailed) object_i$tau[uselatent,NCOL(object_i$tau)] else 1
    ylast <- as.vector(tail(object_i$y, 1))

    mythin <- max(thinpara, thinlatent)
    len <- length(usepara)
    volpred <- matrix(as.numeric(NA), nrow=len, ncol=steps)
    hpred <- matrix(as.numeric(NA), nrow=len, ncol=steps)
    ypred <- matrix(as.numeric(NA), nrow=len, ncol=steps+arorder)

    ypred[,seq_len(arorder)] <- matrix(as.vector(tail(object_i$y, arorder)), nrow=len, ncol=arorder, byrow=TRUE)  # AR(x) helper columns

    betacoeff <- if (contains_beta(object_i)) {
      if (arorder > 0) object_i$beta[usepara, c(1, rev(seq_len(NCOL(object_i$beta)-1))+1), drop=FALSE]
      else object_i$beta[usepara, , drop=FALSE]
    } else matrix(0)

    resilast <- if (object_i$meanmodel == "none") {  # last fitted residual
      ylast*exp(-hlast/2)/sqrt(taulast)
    } else {  # if mean regression
      (ylast - colSums(object_i$designmatrix[NROW(object_i$designmatrix),]*t(betacoeff)))*exp(-hlast/2)/sqrt(taulast)  # recycles the last row of the design matrix
    }

    hpred[,1] <- mu+phi*(hlast-mu) + sigma*(rho*resilast + sqrt(1-rho^2)*rnorm(len))
    tau <- if (heavy_tailed) 1/rgamma(len, shape=.5*nu, rate=.5*(nu-2)) else 1
    volpred[,1] <- sqrt(tau) * exp(hpred[,1]/2)
    if (steps > 1) {
      for (i in seq.int(from=2, to=steps)) {
        resi <- rnorm(len)
        incr <- rho*resi + sqrt(1-rho^2)*rnorm(len)
        ypred[,i-1+arorder] <- rowSums(regressors(ypred, newdata, i-1) * betacoeff) + volpred[,i-1]*resi
        hpred[,i] <- mu + phi*(hpred[,i-1] - mu) + sigma*incr
        tau <- if (heavy_tailed) 1/rgamma(len, shape=.5*nu, rate=.5*(nu-2)) else 1
        volpred[,i] <- sqrt(tau) * exp(hpred[,i]/2)
      }
    }
    ypred[,steps+arorder] <- rowSums(regressors(ypred, newdata, steps) * betacoeff) + volpred[,steps]*rnorm(len)

    ypred <- ypred[, setdiff(seq_len(NCOL(ypred)), seq_len(arorder)), drop=FALSE]  # remove temporary AR(x) helper columns
    lastname <- tail(colnames(object_i$latent), 1)
    lastnumber <- as.integer(gsub("h_", "", lastname))
    colnames(volpred) <- paste0("sd_", seq(lastnumber + 1, lastnumber + steps))
    colnames(hpred) <- paste0("h_", seq(lastnumber + 1, lastnumber + steps))
    colnames(ypred) <- paste0("y_", seq(lastnumber + 1, lastnumber + steps))
    ret$vol[[chain]] <- coda::mcmc(volpred, start=mythin, end=len*mythin, thin=mythin)
    ret$h[[chain]] <- coda::mcmc(hpred, start=mythin, end=len*mythin, thin=mythin)
    ret$y[[chain]] <- coda::mcmc(ypred, start=mythin, end=len*mythin, thin=mythin)
  }

  ret$vol <- mcmc.list(ret$vol)
  ret$h <- mcmc.list(ret$h)
  ret$y <- mcmc.list(ret$y)
  class(ret) <- c("svpredict")
  ret
}

#' @export
print.svpredict <- function (x, ...) {
  cat("'svpredict' object containing predicted values for variables:\n")
  cat("  - predicted observations:        ", paste(coda::varnames(x$y[[1]]), collapse = ", "), "\n", sep = "")
  cat("  - predicted standard deviations: ", paste(coda::varnames(x$vol[[1]]), collapse = ", "), "\n", sep = "")
  cat("  - predicted latent variables:    ", paste(coda::varnames(x$h[[1]]), collapse = ", "), "\n", sep = "")
  invisible(x)
}

