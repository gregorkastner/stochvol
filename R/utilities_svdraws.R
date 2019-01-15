#' Computes (de-meaned) log returns.
#' 
#' Small utlity function returning either \code{diff(log(x))} in case the
#' argument \code{demean} is set to \code{FALSE}, or \code{diff(log(x)) -
#' mean(diff(log(x)))} in case that \code{demean} is \code{TRUE}.
#' 
#' 
#' @param x Real-valued vector.
#' @param demean A single logical value indicating whether the returns should
#' be de-meaned. Defaults to \code{FALSE}.
#' @return A vector of length \code{length(x) - 1}, containing (de-meaned)
#' returns.
#' @author Gregor Kastner \email{gregor.kastner@@wu.ac.at}
#' @keywords utilities
#' @export
logret <- function(x, demean = FALSE) {
  logretx <- tail(diff(log(x)), length(x)-1)
  if (isTRUE(demean)) logretx <- logretx - mean(logretx)
  logretx
}


#' @export
para <- function(x) {
 x$para
}


#' @export
latent <- function(x) {
 x$latent
}


#' @export
latent0 <- function(x) {
 x$latent0
}


#' @export
priors <- function(x) {
 x$priors
}


#' @export
thinning <- function(x) {
 x$thinning
}


#' @export
runtime <- function(x) {
 x$runtime
}
# TODO common docs for the above functions


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
#' \code{thintime}.} \item{summary}{\code{list} containing a collection of
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
#' @author Gregor Kastner \email{gregor.kastner@@wu.ac.at}
#' @seealso \code{\link{svsample}}
#' @keywords utilities
#' @examples
#' 
#' ## Here is a baby-example to illustrate the idea.
#' ## Simulate an SV time series of length 51 with default parameters:
#' sim <- svsim(51)
#' 
#' ## Draw from the posterior (but save only every fifth point in time):
#' res <- svsample(sim$y, draws = 7000, thintime = 5, priorphi = c(10, 1.5))
#' 
#' ## Check out the results:
#' summary(res)
#' plot(res)
#' 
#' ## Look at other quantiles and calculate ESS of latents:
#' newquants <- c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
#' res <- updatesummary(res, quantiles = newquants, esslatent = TRUE)
#' 
#' ## See the difference?
#' summary(res)
#' plot(res)
#' 
#' @export
updatesummary <- function(x, quantiles = c(.05, .5, .95), esspara = TRUE, esslatent = FALSE) {

  # Check if conditional t errors are used
  terr <- "nu" %in% colnames(x$para)

  summaryfunction <- function(x, quants = quantiles, ess = TRUE) {
    if (ess) {
      c(mean = mean(x), sd = sd(x), quantile(x, quantiles),
        ESS = as.numeric(effectiveSize(x)))
    } else {
      c(mean = mean(x), sd = sd(x), quantile(x, quantiles))
    }
  }

  res <- list()

  res$para <- t(apply(x$para, 2, summaryfunction, ess = esspara))
  res$para <- rbind(res$para, "exp(mu/2)" = c(summaryfunction(exp(x$para[,"mu"]/2), ess=FALSE), res$para["mu", "ESS"]))
  res$para <- rbind(res$para, "sigma^2" = c(summaryfunction(x$para[,"sigma"]^2, ess=FALSE), res$para["sigma", "ESS"]))

  res$latent <- t(apply(x$latent, 2, summaryfunction, ess = esslatent))
  vol <- exp(x$latent/2)
  res$latent <- cbind(res$latent, "mean(exp(h_t/2))" = colMeans(vol))
  res$latent <- cbind(res$latent, "sd(exp(h_t/2))" = apply(vol, 2, sd))

  if (!is.null(x$latent0))
    res$latent0 <- c(summaryfunction(x$latent0, ess = esslatent), "mean(exp(h_t/2))" = mean(exp(x$latent0/2)), "sd(exp(h_t/2))" = sd(exp(x$latent0/2)))

  if (terr && x$thinning$para == x$thinning$latent) {
    vol <- sqrt((exp(x$latent) * x$para[,"nu"]) / (x$para[,"nu"] - 2))
    res$sd <- t(apply(vol, 2, summaryfunction, ess = esslatent))
    rownames(res$sd) <- gsub("h", "sd", rownames(res$latent))
  }

  if (exists("beta", x)) res$beta <- t(apply(x$beta, 2, summaryfunction, ess = esspara))

  x$summary <- res
  x
}

#' @export
summary.svldraws <- function(object, showpara = TRUE, showlatent = TRUE, ...) {
  ret <- vector("list")
  class(ret) <- "summary.svldraws"
  ret$mcp <- mcpar(para(object))
  ret$mcl <- mcpar(latent(object))
  ret$priors <- priors(object)
  if (isTRUE(showpara)) ret$para <- para(object$summary)
  if (isTRUE(showlatent)) ret$latent <- latent(object$summary)
  ret
}

#' @export
summary.svdraws <- function (object, showpara = TRUE, showlatent = TRUE, ...) {
  ret <- summary.svldraws(object, showpara = showpara, showlatent = showlatent)
  class(ret) <- "summary.svdraws"
  if (isTRUE(showlatent)) ret$latent <- rbind("h_0" = latent0(object$summary), ret$latent)
  ret
}

#' @method print summary.svdraws
#' @export
print.summary.svdraws <- function(x, ...) { 
  cat("\nSummary of ", x$mcp[2]-x$mcp[1]+x$mcp[3], " MCMC draws after a burn-in of ", x$mcp[1]-x$mcp[3], ".\n", sep="")
  cat("Prior distributions:\n")
  cat("mu        ~ Normal(mean = ", x$priors$mu[1], ", sd = ", x$priors$mu[2], ")\n", sep="")
  cat("(phi+1)/2 ~ Beta(a0 = ", x$priors$phi[1], ", b0 = ", x$priors$phi[2], ")\n", sep="")
  cat("sigma^2   ~ ", x$priors$sigma, " * Chisq(df = 1)\n", sep="")

  if ("nu" %in% colnames(x$priors)) cat("nu        ~ Unif(lower = ", x$priors$nu[1], ", upper = ", x$priors$nu[2], ")\n", sep="")
  if ("rho" %in% colnames(x$priors)) cat("(phi+1)/2 ~ Beta(a0 = ", x$priors$phi[1], ", b0 = ", x$priors$phi[2], ")\n", sep="")

  if (exists("para", x)) {
    cat("\nPosterior draws of parameters (thinning = ", x$mcp[3], "):\n", sep='')
    print(x$para, digits=2, ...)
  }

  if (exists("latent", x)) {
    cat("\nPosterior draws of initial and contemporaneous latents (thinning = ", x$mcl[3], "):\n", sep='')
    print(x$latent, digits=2, ...)
  }
  invisible(x)
}

#' @method print summary.svldraws
#' @export
print.summary.svldraws <- function (x, ...) {
  ret <- print.summary.svdraws(x = x, ...)
  invisible(ret)
}

#' @export
print.svdraws <- function(x, showpara = TRUE, showlatent = TRUE, ...) {
  if (isTRUE(showpara)) {
    cat("\n*** Posterior draws of parameters ***\n")
    print(para(x), ...)
  }

  if (isTRUE(showlatent)) {
    if (!is.null(latent0(x))) {
      cat("\n*** Posterior draws of initial latent variable h_0 ***\n")
      print(latent0(x), ...)
    }
    cat("\n*** Posterior draws of contemporaneous latent variables h_t ***\n")
    print(latent(x), ...)
  }
  invisible(x)
}

#' @export
print.svldraws <- function (x, showpara = TRUE, showlatent = TRUE, ...) {
  print.svdraws(x = x, showpara = showpara, showlatent = showlatent, ...)
}

#' @export
residuals.svdraws <- function(object, type = "mean", ...) {
  if (!inherits(object, "svdraws")) stop("This function expects an 'svdraws' object.")
  if (!type %in% c("mean", "median")) stop("Argument 'type' must currently be either 'mean' or 'median'.")

  if (object$thinning$time != 1) warning("Not every point in time has been stored ('thintime' was set to a value unequal to 1 during sampling), thus only some residuals have been extracted.")

  if (type == "mean") {
    res <- rowMeans(as.numeric(object$y)[seq(1, length(object$y), by=object$thinning$time)] / exp(t(object$latent)/2))
  }

  if (type == "median") {
    res <- apply(as.numeric(object$y)[seq(1, length(object$y), by=object$thinning$time)] / exp(t(object$latent)/2), 1, median)
  }

  names(res) <- sub("h", "r", colnames(object$latent))
  class(res) <- "svresid"
  attr(res, "type") <- type

  # Also return posterior mean/median of df parameter if terr = TRUE
  if ("nu" %in% colnames(object$para)) attr(res, "nu") <- get(type)(object$para[,"nu"])

  res
}

#' @export
residuals.svldraws <- function (object, type = "mean", ...) {
  res <- residuals.svdraws(object = object, type = type, ...)
  class(res) <- "svlresid"
  res
}


#' Prediction of Future Log-Volatilities
#' 
#' Simulates draws from the predictive density of the latent log-volatility
#' process.
#' 
#' 
#' @param object \code{svdraws} object.
#' @param steps single number, coercible to integer. Denotes the number of
#' steps to forecast.
#' @param ...  currently ignored.
#' @return Returns an object of class \code{c("svpredict", "mcmc")} containing
#' simulations from the predictive density of \code{h_(n+1),...,h_(n+steps)}.
#' @note You can use the usual \code{coda} methods for \code{mcmc} objects to
#' print, plot, or summarize the predictions, or use them within
#' \code{\link{volplot}} or \code{\link{plot.svdraws}}.
#' @author Gregor Kastner \email{gregor.kastner@@wu.ac.at}
#' @seealso \code{\link{plot.svdraws}}, \code{\link{volplot}}.
#' @keywords ts
#' @examples
#' 
#' ## Simulate a short and highly persistent SV process 
#' sim <- svsim(100, mu = -10, phi = 0.99, sigma = 0.2)
#' 
#' ## Obtain 5000 draws from the sampler (that's not a lot)
#' draws <- svsample(sim$y, draws = 5000, burnin = 100,
#' 		  priormu = c(-10, 1), priorphi = c(20, 1.5), priorsigma = 0.2)
#' 
#' ## Predict 10 days ahead
#' fore <- predict(draws, 10)
#' 
#' ## Check out the results
#' summary(fore)
#' plot(draws, forecast = fore)
#' @export
predict.svdraws <- function(object, steps = 1L, newdata = NULL, ...) {
  if (!(inherits(object, "svdraws"))) stop("Argument 'object' must be of class 'svdraws' or 'svldraws'.")

  steps <- as.integer(steps)
  if (steps < 1) stop("Argument 'steps' must be greater or equal to 1.")
  # Error checking for the mean model
  arorder <- 0  # AR(0) means either a constant mean or a no mean model here
  if (object$meanmodel == "none") {
    if (!is.null(newdata)) warning("No regression coefficients were used when estimating the model. Omitting 'newdata'.")
    regressors <- function (y, newdata, stepind) matrix(0)
  } else if (object$meanmodel == "matrix") {
    if (is.null(newdata)) stop("Regressors needed for prediction. Please provide regressors through parameter 'newdata'.")
    newdata <- as.matrix(newdata)
    if (is.null(object$beta) || NCOL(object$beta) != NCOL(newdata)) stop(paste0("The number of fitted regression coefficients (", NCOL(object$beta), ") does not equal the number of given regressors (", NCOL(newdata), ")."))
    if (NROW(newdata) != steps) stop("The size of the design matrix does not match the number of steps to predict.")
    regressors <- function (y, newdata, stepind) matrix(newdata[stepind, ], nrow = len, ncol = NCOL(newdata), byrow = TRUE)  # matches the format in the ar* case
  } else if (object$meanmodel == "constant") {
    if (!is.null(newdata)) warning("Constant mean was assumed when estimating the model. Omitting 'newdata'.")
    regressors <- function (y, newdata, stepind) matrix(1)
  } else if (any(grep("ar[1-9][0-9]*", object$meanmodel))) {
    arorder <- as.integer(gsub("ar", "", object$meanmodel))
    if (!is.null(newdata)) warning(paste0("An AR(", arorder, ") mean was assumed when estimating the model. Omitting 'newdata'."))
    regressors <- function (y, newdata, stepind) cbind(1, y[,seq.int(stepind, stepind-1+arorder),drop=FALSE])
  } else {
    stop("Unknown mean model. Please contact the developer.")
  }

  thinlatent <- object$thinning$latent
  thinpara <- object$thinning$para
  if (thinpara != thinlatent) {
    warning("Thinning of parameters is different from thinning of latent variables. Trying to sort this out.")
    if (thinpara %% thinlatent == 0) {
      usepara <- 1:(dim(object$para)[1])
      uselatent <- seq(thinpara/thinlatent, dim(object$latent)[1], by=thinpara/thinlatent)
    } else if (thinlatent %% thinpara == 0) {
      uselatent <- 1:(dim(object$latent)[1])
      usepara <- seq(thinlatent/thinpara, dim(object$para)[1], by=thinlatent/thinpara)
    } else stop("Incompatible thinning parameters. Prediction currently not implemented.")
  } else {
    usepara <- uselatent <- seq.int(dim(object$para)[1])
  }

  mu <- object$para[,"mu"][usepara]
  phi <- object$para[,"phi"][usepara]
  sigma <- object$para[,"sigma"][usepara]
  hlast <- object$latent[,dim(object$latent)[2]][uselatent]
  ylast <- tail(object$y, 1)

  mythin <- max(thinpara, thinlatent)
  len <- length(sigma)
  volpred <- matrix(as.numeric(NA), nrow=len, ncol=steps)
  ypred <- matrix(as.numeric(NA), nrow=len, ncol=steps+arorder)

  ypred[,seq_len(arorder)] <- matrix(tail(object$y, arorder), nrow=len, ncol=arorder, byrow=TRUE)  # AR(x) helper columns

  rho <- if ("rho" %in% colnames(object$para)) object$para[, "rho"][usepara] else 0
  nu <- if ("nu" %in% colnames(object$para)) object$para[, "nu"][usepara] else Inf
  betacoeff <- if (exists("beta", object)) object$beta[usepara, c(1, rev(seq_len(NCOL(object$beta)-1))+1), drop=FALSE] else matrix(0)

  resilast <- if (object$meanmodel == "none") {  # last fitted residual
    ylast*exp(-hlast/2)
  } else {  # if mean regression
    (ylast - colSums(object$priors$designmatrix[NROW(object$priors$designmatrix),]*t(betacoeff)))*exp(-hlast/2)  # recycles the last row of the design matrix
  }

  volpred[,1] <- mu+phi*(hlast-mu) + sigma*(rho*resilast + sqrt(1-rho^2)*rnorm(len))
  if (steps > 1) {
    resi <- rt(len, df=nu)  # either rho == 0 or nu == Inf
    incr <- rho*resi + sqrt(1-rho^2)*rnorm(len)
    for (i in seq.int(from=2, to=steps)) {
      ypred[,i-1+arorder] <- rowSums(regressors(ypred, newdata, i-1) * betacoeff) + exp(-volpred[,i-1]/2)*resi
      volpred[,i] <- mu + phi*(volpred[,i-1] - mu) + sigma*incr
    }
  }
  ypred[,steps+arorder] <- rowSums(regressors(ypred, newdata, steps) * betacoeff) + exp(-volpred[,steps]/2)*rnorm(len)

  ypred <- ypred[, setdiff(seq_len(NCOL(ypred)), seq_len(arorder)), drop=FALSE]  # remove temporary AR(x) helper columns
  lastname <- tail(colnames(object$latent), 1)
  lastnumber <- as.integer(gsub("h_", "", lastname))
  colnames(volpred) <- paste0("h_", seq(lastnumber + 1, lastnumber + steps))
  colnames(ypred) <- paste0("y_", seq(lastnumber + 1, lastnumber + steps))
  ret <- list(h = coda::mcmc(volpred, start=mythin, end=len*mythin, thin=mythin),
              y = coda::mcmc(ypred, start=mythin, end=len*mythin, thin=mythin))

  class(ret) <- c("svpredict")
  ret
}

#' @export
predict.svldraws <- function (object, steps = 1L, newdata = NULL, ...) {
  ret <- predict.svdraws(object = object, steps = steps, newdata = newdata, ...)
  class(ret) <- c("svlpredict", "svpredict")
  ret
}

# used to forecast AR-SV models (needs more testing!)


#' Dynamic prediction for the AR-SV model
#' 
#' Simulates draws from the posterior predictive density of a fitted AR-SV
#' model.
#' 
#' 
#' @param object \code{svdraws} object as returned from \code{\link{svsample}}.
#' @param volpred \code{svpredict} object as returned from
#' \code{\link{predict.svdraws}}.
#' @return Returns an object of class \code{c("distpredict", "mcmc")}
#' containing simulations from the posterior predictive density of
#' \code{y_(n+1),...,y_(n+steps)}.
#' @note You can use the usual \code{coda} methods for \code{mcmc} objects to
#' print, plot, or summarize the predictions.
#' @author Gregor Kastner \email{gregor.kastner@@wu.ac.at}
#' @seealso \code{\link{predict.svdraws}}.
#' @keywords ts
#' @examples
#' \dontrun{
#' data(exrates)
#' y <- exrates$USD
#' 
#' ## Fit AR(1)-SV model to EUR-USD exchange rates
#' res <- svsample(y, designmatrix = "ar1")
#' 
#' ## Use predict.svdraws to obtain predictive volatilities
#' ahead <- 100
#' predvol <- predict(res, steps = ahead)
#' 
#' ## Use arpredict to obtain draws from the posterior predictive
#' preddraws <- arpredict(res, predvol)
#' 
#' ## Calculate predictive quantiles
#' predquants <- apply(preddraws, 2, quantile, c(.1, .5, .9))
#' 
#' ## Visualize
#' ts.plot(y, xlim = c(length(y) - ahead, length(y) + ahead),
#' 	ylim = range(predquants))
#' for (i in 1:3) {
#'  lines((length(y) + 1):(length(y) + ahead), predquants[i,],
#'        col = 3, lty = c(2, 1, 2)[i])
#' }
#' }
#' @export arpredict
arpredict <- function(object, volpred) {
 warning("Function `arpredict` has been deprecated. Please use the new features of `predict` instead.")
 if (!inherits(object, "svdraws")) stop("Argument 'object' must be of class 'svdraws'.")
 if (!inherits(volpred, "svpredict")) stop("Argument 'volpred' must be of class 'svpredict'.")
 if (colnames(object$priors$designmatrix)[1] == "const") dynamic <- TRUE else stop("Probably not an AR-specification.")
 order <- ncol(object$priors$designmatrix) - 1

 len <- nrow(volpred)
 steps <- ncol(volpred)
 fromtoby <- attr(volpred, "mcpar")
 usepara <- seq(from = fromtoby[1], to = fromtoby[2], by = fromtoby[3])

 if (ncol(object$para) == 4) {
  nu <- object$para[,"nu"][usepara]
 } else {
  nu <- Inf  # corresponds to conditional normality
 }

 if (ncol(object$beta) > 1) {
  betarev <- object$beta[usepara,c(1,ncol(object$beta):2)]
  lastX <- matrix(c(1, object$y[length(object$y) - order:1 + 1]), nrow = 1)
 } else {
  betarev <- object$beta[usepara,,drop=FALSE]
  lastX <- matrix(1, nrow = 1)
 }

 meanpred <- mcmc(matrix(as.numeric(NA), nrow = len, ncol = steps),
		  start = fromtoby[1], end = fromtoby[2], thin = fromtoby[3])
 meanpred[,1] <- tcrossprod(lastX, betarev) + exp(volpred[,1]/2)*rt(len, df = nu)
 if (steps > 1) {
  lastX <- matrix(rep(lastX, len), nrow = len, byrow = TRUE)
  for (i in (seq.int(steps-1) + 1)) {
   if (ncol(object$beta) > 1)
    lastX <- cbind(lastX[,-2, drop = FALSE], meanpred[,i-1])
   meanpred[,i] <- rowSums(lastX*betarev) + exp(volpred[,i]/2)*rt(len, df = nu)
  }
 }
 class(meanpred) <- c("distpredict", "mcmc")
 meanpred
}
