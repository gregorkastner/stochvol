#' Probability Density Function Plot for the Parameter Posteriors
#' 
#' Displays a plot of the density estimate for the posterior distribution of
#' the parameters \code{mu}, \code{phi}, \code{sigma} (and potentially
#' \code{nu}), computed by the \code{\link[stats]{density}} function.
#' 
#' \code{paradensplot} is modeled after \code{\link[coda]{densplot}} in the
#' \code{coda} package, with some modifications for parameters that have
#' (half-)bounded support.
#' 
#' @param x \code{svdraws} object.
#' @param showobs logical value, indicating whether the observations should be
#' displayed along the x-axis. If many draws have been obtained, the default
#' (\code{TRUE}) can render the plotting to be quite slow, and you might want
#' to try setting \code{showobs} to \code{FALSE}.
#' @param showprior logical value, indicating whether the prior distribution
#' should be displayed. The default value is \code{TRUE}.
#' @param showxlab logical value, indicating whether the x-axis should be
#' labelled with the number of iterations and the bandwith obtained from
#' \code{\link[stats]{density}}. The default value is \code{TRUE}.
#' @param mar numerical vector of length 4, indicating the plot margins. See
#' \code{\link[graphics]{par}} for details. The default value is \code{c(1.9,
#' 1.9, 1.9, 0.5)}, which is slightly smaller than the R-defaults.
#' @param mgp numerical vector of length 3, indicating the axis and label
#' positions. See \code{\link[graphics]{par}} for details. The default value is
#' \code{c(2, 0.6, 0)}, which is slightly smaller than the R-defaults.
#' @param simobj object of class \code{svsim} as returned by the SV simulation
#' function \code{\link{svsim}}. If provided, ``true'' data generating values
#' will be added to the plots.
#' @param \dots further arguments are passed on to the invoked \code{plot}
#' function.
#' @return Called for its side effects. Returns argument \code{x} invisibly.
#' @note You can call this function directly, but it is more commonly called by
#' the \code{\link{plot.svdraws}} method.
#' @author Gregor Kastner \email{gregor.kastner@@wu.ac.at}
#' @seealso \code{\link{paratraceplot}}, \code{\link{volplot}},
#' \code{\link{plot.svdraws}}
#' @keywords hplot
#' @export
paradensplot <- function (x, ...) {
  UseMethod("paradensplot")
}

paradensplot.svdraws <- function(x, showobs = TRUE, showprior = TRUE, showxlab = TRUE,
                                 mar = c(1.9, 1.9, 1.9, .5), mgp = c(2, .6, 0), simobj = NULL, ...) {
  if (!is(x, "svdraws")) stop("This function expects an 'svdraws' object.")
  if (!is.logical(showobs)) stop("If provided, argument 'showobs' must be TRUE or FALSE.")
  if (!is.null(simobj)) {
    if (!is(simobj, "svsim")) stop("If provided, simobj must be an 'svsim' object.")
    sim <- TRUE
  } else sim <- FALSE
  oldpar <- par(mar=mar)
  paranames <- c(quote(mu), quote(phi), quote(sigma), quote(nu))
  cutat1 <- c(FALSE, TRUE, FALSE, FALSE)
  for (i in 1:ncol(x$para)) {
    mydensplot(x$para[,i], show.obs=showobs, main=paste("Density of", paranames[i]),
               cutat1=cutat1[i], showxlab=showxlab, mgp = mgp, ...)
    if (isTRUE(showprior)) {
      paras <- x$priors[[i]]
      vals <- seq(from=par('usr')[1], to=par('usr')[2], len=1000)
      if (i == 1) lines(vals, dnorm(vals, paras[1], paras[2]), col=8, lty=2)
      else if (i == 2) lines(vals, .5*dbeta((vals+1)/2, paras[1], paras[2]), col=8, lty=2)
      else if (i == 3) lines(vals, 2*dnorm(vals, 0, sqrt(paras[1])), col=8, lty=2)
      else if (i == 4) lines(vals, dunif(vals, x$priors$nu[1], x$priors$nu[2]), col = 8, lty = 2)
      if (sim && (i <= 3 || length(simobj$para) == 4)) {
        points(simobj$para[i], 0, col = 3, cex = 2, pch = 16)
      }
    }
  }
  par(oldpar)
  invisible(x)
}

paradensplot.svldraws <- function(x, showobs = TRUE, showprior = TRUE, showxlab = TRUE,
                                  mar = c(1.9, 1.9, 1.9, .5), mgp = c(2, .6, 0), simobj = NULL, ...) {
  if (!is(x, "svldraws")) stop("This function expects an 'svldraws' object.")
  if (!is.logical(showobs)) stop("If provided, argument 'showobs' must be TRUE or FALSE.")
  if (!is.null(simobj)) {
    if (!inherits(simobj, "svsim")) stop("If provided, simobj must be an 'svsim' object.")
    sim <- TRUE
  } else sim <- FALSE
  oldpar <- par(mar=mar)
  paranames <- c(quote(mu), quote(phi), quote(sigma), quote(rho))
  cutat1 <- c(FALSE, TRUE, FALSE, FALSE)
  for (i in 1:ncol(x$para)) {
    mydensplot(x$para[,i], show.obs=showobs, main=paste("Density of", paranames[i]),
               cutat1=cutat1[i], showxlab=showxlab, mgp = mgp, ...)
    if (isTRUE(showprior)) {
      paras <- x$priors[[i]]
      vals <- seq(from=par('usr')[1], to=par('usr')[2], len=1000)
      if (i == 1) lines(vals, dnorm(vals, paras[1], paras[2]), col=8, lty=2)
      else if (i == 2) lines(vals, .5*dbeta((vals+1)/2, paras[1], paras[2]), col=8, lty=2)
      else if (i == 3) lines(vals, 2*dnorm(vals, 0, sqrt(paras[1])), col=8, lty=2)
      else if (i == 4) lines(vals, .5*dbeta((vals+1)/2, paras[1], paras[2]), col=8, lty=2)
      if (sim) {
        points(simobj$para[i], 0, col = 3, cex = 2, pch = 16)
      }
    }
  }
  par(oldpar)
  invisible(x)
}



#' Trace Plot of MCMC Draws from the Parameter Posteriors
#' 
#' Displays a plot of iterations vs. sampled values the parameters \code{mu},
#' \code{phi}, \code{sigma} (and potentially \code{nu}), with a separate plot
#' per variable.
#' 
#' \code{paratraceplot} is modeled after \code{\link[coda]{traceplot}} in the
#' \code{coda} package, with very minor modifications.
#' 
#' @param x \code{svdraws} object.
#' @param mar numerical vector of length 4, indicating the plot margins. See
#' \code{\link[graphics]{par}} for details. The default value is \code{c(1.9,
#' 1.9, 1.9, 0.5)}, which is slightly smaller than the R-defaults.
#' @param mgp numerical vector of length 3, indicating the axis and label
#' positions. See \code{\link[graphics]{par}} for details. The default value is
#' \code{c(2, 0.6, 0)}, which is slightly smaller than the R-defaults.
#' @param simobj object of class \code{svsim} as returned by the SV simulation
#' function \code{\link{svsim}}. If provided, ``true'' data generating values
#' will be added to the plots.
#' @param \dots further arguments are passed on to the invoked \code{matplot}
#' function.
#' @return Called for its side effects. Returns argument \code{x} invisibly.
#' @note You can call this function directly, but it is more commonly called by
#' the \code{\link{plot.svdraws}} method.
#' @author Gregor Kastner \email{gregor.kastner@@wu.ac.at}
#' @seealso \code{\link{paradensplot}}, \code{\link{volplot}},
#' \code{\link{plot.svdraws}}
#' @keywords hplot
#' @export
paratraceplot <- function (x, ...) {
  UseMethod("paratraceplot")
}

paratraceplot.svdraws <- function(x, mar = c(1.9, 1.9, 1.9, .5), mgp = c(2, .6, 0), simobj = NULL, ...) {
  if (!is(x, "svdraws")) stop("This function expects an 'svdraws' object.")
  if (!is.null(simobj)) {
    if (!is(simobj, "svsim")) stop("If provided, simobj must be an 'svsim' object.")
    sim <- TRUE
  } else sim <- FALSE
  oldpar <- par(mar=mar)
  paranames <- c(quote(mu), quote(phi), quote(sigma), quote(nu))
  for (i in 1:ncol(x$para)) {
    mytraceplot(x$para[,i], xlab="", mgp = mgp,
                main=paste("Trace of ", paranames[i], " (thinning = ", x$thinning$para,")", sep=''), ...)
    if (sim && (i <= 3 || length(simobj$para) == 4)) {
      abline(h = simobj$para[i], col = 3, lty = 2)
    }
  }
  par(oldpar)
  invisible(x)
}

paratraceplot.svldraws <- function(x, mar = c(1.9, 1.9, 1.9, .5), mgp = c(2, .6, 0), simobj = NULL, ...) {
  if (!is(x, "svldraws")) stop("This function expects an 'svldraws' object.")
  if (!is.null(simobj)) {
    if (!is(simobj, "svlsim")) stop("If provided, simobj must be an 'svlsim' object.")
    sim <- TRUE
  } else sim <- FALSE
  oldpar <- par(mar=mar)
  paranames <- c(quote(mu), quote(phi), quote(sigma), quote(rho))
  for (i in 1:ncol(x$para)) {
    mytraceplot(x$para[,i], xlab="", mgp = mgp,
                main=paste("Trace of ", paranames[i], " (thinning = ", x$thinning$para,")", sep=''), ...)
    if (sim) {
      abline(h = simobj$para[i], col = 3, lty = 2)
    }
  }
  par(oldpar)
  invisible(x)
}



#' Plotting Quantiles of the Latent Volatilities
#' 
#' Displays quantiles of the posterior distribution of the volatilities over
#' time as well as predictive distributions of future volatilities.
#' 
#' 
#' @param x \code{svdraws} object.
#' @param forecast nonnegative integer or object of class \code{svpredict}, as
#' returned by \code{\link{predict.svdraws}}. If an integer greater than 0 is
#' provided, \code{\link{predict.svdraws}} is invoked to obtain the
#' \code{forecast}-step-ahead prediction. The default value is \code{0}.
#' @param dates vector of length \code{ncol(x$latent)}, providing optional
#' dates for labeling the x-axis. The default value is \code{NULL}; in this
#' case, the axis will be labeled with numbers.
#' @param show0 logical value, indicating whether the initial volatility
#' \code{exp(h_0/2)} should be displayed. The default value is \code{FALSE}.
#' @param col vector of color values (see \code{\link[graphics]{par}}) used for
#' plotting the quantiles. The default value \code{NULL} results in gray lines
#' for all quantiles expect the median, which is displayed in black.
#' @param forecastlty vector of line type values (see
#' \code{\link[graphics]{par}}) used for plotting quantiles of predictive
#' distributions. The default value \code{NULL} results in dashed lines.
#' @param tcl The length of tick marks as a fraction of the height of a line of
#' text. See \code{\link[graphics]{par}} for details. The default value is
#' \code{-0.4}, which results in slightly shorter tick marks than usual.
#' @param mar numerical vector of length 4, indicating the plot margins. See
#' \code{\link[graphics]{par}} for details. The default value is \code{c(1.9,
#' 1.9, 1.9, 0.5)}, which is slightly smaller than the R-defaults.
#' @param mgp numerical vector of length 3, indicating the axis and label
#' positions. See \code{\link[graphics]{par}} for details. The default value is
#' \code{c(2, 0.6, 0)}, which is slightly smaller than the R-defaults.
#' @param simobj object of class \code{svsim} as returned by the SV simulation
#' function \code{\link{svsim}}. If provided, ``true'' data generating values
#' will be added to the plot(s).
#' @param \dots further arguments are passed on to the invoked \code{ts.plot}
#' function.
#' @return Called for its side effects. Returns argument \code{x} invisibly.
#' @note In case you want different quantiles to be plotted, use
#' \code{\link{updatesummary}} on the \code{svdraws} object first. An example
#' of doing so is given below.
#' @author Gregor Kastner \email{gregor.kastner@@wu.ac.at}
#' @seealso \code{\link{updatesummary}}, \code{\link{paratraceplot}},
#' \code{\link{paradensplot}}, \code{\link{plot.svdraws}}.
#' @keywords hplot ts
#' @examples
#' 
#' ## Simulate a short and highly persistent SV process 
#' sim <- svsim(100, mu = -10, phi = 0.99, sigma = 0.2)
#' 
#' ## Obtain 5000 draws from the sampler (that's not a lot)
#' draws <- svsample(sim$y, draws = 5000, burnin = 100,
#' 		  priormu = c(-10, 1), priorphi = c(20, 1.5),
#' 		  priorsigma = 0.2)
#' 
#' ## Plot the latent volatilities and some forecasts
#' volplot(draws, forecast = 10)
#' 
#' ## Re-plot with different quantiles
#' newquants <- c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
#' draws <- updatesummary(draws, quantiles=newquants)
#' 
#' volplot(draws, forecast = 10)
#' 
#' @export
volplot <- function(x, forecast = 0, dates = NULL, show0 = FALSE,
                    col = NULL, forecastlty = NULL, tcl = -.4,
                    mar = c(1.9, 1.9, 1.9, .5), mgp = c(2, .6, 0), simobj = NULL, ...) {
  if (!is(x, "svdraws") && !is(x, "svldraws")) stop("This function expects an 'svdraws' or an 'svldraws' object.")
  if (!is.null(simobj)) {
    if (!is(simobj, "svsim") && !is(simobj, "svlsim")) stop("If provided, simobj must be an 'svsim' or an 'svlsim' object.")
    sim <- TRUE
  } else sim <- FALSE
  oldpar <- par(mar = mar)
  where <- grep("%", dimnames(x$summary$latent)[[2]])
  volquants <- t(100*exp(x$summary$latent[,where,drop=FALSE]/2))  # monotone transformation!
  nvolquants <- dim(volquants)[1]
  timelen <- dim(volquants)[2]
  if (is.null(nvolquants) | all(is.na(volquants))) stop("No quantiles to plot.")
  if (is.null(col)) {
    cols <- rep(8, nvolquants)
    cols[dimnames(volquants)[[1]] == "50%"] <- 1
  } else cols <- col
  if (is.null(forecastlty)) forecastlty <- 2

  if (is(forecast, "svpredict") || is(forecast, "svlpredict") || (is.numeric(forecast) && length(forecast) == 1 && all(forecast != 0))) { # also draw future values
    thintime <- x$thinning$time

    if (thintime != 1) {
      lasth <- as.integer(gsub("h_", "", dimnames(x$latent)[[2]][dim(x$latent)[2]]))
      if (length(x$y) > lasth) {
        warning(paste("Thinning for time 'thintime' has not been set to one during sampling. This means we are forecasting conditional on h_", lasth, " and not on h_", length(x$y), ".", sep=''))
      }
    }

    if(is.numeric(forecast) && length(forecast) == 1 && all(forecast >= 1)) {
      #   warning("Calling prediction method.")
      forecast <- predict(x, forecast)
    }
    if(!is(forecast, "svpredict") || !is(forecast, "svlpredict")) stop("Argument 'forecast' must be a single nonnegative integer, or of class type 'svpredict' or 'svlpredict'.")
    if(is(forecast, "svlpredict")) {
      show0 <- FALSE
    }

    futlen <- dim(forecast)[2]

    xs <- matrix(rep(seq(timelen, timelen + futlen/thintime, len=futlen+1), nvolquants), nrow=futlen+1)
    quants <- as.numeric(gsub("%", "", dimnames(volquants)[[1]]))/100
    ys <- rbind(volquants[,timelen], t(matrix(apply(100*exp(forecast/2), 2, quantile, quants), nrow=nvolquants)))

    if (futlen/thintime > .01*timelen) {  # increase xlim to give space for forecast
      if (thintime == 1) {
        xlim <- c(0, timelen+futlen/thintime)
      } else {
        xlim <- c(1, timelen+futlen/thintime)
      }
    } else {
      xlim <- NULL
    }
  } else xlim <- NULL

  if(exists("sd", x$summary)) {
    mymain <- paste("Estimated conditional volatilities in percent (", paste(dimnames(volquants)[[1]], collapse=' / '),
                    " posterior quantiles)", sep = '')
  } else {
    mymain <- paste("Estimated volatilities in percent (", paste(dimnames(volquants)[[1]], collapse=' / '),
                    " posterior quantiles)", sep = '')
  }

  ts.plot(t(volquants), gpars=list(xlim=xlim, col=cols, xlab='', xaxt='n', mgp=mgp, tcl=tcl,
                                   main = mymain, ...))

  if (sim) {
    lines(100*simobj$vol, col = 3)
  }

  if (is(forecast, "svpredict") || is(forecast, "svlpredict")) {
    for (i in 1:nvolquants) lines(xs[,i], ys[,i], lty=forecastlty, col=cols[i])
  }

  ax <- axis(1, tick=FALSE, labels=FALSE)  # just automagic axis ticks, don't draw yet

  if (show0) { # also draw latent0:
    thintime <- x$thin$time
    xs <- matrix(rep(c(1-1/thintime,1), nvolquants), nrow=2)
    where <- grep("%", names(x$summary$latent0))
    ys <- rbind(100*exp(x$summary$latent0[where]/2), volquants[,1])
    for (i in 1:nvolquants) lines(xs[,i], ys[,i], lty=forecastlty, col=cols[i])
  }

  if (is.null(dates)) {
    dates <- c(0L, as.integer(gsub("h_", "", dimnames(x$latent)[[2]])))
    if (max(ax) > length(dates)) {  # means we are probably forecasting and need extra axis labels
      dates <- c(dates, seq(length(dates), max(ax), by=dates[2]-dates[1]))
    }
  } else {
    if (is(dates, "Date")) dates <- as.character(dates)
    if (length(dates) != ncol(x$latent)) {
      stop("Length of argument 'dates' differs from ncol(x$latent).")
    }
    dates <- c('', dates)
    ax <- ax[ax != 0]  # avoid "zero" tick
  }
  axis(1, at=ax, labels=dates[ax+1], mgp=mgp, tcl=tcl)

  if(exists("sd", x$summary)) {  # only for t distributed residuals
    where <- grep("%", dimnames(x$summary$latent)[[2]])
    ts.plot(100*x$summary$sd[,where], gpars=list(xlim=xlim, col=cols, xlab='', xaxt='n', mgp=mgp, tcl=tcl,
                                                 main = paste("Estimated volatilities in percent (",
                                                              paste(dimnames(volquants)[[1]], collapse=' / '),
                                                              " posterior quantiles)", sep=''), ...))

    if (sim) {
      standardizer <- sqrt(simobj$para$nu / (simobj$para$nu - 2))
      lines(100*simobj$vol*standardizer, col = 3)
    }

    if (is(forecast, "svpredict")) {
      standardizer <- sqrt(x$para[,"nu"] / (x$para[,"nu"] - 2))
      where <- grep("%", dimnames(x$summary$latent)[[2]])
      ys <- rbind(100*x$summary$sd[timelen,where,drop=FALSE],
                  t(matrix(apply(100*exp(forecast/2)*standardizer, 2, quantile, quants), nrow=nvolquants)))

      for (i in 1:nvolquants) lines(xs[,i], ys[,i], lty=forecastlty, col=cols[i])
    }
    axis(1, at=ax, labels=dates[ax+1], mgp=mgp, tcl=tcl)
  }

  par(oldpar)
  invisible(x)
}



#' Graphical Summary of the Posterior Distribution
#' 
#' \code{plot.svdraws} generates some plots visualizing the posterior
#' distribution and can also be used to display predictive distributions of
#' future volatilities.
#' 
#' This function sets up the page layout and calls \code{\link{volplot}},
#' \code{\link{paratraceplot}} and \code{\link{paradensplot}}.
#' 
#' @param x \code{svdraws} object.
#' @param forecast nonnegative integer or object of class \code{svpredict}, as
#' returned by \code{\link{predict.svdraws}}. If an integer greater than 0 is
#' provided, \code{\link{predict.svdraws}} is invoked to obtain the
#' \code{forecast}-step-ahead prediction. The default value is \code{0}.
#' @param dates vector of length \code{ncol(x$latent)}, providing optional
#' dates for labelling the x-axis. The default value is \code{NULL}; in this
#' case, the axis will be labelled with numbers.
#' @param show0 logical value, indicating whether the initial volatility
#' \code{exp(h_0/2)} should be displayed. The default value is \code{FALSE}.
#' @param showobs logical value, indicating whether the observations should be
#' displayed along the x-axis. If many draws have been obtained, the default
#' (\code{TRUE}) can render the plotting to be quite slow, and you might want
#' to try setting \code{showobs} to \code{FALSE}.
#' @param showprior logical value, indicating whether the prior distribution
#' should be displayed. The default value is \code{TRUE}.
#' @param col vector of color values (see \code{\link[graphics]{par}}) used for
#' plotting the quantiles. The default value \code{NULL} results in gray lines
#' for all quantiles expect the median, which is displayed in black.
#' @param forecastlty vector of line type values (see
#' \code{\link[graphics]{par}}) used for plotting quantiles of predictive
#' distributions. The default value \code{NULL} results in dashed lines.
#' @param tcl The length of tick marks as a fraction of the height of a line of
#' text. See \code{\link[graphics]{par}} for details. The default value is
#' \code{-0.4}, which results in slightly shorter tick marks than usual.
#' @param mar numerical vector of length 4, indicating the plot margins. See
#' \code{\link[graphics]{par}} for details. The default value is \code{c(1.9,
#' 1.9, 1.9, 0.5)}, which is slightly smaller than the R-defaults.
#' @param mgp numerical vector of length 3, indicating the axis and label
#' positions. See \code{\link[graphics]{par}} for details. The default value is
#' \code{c(2, 0.6, 0)}, which is slightly smaller than the R-defaults.
#' @param simobj object of class \code{svsim} as returned by the SV simulation
#' function \code{\link{svsim}}. If provided, the ``true'' data generating
#' values will be added to the plots.
#' @param \dots further arguments are passed on to the invoked plotting
#' functions.
#' @return Called for its side effects. Returns argument \code{x} invisibly.
#' @note In case you want different quantiles to be plotted, use
#' \code{\link{updatesummary}} on the \code{svdraws} object first. An example
#' of doing so is given in the Examples section.
#' @author Gregor Kastner \email{gregor.kastner@@wu.ac.at}
#' @seealso \code{\link{updatesummary}}, \code{\link{volplot}},
#' \code{\link{paratraceplot}}, \code{\link{paradensplot}}.
#' @keywords hplot
#' @examples
#' 
#' ## Simulate a short and highly persistent SV process 
#' sim <- svsim(100, mu = -10, phi = 0.99, sigma = 0.2)
#' 
#' ## Obtain 5000 draws from the sampler (that's not a lot)
#' draws <- svsample(sim$y, draws = 5000, burnin = 1000,
#'   priormu = c(-10, 1), priorphi = c(20, 1.5), priorsigma = 0.2)
#' 
#' ## Plot the latent volatilities and some forecasts
#' plot(draws, forecast = 10)
#' 
#' ## Re-plot with different quantiles
#' newquants <- c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
#' draws <- updatesummary(draws, quantiles = newquants)
#' 
#' plot(draws, forecast = 20, showobs = FALSE, col = seq(along = newquants),
#'      forecastlty = 3, showprior = FALSE)
#' 
#' @export
plot.svdraws <- function(x, forecast = NULL, dates = NULL,
			 show0 = FALSE, showobs = TRUE, showprior = TRUE, col = NULL,
			 forecastlty = NULL, tcl = -0.4,
			 mar = c(1.9, 1.9, 1.7, .5), mgp = c(2, .6, 0),
			 simobj = NULL, ...) {
  oldpar <- par(mfrow=c(1,1))
  if (ncol(x$para) == 4) {
    layout(matrix(c(1, 1, 1, 1, 2, 2, 2, 2, 3, 4, 5, 6, 7, 8, 9, 10), 4, byrow = TRUE))
  } else {
    layout(matrix(c(1, 1, 1, 2, 3, 4, 5, 6, 7), 3, byrow = TRUE))
  }
  volplot(x, dates = dates, show0 = show0, forecast = forecast,
          forecastlty = forecastlty, col = col, tcl = tcl, mar = mar,
          mgp = mgp, simobj = simobj, ...)
  paratraceplot(x, mar = mar, mgp = mgp, simobj = simobj, ...)
  paradensplot(x, showobs = showobs, showprior = showprior,
               showxlab = FALSE, mar = mar, mgp = mgp, simobj = simobj, ...)
  par(oldpar)
  invisible(x)
}


# TODO joint docs with plot.svdraws

#' @export
plot.svldraws <- function (x, forecast = NULL, dates = NULL,
			 show0 = FALSE, showobs = TRUE, showprior = TRUE, col = NULL,
			 forecastlty = NULL, tcl = -0.4,
			 mar = c(1.9, 1.9, 1.7, .5), mgp = c(2, .6, 0),
			 simobj = NULL, ...) {
  oldpar <- par(mfrow=c(1,1))
  layout(matrix(c(1, 1, 1, 1, 2, 3, 4, 5, 6, 7, 8, 9), ncol=4, byrow = TRUE))
  volplot(x, dates = dates, show0 = show0, forecast = forecast,
          forecastlty = forecastlty, col = col, tcl = tcl, mar = mar,
          mgp = mgp, simobj = simobj, ...)
  paratraceplot(x, mar = mar, mgp = mgp, simobj = simobj, ...)
  paradensplot(x, showobs = showobs, showprior = showprior,
               showxlab = FALSE, mar = mar, mgp = mgp, simobj = simobj, ...)
  par(oldpar)
  invisible(x)
}

# modified density plot (from coda package)
mydensplot <- function(x, show.obs = TRUE, bwf, main = "", ylim, cutat1=FALSE, showxlab=TRUE, mgp = c(2,.6,0), tcl=-.4, ...) {
  for (i in 1:nvar(x)) {
    x.i <- as.matrix(x)[, i, drop = TRUE]
    if (missing(bwf)) 
      bwf <- function(xx) {
        xx <- xx[!is.na(as.vector(xx))]
        return(1.06 * min(sd(xx), IQR(xx)/1.34) * length(xx)^-0.2)
      }
    bw <- bwf(x.i)
    width <- 4 * bw
    if (max(abs(x.i - floor(x.i))) == 0 || bw == 0) {
      hist(x.i, prob = TRUE, main = main, ...)
    } else {
      density.scale <- "open"
      if (isTRUE(cutat1)) {
        if (1-max(x.i) < 2*bw) {
          density.scale <- "cutat1"
          x.i <- c(x.i, 2 - x.i)
          if (1+min(x.i) < 2*bw) {
            density.scale <- "cutatboth"
            x.i <- c(x.i, -2 - x.i, 2 - x.i)
          }
        } else if (1+min(x.i) < 2*bw) {
          density.scale <- "cutat-1"
          x.i <- c(x.i, -2 - x.i)
        }
      } else if (max(x.i) <= 1 && 1 - max(x.i) < 2 * bw) {
        if (min(x.i) >= 0 && min(x.i) < 2 * bw) {
          density.scale <- "proportion"
          x.i <- c(x.i, -x.i, 2 - x.i)
        }
      } else if (min(x.i) >= 0 && min(x.i) < 2 * bw) {
        density.scale <- "positive"
        x.i <- c(x.i, -x.i)
      }
      dens <- density(x.i, width = width)

      if (density.scale == "proportion") {
        dens$y <- 3 * dens$y[dens$x >= 0 & dens$x <= 1]
        dens$x <- dens$x[dens$x >= 0 & dens$x <= 1]
      }
      else if (density.scale == "positive") {
        dens$y <- 2 * dens$y[dens$x >= 0]
        dens$x <- dens$x[dens$x >= 0]
      }
      else if (density.scale == "cutat1") {
        dens$y <- 2 * dens$y[dens$x <= 1]
        dens$x <- dens$x[dens$x <= 1]
      }
      else if (density.scale == "cutat-1") {
        dens$y <- 2 * dens$y[dens$x >= -1]
        dens$x <- dens$x[dens$x >= -1]
      }
      else if (density.scale == "cutatboth") {
        dens$y <- 3 * dens$y[dens$x >= -1 & dens$x <= 1]
        dens$x <- dens$x[dens$x >= -1 & dens$x <= 1]
      }

      if (missing(ylim)) 
        ylim <- c(0, max(dens$y))

      plot(dens, ylab = "", main = main, type = "l", 
           ylim = ylim, xlab="", mgp = mgp, tcl = tcl, ...)
      if(isTRUE(showxlab)) {
        if (is.R()) {
          mtext(paste("N =", niter(x), "  Bandwidth =",
                      formatC(dens$bw)), side=1, line=2.7, cex=.7)
        } else {
          mtext(paste("N =", niter(x), "  Bandwidth =",
                      formatC(bw)), side=1, line=2.7, cex=.7)
        }
      }
      if (show.obs) 
        lines(x.i[1:niter(x)], rep(max(dens$y)/100, niter(x)), 
              type = "h")
    }
    if (!is.null(varnames(x)) && is.null(list(...)$main)) 
      title(paste("Density of", varnames(x)[i]))
  }
  invisible(x)
}

# modified traceplot (from coda)
mytraceplot <- function (x, smooth = FALSE, col = 1:6, type = "l", ylab = "", xlab = "Iterations", mgp = c(2,.6,0), tcl = -.4, ...) {
  x <- mcmc.list(x)
  for (j in 1:nvar(x)) {
    xp <- as.vector(time(x))
    yp <- if (nvar(x) > 1) x[, j, drop = TRUE] else x
    yp <- do.call("cbind", yp)
    matplot(xp, yp, xlab = xlab, ylab = ylab, type = type, 
            col = col, mgp = mgp, tcl = tcl, ...)
    if (!is.null(varnames(x)) && is.null(list(...)$main)) 
      title(paste("Trace of ", varnames(x)[j], " (thin = ", attr(x, "thinning")$thinpara,")", sep=''))
    if (smooth) {
      scol <- rep(col, length = nchain(x))
      for (k in 1:nchain(x))
        lines(lowess(xp, yp[, k]), col = scol[k])
    }
  }
  invisible(x)
}

#' @export
plot.svresid <- function(x, origdata = NA,
                         mains = c("Residual plot", "Q-Q plot"),
                         mar = c(2.9, 2.7, 2.2, .5),
                         mgp = c(1.7, .6, 0), ...) {

  if (any(is.na(origdata))) {
    oldpar <- par(mfrow=c(1, 2), mar=mar, mgp=mgp)
  } else {
    oldpar <- par(mfrow=c(2, 2), mar=mar, mgp=mgp)
    plot.default(origdata, ylab='Original values', xlab='Time', xaxt='n', ylim=c(-1,1)*max(abs(origdata)), main="Original data", ...)
    where <- seq(1, length(origdata), length=min(7, length(origdata)))
    axis(1, at = where, labels = names(origdata)[where])
    qqnorm(origdata, main=paste(mains[2], "for original data"))
    qqline(origdata, probs = c(0.01, 0.99))
  }

  # Cater for conditional t-distributions
  if (!is.null(attr(x, "nu"))) {
    terr <- TRUE
    nu <- attr(x, "nu")
    xlab <- paste("Theoretical quantiles from a t-distribution with", round(nu, 2), "df")
  } else {
    terr <- FALSE
    nu <- Inf
    xlab <- "Theoretical quantiles from a standard normal distribution"
  }

  plot.default(x, ylab = paste("M", substring(attr(x, "type"), 2), ' standardized residuals', sep = ""),
               xlab='Time', xaxt='n', ylim=c(-1,1)*max(abs(x)),
               main=mains[1], ...)

  if (!terr) abline(h=qnorm(c(.025, .975)), lty=2)
  where <- seq(1, length(x), length=min(7, length(x)))
  axis(1, at = where, labels = gsub("r_", "", names(x)[where]))
  qqplot(qt(ppoints(length(x)), df = nu), x,
         main=paste(mains[2], "for", attr(x, "type"), "standardized residuals"),
         xlab = xlab, ylab = "Sample quantiles")
  qqline(x, probs = c(0.01, 0.99), distribution = function(x) qt(x, df = nu))
  par(oldpar)
  invisible(x)
}


#' @export
plot.svlresid <- function(x, origdata = NA,
                          mains = c("Residual plot", "Q-Q plot"),
                          mar = c(2.9, 2.7, 2.2, .5),
                          mgp = c(1.7, .6, 0), ...) {
  plot.svresid(x = x, origdata = origdata,
               mains = mains,
               mar = mar,
               mgp = mgp, ...)
}
