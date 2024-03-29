#  #####################################################################################
#  R package stochvol by
#     Gregor Kastner Copyright (C) 2013-2018
#     Gregor Kastner and Darjus Hosszejni Copyright (C) 2019-
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

#' Graphical Summary of the Posterior Predictive Distribution
#'
#' \code{plot.svpredict} and \code{plot.svlpredict} generate some plots
#' visualizing the posterior predictive distribution of future volatilites and
#' future observations.
#'
#' @param x \code{svpredict} or \code{svlpredict} object.
#' @param quantiles Which quantiles to plot? Defaults to
#' \code{c(.05, .25, .5, .75, .95)}.
#' @param \dots further arguments are passed on to the invoked
#' \code{\link[stats]{ts.plot}} or \code{\link[graphics]{boxplot}} function.
#' @return Called for its side effects. Returns argument \code{x} invisibly.
#' @note Note that \code{svpredict} or \code{svlpredict} objects can also be
#' used within \code{\link{plot.svdraws}} for a possibly more useful
#' visualization. See the examples in \code{\link{predict.svdraws}} and
#' those below for use cases.
#' @family plotting
#' @keywords hplot
#' @examples
#'
#' ## Simulate a short and highly persistent SV process
#' sim <- svsim(100, mu = -10, phi = 0.99, sigma = 0.1)
#'
#' ## Obtain 5000 draws from the sampler (that's not a lot)
#' draws <- svsample(sim$y, draws = 5000, burnin = 1000)
#'
#' ## Predict 10 steps ahead
#' pred <- predict(draws, 10)
#'
#' ## Visualize the predicted distributions
#' plot(pred)
#'
#' ## Plot the latent volatilities and some forecasts
#' plot(draws, forecast = pred)
#'
#' @family plotting
#' @keywords hplot
#' @export
plot.svpredict <- function(x, quantiles = c(.05, .25, .5, .75, .95), ...) {
  nvar <- coda::nvar(x$h[[1]])
  oldpar <- par(mfrow = c(2, 1), mgp = c(1.8, .8, 0), mar = c(3, 3, 3, 1))
  if (nvar == 1L) {
    boxplot(simplify2array(x$vol, higher = FALSE),  # previous solution only worked with length(quantiles) == 5
            names = paste("Next period / chain", seq_along(x$vol)),
            show.names = TRUE, ...)
    title("Predicted volatility")
  } else {
    ts.plot(matrix(aperm(simplify2array(lapply(x$vol, apply, 2, quantile, quantiles)), c(2, 1, 3)), nrow = nvar),
            xlab = "Periods ahead", ...)
    title(paste0("Predicted volatility (", paste0(100*quantiles, collapse = '% / '),
                 "% quantiles)"))
  }
  if (nvar == 1L) {
    boxplot(simplify2array(x$y, higher = FALSE),
            names = paste("Next period / chain", seq_along(x$y)),
            show.names = TRUE, ...)
    title("Predicted volatility")
  } else {
    ts.plot(matrix(aperm(simplify2array(lapply(x$y, apply, 2, quantile, quantiles)), c(2, 1, 3)), nrow = nvar),
            xlab = "Periods ahead", ...)
    title(paste0("Predicted data (", paste0(100*quantiles, collapse = '% / '),
                 "% quantiles)"))
  }
  par(oldpar)
  invisible(x)
}

#' Probability Density Function Plot for the Parameter Posteriors
#'
#' Displays a plot of the density estimate for the posterior distribution of
#' the parameters \code{mu}, \code{phi}, \code{sigma} (and potentially
#' \code{nu} or \code{rho}), computed by the \code{\link[stats]{density}} function.
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
#' @family plotting
#' @keywords hplot
#' @export
paradensplot <- function(x, showobs = TRUE, showprior = TRUE, showxlab = TRUE,
                         mar = c(1.9, 1.9, 1.9, .5), mgp = c(2, .6, 0), simobj = NULL, ...) {
  if (!inherits(x, "svdraws")) stop("This function expects an 'svdraws' object.")
  if (!is.logical(showobs)) stop("If provided, argument 'showobs' must be TRUE or FALSE.")
  sim <- !is.null(simobj)
  if (sim && !inherits(simobj, "svsim")) {
    stop("If provided, simobj must be an 'svsim' object.")
  }
  oldpar <- par(mar = mar)
  para_names <- c(mu = quote(mu), phi = quote(phi), sigma = quote(sigma), rho = quote(rho), nu = quote(nu))
  para_prior_names <- c(mu = "mu", phi = "phi", sigma = "sigma2", rho = "rho", nu = "nu")
  params <- sampled_parameters(x)
  for (para_name in params) {
    prior <- x$priors[[para_prior_names[para_name] ]]
    cutat <- x$para_inv_transform[[para_name]](range(prior))
    mydensplot(x$para[,para_name], show.obs=showobs, main=paste("Density of", para_names[para_name]),
               cutat=cutat,
               showxlab=showxlab, mgp = mgp, ...)
    if (isTRUE(showprior)) {
      eps <- diff(par('usr'))[1] / 999
      vals <- seq(from=max(c(par('usr')[1], cutat[1]+eps)), to=min(c(par('usr')[2], cutat[2]-eps)), len=1000)
      lines(vals,
            x$para_transform_det[[para_name]](vals) * density(prior)(x$para_transform[[para_name]](vals)),
            col = 8, lty = 2)
      if (sim) {
        points(simobj$para[[para_name]], 0, col = 3, cex = 2, pch = 16)
      }
    }
  }
  par(oldpar)
  invisible(x)
}

#' Trace Plot of MCMC Draws from the Parameter Posteriors
#'
#' Displays a plot of iterations vs. sampled values the parameters \code{mu},
#' \code{phi}, \code{sigma} (and potentially \code{nu} or \code{rho}), with a separate plot
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
#' @family plotting
#' @keywords hplot
#' @export
paratraceplot.svdraws <- function(x, mar = c(1.9, 1.9, 1.9, .5), mgp = c(2, .6, 0), simobj = NULL, ...) {
  if (!inherits(x, "svdraws")) stop("This function expects an 'svdraws' object.")
  if (!is.null(simobj)) {
    if (!inherits(simobj, "svsim")) stop("If provided, simobj must be an 'svsim' object.")
    sim <- TRUE
  } else sim <- FALSE
  oldpar <- par(mar=mar)
  paranames <- c(mu=quote(mu), phi=quote(phi), sigma=quote(sigma), nu=quote(nu), rho=quote(rho))
  params <- sampled_parameters(x)
  for (para_name in params) {
    mytraceplot(para(x, "all")[,para_name], xlab="", mgp = mgp,
                main=paste0("Trace of ", paranames[para_name], " (thin = ", x$thinning$para,")"), ...)
    if (sim && para_name %in% names(simobj$para)) {
      abline(h = simobj$para[[para_name]], col = 3, lty = 2)
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
#' Only available for inputs \code{x} of class \code{svdraws}.
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
#' @param \dots further arguments are passed on to the invoked \code{\link[stats]{ts.plot}}
#' function.
#' @param newdata corresponds to parameter \code{newdata} in \code{\link{predict.svdraws}}.
#' \emph{Only if \code{forecast} is a positive integer and \code{\link{predict.svdraws}}
#' needs a \code{newdata} object.} Corresponds to input
#' parameter \code{designmatrix} in \code{\link{svsample}}.
#' A matrix of regressors with number of rows equal to parameter \code{forecast}.
#' @return Called for its side effects. Returns argument \code{x} invisibly.
#' @note In case you want different quantiles to be plotted, use
#' \code{\link{updatesummary}} on the \code{svdraws} object first. An example
#' of doing so is given below.
#' @author Gregor Kastner \email{gregor.kastner@@wu.ac.at}
#' @seealso \code{\link{updatesummary}}, \code{\link{predict.svdraws}}
#' @keywords hplot ts
#' @family plotting
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
#' draws <- updatesummary(draws, quantiles = newquants)
#'
#' volplot(draws, forecast = 10)
#'
#' @export
volplot <- function(x, forecast = 0, dates = NULL, show0 = FALSE,
                    forecastlty = NULL, tcl = -.4,
                    mar = c(1.9, 1.9, 1.9, .5), mgp = c(2, .6, 0), simobj = NULL,
                    newdata = NULL, ...) {
  if (!inherits(x, "svdraws")) stop("This function expects an 'svdraws' object.")
  if (x$thinning$time != "all") stop("This function requires that all volatilities have been stored during sampling.")
  if (!is.null(simobj)) {
    if (!inherits(simobj, "svsim")) stop("If provided, simobj must be an 'svsim' object.")
    sim <- TRUE
  } else sim <- FALSE
  oldpar <- par(mar = mar)
  to_plot <- x$summary$sd
  where <- grep("%", colnames(to_plot))
  volquants <- t(100*to_plot[,where,drop=FALSE])
  nvolquants <- NROW(volquants)
  timelen <- NCOL(volquants)
  if (is.null(nvolquants) || all(is.na(volquants))) stop("No quantiles to plot.")
  alphas <- 1 - 1.2 * abs(0.5 - as.numeric(gsub("%", "", colnames(to_plot)[grep("%", colnames(to_plot))]))/100)
  # set alpha for colors
  cols <- cbind(1, alphas)
  cols <- apply(cols, 1, function (x) {
                  rgb_cols <- as.vector(col2rgb(x[1]))/255
                  rgb(rgb_cols[1], rgb_cols[2], rgb_cols[3], alpha=x[2])})
  if (is.null(forecastlty)) forecastlty <- 2

  if (inherits(forecast, "svpredict") || (is.numeric(forecast) && length(forecast) == 1 && all(forecast != 0))) { # also draw future values
    lasth <- as.integer(gsub("h_", "", tail(coda::varnames(x$latent[[1]]), 1)))
    if (length(x$y) > lasth) {  # should never happen
      stop("The last log variance, h_n, has not been stored during sampling. Aborting.")
    }

    if(is.numeric(forecast) && length(forecast) == 1 && isTRUE(forecast >= 1)) {
      forecast <- predict(x, forecast, newdata)
    }
    if(!inherits(forecast, "svpredict")) stop("Argument 'forecast' must be a single nonnegative integer, or of class type 'svpredict' or 'svlpredict'.")
    if(inherits(forecast, "svlpredict")) {
      if (show0) warning("Initial volatility not available for the SV model with leverage. Setting 'show0' to 'FALSE'.")
      show0 <- FALSE
    }

    volpred <- join_mcmclist(forecast$vol)
    futlen <- NCOL(volpred)

    xs <- matrix(rep(seq(timelen, timelen + futlen, len=futlen+1), nvolquants), nrow=futlen+1)
    quants <- as.numeric(gsub("%", "", rownames(volquants)))/100
    ys <- rbind(volquants[,timelen], t(matrix(apply(100*volpred, 2, quantile, quants), nrow=nvolquants)))

    if (futlen > .01*timelen) {  # increase xlim to give space for forecast
      xlim <- c(0, timelen + futlen)
    } else {
      xlim <- NULL
    }
  } else xlim <- NULL

  ts.plot(t(volquants), gpars=list(xlim=xlim, col=cols, xlab='', xaxt='n', mgp=mgp, tcl=tcl,
                                   main = paste("Estimated volatilities in percent (", paste(dimnames(volquants)[[1]], collapse=' / '),
                                                " posterior quantiles)", sep = ''),
                                   ...))

  if (sim) {
    lines(100*simobj$vol, col = 3)
  }

  if (inherits(forecast, "svpredict")) {
    for (i in seq_len(nvolquants)) lines(xs[,i], ys[,i], lty=forecastlty, col=cols[i])
  }

  ax <- axis(1, tick=FALSE, labels=FALSE)  # just automagic axis ticks, don't draw yet

  if (show0) { # also draw latent0:
    xs <- matrix(rep(c(0, 1), nvolquants), nrow=2)
    where <- grep("%", names(x$summary$latent0))
    ys <- rbind(100*exp(x$summary$latent0[where]/2), volquants[,1])
    for (i in seq_len(nvolquants)) lines(xs[,i], ys[,i], lty=forecastlty, col=cols[i])
  }

  if (is.null(dates)) {
    dates <- c(0L, as.integer(gsub("h_", "", coda::varnames(latent(x)))))
    if (max(ax) > length(dates)) {  # means we are probably forecasting and need extra axis labels
      dates <- c(dates, seq(length(dates), max(ax), by=dates[2]-dates[1]))
    }
  } else {
    if (inherits(dates, "Date")) dates <- as.character(dates)
    if (length(dates) != NCOL(latent(x, 1))) {
      stop("Length of argument 'dates' differs from NCOL(latent(x, 1)).")
    }
    dates <- c('', dates)
    ax <- ax[ax != 0]  # avoid "zero" tick
  }
  axis(1, at=ax, labels=dates[ax+1], mgp=mgp, tcl=tcl)

  par(oldpar)
  invisible(x)
}


#' Graphical Summary of the Posterior Distribution
#'
#' \code{plot.svdraws} and \code{plot.svldraws} generate some plots visualizing the posterior
#' distribution and can also be used to display predictive distributions of
#' future volatilities.
#'
#' These functions set up the page layout and call \code{\link{volplot}},
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
#' Only available for inputs \code{x} of class \code{svdraws}.
#' @param showobs logical value, indicating whether the observations should be
#' displayed along the x-axis. If many draws have been obtained, the default
#' (\code{TRUE}) can render the plotting to be quite slow, and you might want
#' to try setting \code{showobs} to \code{FALSE}.
#' @param showprior logical value, indicating whether the prior distribution
#' should be displayed. The default value is \code{TRUE}.
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
#' @param newdata corresponds to parameter \code{newdata} in \code{\link{predict.svdraws}}.
#' \emph{Only if \code{forecast} is a positive integer and \code{\link{predict.svdraws}}
#' needs a \code{newdata} object.} Corresponds to input
#' parameter \code{designmatrix} in \code{\link{svsample}}.
#' A matrix of regressors with number of rows equal to parameter \code{forecast}.
#' @param \dots further arguments are passed on to the invoked plotting
#' functions.
#' @return Called for its side effects. Returns argument \code{x} invisibly.
#' @note In case you want different quantiles to be plotted, use
#' \code{\link{updatesummary}} on the \code{svdraws} object first. An example
#' of doing so is given in the Examples section.
#' @author Gregor Kastner \email{gregor.kastner@@wu.ac.at}
#' @seealso \code{\link{updatesummary}}, \code{\link{predict.svdraws}}
#' @family plotting
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
#' plot(draws, forecast = 20, showobs = FALSE,
#'      forecastlty = 3, showprior = FALSE)
#'
#' @export
plot.svdraws <- function(x, forecast = NULL, dates = NULL,
			 show0 = FALSE, showobs = TRUE, showprior = TRUE,
			 forecastlty = NULL, tcl = -0.4,
			 mar = c(1.9, 1.9, 1.7, .5), mgp = c(2, .6, 0),
			 simobj = NULL, newdata = NULL, ...) {
  # Helper values
  params <- sampled_parameters(x)
  heavy_tailed_sv <- "nu" %in% params
  plot_volatility_series <- thinning(x)$time == "all"
  npara <- length(params)
  if (is.null(simobj) && !is.null(x$simobj)) {
    message("Simulation object extracted from input")
    simobj <- x$simobj
  }

  # Set layout
  index <- 1L
  chart_indices <- c()
  if (plot_volatility_series) {  # volatility chart(s)
    chart_indices <- c(chart_indices, rep_len(index, length.out = npara))
    index <- index + 1L
  }
  if (!x$resampled) {
    chart_indices <- c(chart_indices, index - 1L + seq_len(npara))  # parameter trace plots
    index <- index + npara
  }
  chart_indices <- c(chart_indices, index - 1L + seq_len(npara))  # parameter density plots
  index <- index + npara
  oldpar <- par(mfrow=c(1,1))  # 'layout' and 'par' do not play well together
  layout(matrix(chart_indices, ncol = npara, byrow = TRUE))

  # Plotting
  if (plot_volatility_series) {
    volplot(x, dates = dates, show0 = show0, forecast = forecast,
            forecastlty = forecastlty, tcl = tcl, mar = mar,
            mgp = mgp, simobj = simobj, newdata = newdata, ...)
  }
  if (!x$resampled) {
    paratraceplot(x, mar = mar, mgp = mgp, simobj = simobj, ...)
  } else {
    message("Traceplots are not shown after re-sampling")
  }
  paradensplot(x, showobs = showobs, showprior = showprior,
               showxlab = FALSE, mar = mar, mgp = mgp, simobj = simobj, ...)
  par(oldpar)
  invisible(x)
}

# modified density plot (from coda package)
# x is an mcmc.list with a single variable
mydensplot <- function(x, show.obs = TRUE, bwf, main = "", ylim, cutat=c(-Inf, Inf), showxlab=TRUE, mgp = c(2,.6,0), tcl=-.4, ...) {
  varname <- coda::varnames(x[[1]])
  draws <- coda::niter(x[[1]])
  for (i in seq_along(x)) {
    x_i <- as.vector(x[[i]])
    range_x_i <- range(x_i)
    if (range_x_i[1] < cutat[1] || range_x_i[2] > cutat[2]) {
      stop("Argument 'cutat' does not include range of variable.")
    }
    if (missing(bwf)) {
      bwf <- function(xx) {
        xx <- xx[!is.na(as.vector(xx))]
        return(1.06 * min(sd(xx), IQR(xx)/1.34) * length(xx)^-0.2)
      }
    }
    bw <- bwf(x_i)
    width <- 4 * bw
    if (max(abs(x_i - floor(x_i))) == 0 || bw == 0) {
      hist(x_i, prob = TRUE, main = main, col = i, ...)
    } else {
      density_scale <- "open"
      cut_at_bottom <- isTRUE(is.finite(cutat[1]) && range_x_i[1]-cutat[1] < 2*bw)
      cut_at_top <- isTRUE(is.finite(cutat[2]) && cutat[2]-range_x_i[2] < 2*bw)
      x_i <- if (cut_at_bottom && cut_at_top) {
        c(x_i, 2*cutat[1] - x_i, 2*cutat[2] - x_i)
      } else if (cut_at_bottom) {
        c(x_i, 2*cutat[1] - x_i)
      } else if (cut_at_top) {
        c(x_i, 2*cutat[2] - x_i)
      } else {
        x_i
      }
      upscale <- 1 + cut_at_bottom + cut_at_top
      dens <- density(x_i, width = width)
      index <- dens$x >= cutat[1] & dens$x <= cutat[2]
      dens$x <- dens$x[index]
      dens$y <- upscale * dens$y[index]

      if (missing(ylim)) {
        ylim <- c(0, max(dens$y))
      }
      if (i == 1) {
        plot(dens, ylab = "", main = main, type = "l",
             ylim = ylim, xlab="", mgp = mgp, tcl = tcl, col = i, ...)
      } else {
        lines(dens$x, dens$y, type = "l", col = i)
      }
      if (show.obs) {
        lines(x_i[seq_len(draws)], rep(max(dens$y)/100, draws),
              type = "h", col = i)
      }
    }
  }
  if(isTRUE(showxlab)) {
    mtext(paste("N =", coda::niter(x[[1]]), "  Bandwidth =",
                formatC(dens$bw)), side=1, line=2.7, cex=.7)
  }
  invisible(x)
}

# modified traceplot (from coda)
# x is an mcmc.list with a single variable
mytraceplot <- function (x, type = "l", ylab = "", xlab = "Iterations", mgp = c(2,.6,0), tcl = -.4, ...) {
  varname <- coda::varnames(x[[1]])
  xp <- as.vector(time(x[[1]]))
  yp <- simplify2array(x)
  graphics::matplot(xp, yp, xlab = xlab, ylab = ylab, type = type,
                    col = seq_along(x), mgp = mgp, tcl = tcl, ...)
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
    correction <- sqrt(nu / (nu - 2))
    xlab <- paste("Theoretical quantiles from a t-distribution with", round(nu, 2), "df")
  } else {
    terr <- FALSE
    nu <- Inf
    correction <- 1
    xlab <- "Theoretical quantiles from a standard normal distribution"
  }

  plot.default(x, ylab = paste("M", substring(attr(x, "type"), 2), ' standardized residuals', sep = ""),
               xlab='Time', xaxt='n', ylim=c(-1,1)*max(abs(x)),
               main=mains[1], ...)

  if (!terr) abline(h=qnorm(c(.025, .975)), lty=2)
  where <- seq(1, length(x), length=min(7, length(x)))
  axis(1, at = where, labels = gsub("r_", "", names(x)[where]))
  qqplot(qt(ppoints(length(x)), df = nu), correction * x,
         main=paste(mains[2], "for", attr(x, "type"), "standardized residuals"),
         xlab = xlab, ylab = "Sample quantiles")
  qqline(x, probs = c(0.01, 0.99), distribution = function(x) qt(x, df = nu))
  par(oldpar)
  invisible(x)
}

