#' Simulating a Stochastic Volatility Process
#' 
#' \code{svsim} is used to produce realizations of a stochastic volatility (SV)
#' process.
#' 
#' This function draws an initial log-volatility \code{h_0} from the stationary
#' distribution of the AR(1) process defined by \code{phi}, \code{sigma}, and \code{mu}.
#' Then the function jointly simulates the log-volatility series
#' \code{h_1,...,h_n} with the given AR(1) structure, and the ``log-return'' series
#' \code{y_1,...,y_n} with mean 0 and standard deviation \code{exp(h/2)}.
#' Additionally, for each index \code{i}, \code{y_i} can be set to have a conditionally heavy-tailed
#' residual (through \code{nu}) and/or to be correlated with \code{(h_{i+1}-h_i)}
#' (through \code{rho}, the so-called leverage effect, resulting in asymmetric ``log-returns'').
#' 
#' @param len length of the simulated time series.
#' @param mu level of the latent log-volatility AR(1) process. The defaults
#' value is \code{-10}.
#' @param phi persistence of the latent log-volatility AR(1) process. The
#' default value is \code{0.98}.
#' @param sigma volatility of the latent log-volatility AR(1) process. The
#' default value is \code{0.2}.
#' @param nu degrees-of-freedom for the conditional innovations distribution.
#' The default value is \code{Inf}, corresponding to standard normal
#' conditional innovations.
#' @param rho correlation between the observation and the increment of the
#' log-volatility. The default value is \code{0}, corresponding to the basic
#' SV model with symmetric ``log-returns''.
#' @return The output is a list object of class \code{svsim} containing
#' \item{y}{a vector of length \code{len} containing the simulated data,
#' usually interpreted as ``log-returns''.}
#' \item{vol}{a vector of length
#' \code{len} containing the simulated instantaneous volatilities
#' \code{exp(h_t/2)}.}
#' \item{vol0}{The initial volatility \code{exp(h_0/2)},
#' drawn from the stationary distribution of the latent AR(1) process.}
#' \item{para}{a named list with five elements \code{mu}, \code{phi},
#' \code{sigma}, \code{nu}, and \code{rho}, containing
#' the corresponding arguments.}
#' @note The function generates the ``log-returns'' by
#' \code{y <- exp(-h/2)*rt(h, df=nu)}. That means that in the case of \code{nu < Inf}
#' the (conditional) volatility is \code{sqrt(nu/(nu-2))*exp(h/2)}, and that corrected value
#' is shown in the \code{print}, \code{summary} and \code{plot} methods.
#' 
#' To display the output use \code{print}, \code{summary} and \code{plot}. The
#' \code{print} method simply prints the content of the object in a moderately
#' formatted manner. The \code{summary} method provides some summary statistics
#' (in \%), and the \code{plot} method plots the the simulated 'log-returns'
#' \code{y} along with the corresponding volatilities \code{vol}.
#' @author Gregor Kastner \email{gregor.kastner@@wu.ac.at}
#' @seealso \code{\link{svsample}}
#' @keywords datagen ts
#' @examples
#' 
#' ## Simulate a highly persistent SV process of length 500
#' sim <- svsim(500, phi = 0.99, sigma = 0.1)
#' 
#' print(sim)
#' summary(sim)
#' plot(sim)
#' 
#' ## Simulate an SV process with leverage
#' sim <- svsim(200, phi = 0.94, sigma = 0.15, rho = -0.6)
#' 
#' print(sim)
#' summary(sim)
#' plot(sim)
#' 
#' ## Simulate an SV process with conditionally heavy-tails
#' sim <- svsim(250, phi = 0.91, sigma = 0.05, nu = 5)
#' 
#' print(sim)
#' summary(sim)
#' plot(sim)
#' @export
svsim <- function(len, mu = -10, phi = 0.98, sigma = 0.2, nu = Inf, rho = 0) {

  # Some error checking
  if (any(is.na(len)) || !is.numeric(len) || length(len) != 1 || any(len < 1)) {
    stop("Argument 'len' (length of simulated series) must be a single number >= 2.")
  } else {
    len <- as.integer(len)
  }

  if (!is.numeric(mu) || length(mu) != 1) {
    stop("Argument 'mu' (level of latent variable) must be a single number.")
  }

  if (!is.numeric(phi) || length(phi) != 1) {
    stop("Argument 'phi' (persistence of latent variable) must be a single number.")
  }

  if (!is.numeric(sigma) || length(sigma) != 1 || sigma <= 0) {
    stop("Argument 'sigma' (volatility of latent variable) must be a single number > 0.")
  }

  if (!is.numeric(nu) || length(nu) != 1 || nu <= 2) {
    stop("Argument 'nu' (degrees of freedom for the conditional error) must be a single number > 2.")
  }

  if (!is.numeric(rho) || length(rho) != 1 || abs(rho) >= 1) {
    stop("Argument 'rho' (correlation between the observations and the volatility increments) must be a single number between -1 and 1 exclusive.")
  }

  h <- rep(as.numeric(NA), len)
  h0 <- rnorm(1, mean=mu, sd=sigma/sqrt(1-phi^2))
  standardizer <- if (is.finite(nu)) sqrt((nu-2)/nu) else 1
  eps <- rt(len, df = nu)
  eta <- rho * eps * standardizer + sqrt(1-rho^2) * rnorm(len)

  # simulate w/ simple loop
  h[1] <- mu + phi*(h0-mu) + sigma*rnorm(1)  # same marginal distribution as h0
  for (i in seq_len(len-1)) {
    h[i+1] <- mu + phi*(h[i]-mu) + sigma*eta[i]
  }
  y <- exp(h / 2) * eps  # "log-returns"

  ret <- list(y = y,
              vol = exp(h/2),
              para = list(mu = mu,
                          phi = phi,
                          sigma = sigma))
  ret$vol0 <- exp(h0/2)
  ret$para$rho <- rho
  ret$para$nu <- nu
  ret$correction <- 1/standardizer  # TODO discuss with Gregor
  class(ret) <- "svsim"
  ret
}

#' @export
print.svsim <- function(x, ...) {
  cat("\nSimulated time series consisting of ", length(x$y), " observations.\n\n",
      "Parameters: level of latent variable                  mu = ", x$para$mu, "\n",
      "            persistence of latent variable           phi = ", x$para$phi, "\n",
      "            standard deviation of latent variable  sigma = ", x$para$sigma, "\n", sep="")
  if ("nu" %in% names(x$para)) cat("            degrees of freedom parameter              nu =", x$para$nu, "\n")
  if ("rho" %in% names(x$para)) cat("            leverage effect parameter                rho =", x$para$rho, "\n")
  cat("\nSimulated initial volatility:", x$vol0*x$correction, "\n")
  cat("\nSimulated volatilities:\n")
  print(x$vol*x$correction, ...)
  cat("\nSimulated data (usually interpreted as 'log-returns'):\n")
  print(x$y, ...)
}

#' @export
plot.svsim <- function(x, mar = c(3, 2, 2, 1), mgp = c(1.8, .6, 0), ...) {
  op <- par(mfrow = c(2, 1), mar = mar, mgp = mgp)
  plot.ts(100*x$y, ylab = "", ...)
  mtext("Simulated data: 'log-returns' (in %)", cex = 1.2, line = .4, font = 2)
  plot.ts(100*x$vol*x$correction, ylab = "", ...)
  mtext("Simulated volatilities (in %)", cex = 1.2, line = .4, font = 2)
  par(op)
}

#' @export
summary.svsim <- function(object, ...) {
  ret <- vector("list")
  class(ret) <- "summary.svsim"
  ret$len <- length(object$y)
  ret$para <- object$para
  ret$vol <- summary(100*object$vol)
  ret$vol.corrected <- summary(100*object$vol*object$correction)
  ret$y <- summary(100*object$y)
  ret$vol0 <- 100*object$vol0
  ret$vol0.corrected <- 100*object$vol0*object$correction
  ret
}

#' @method print summary.svsim
#' @export
print.summary.svsim  <- function(x, ...) {
  cat("\nSimulated time series consisting of ", x$len, " observations.\n",
      "\nParameters: level of latent variable                  mu = ", x$para$mu, 
      "\n            persistence of latent variable           phi = ", x$para$phi,
      "\n            standard deviation of latent variable  sigma = ", x$para$sigma,
      "\n            degrees of freedom parameter              nu = ", x$para$nu,
      "\n            leverage effect parameter                rho = ", x$para$rho, "\n", sep="")
  cat("\nSimulated initial volatility (in %): ")
  cat(x$vol0.corrected, "\n")
  cat("\nSummary of simulated volatilities (in %):\n")
  print(x$vol.corrected)
  cat("\nSummary of simulated data (in %):\n")
  print(x$y)
  invisible(x)
}

