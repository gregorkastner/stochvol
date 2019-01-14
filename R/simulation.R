#' Simulating a Stochastic Volatility Process
#' 
#' \code{svsim} is used to produce realizations of a stochastic volatility (SV)
#' process.
#' 
#' This function draws an initial log-volatility \code{h_0} from the stationary
#' distribution of the AR(1) process and iteratively generates
#' \code{h_1,...,h_n}. Finally, the ``log-returns'' are simulated from a normal
#' distribution with mean 0 and standard deviation \code{exp(h/2)}.
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
#' @return The output is a list object of class \code{svsim} containing
#' \item{y}{a vector of length \code{len} containing the simulated data,
#' usually interpreted as ``log-returns''.} \item{vol}{a vector of length
#' \code{len} containing the simulated instantaneous volatilities
#' \code{exp(h_t/2)}.} \item{vol0}{the initial volatility \code{exp(h_0/2)},
#' drawn from the stationary distribution of the latent AR(1) process.}
#' \item{para}{a named list with three elements \code{mu}, \code{phi},
#' \code{sigma} (and potentially \code{nu}), containing the corresponding
#' arguments.}
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

  if (is.finite(nu) && !isTRUE(all.equal(rho, 0))) {  # to make life simpler
    stop("The case of both t-distributed errors and asymmetry is not yet implemented. Either argument 'rho' has to be 0 or argument 'nu' has to be infinite.")
  }

  if (isTRUE(all.equal(rho, 0))) rho <- 0

  h <- rep(as.numeric(NA), len)
  h0 <- rnorm(1, mean=mu, sd=sigma/sqrt(1-phi^2))
  eps <- rt(len, df = nu)
  eta <- rho * eps + sqrt(1-rho^2) * rnorm(len)  # either rho == 0 or eps ~ N(0, 1)

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
  if (is.finite(nu)) ret$para$nu <- nu
  if (!isTRUE(all.equal(rho, 0))) {
    ret$para$rho <- rho
    class(ret) <- c("svlsim", "svsim")
  } else {
    ret$vol0 <- exp(h0/2)
    class(ret) <- "svsim"
  }
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
  cat("\nSimulated initial conditional volatility:", x$vol0, "\n")
  cat("\nSimulated conditional volatilities:\n")
  print(x$vol, ...)
  cat("\nSimulated data (usually interpreted as 'log-returns'):\n")
  print(x$y, ...)
}

#' @export
plot.svsim <- function(x, mar = c(3, 2, 2, 1), mgp = c(1.8, .6, 0), ...) {
  op <- par(mfrow = c(2, 1), mar = mar, mgp = mgp)
  plot.ts(100*x$y, ylab = "", ...)
  mtext("Simulated data: 'log-returns' (in %)", cex = 1.2, line = .4, font = 2)
  plot.ts(100*x$vol, ylab = "", ...)
  mtext("Simulated conditional volatilities (in %)", cex = 1.2, line = .4, font = 2)
  par(op)
}

#' @export
summary.svsim <- function(object, ...) {
  ret <- vector("list")
  class(ret) <- "summary.svsim"
  ret$len <- length(object$y)
  ret$para <- object$para
  ret$vol <- summary(100*object$vol)
  ret$y <- summary(100*object$y)
  if (exists("vol0", ret)) ret$vol0 <- 100*object$vol0
  ret
}

#' @export
summary.svlsim <- function (object, ...) {
  ret <- summary.svsim(object = object, ...)
  class(ret) <- c("summary.svlsim", "summary.svsim")
}

#' @method print summary.svsim
#' @export
print.summary.svsim  <- function(x, ...) {
  cat("\nSimulated time series consisting of ", x$len, " observations.\n",
      "\nParameters: level of latent variable                  mu = ", x$para$mu, 
      "\n            persistence of latent variable           phi = ", x$para$phi,
      "\n            standard deviation of latent variable  sigma = ", x$para$sigma, "\n", sep="")
  if ("nu" %in% names(x$para)) cat("            degrees of freedom parameter              nu =", x$para$nu, "\n")
  if ("rho" %in% names(x$para)) cat("            leverage effect parameter                 rho =", x$para$rho, "\n")

  if (exists("vol0", x)) {
    cat("\nSimulated initial conditional volatility (in %): ")
    cat(x$vol0, "\n")
  }
  cat("\nSummary of simulated conditional volatilities (in %):\n")
  print(x$vol)
  cat("\nSummary of simulated data (in %):\n")
  print(x$y)
  invisible(x)
}

