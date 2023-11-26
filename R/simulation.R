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
#' \describe{
#' \item{y}{vector of length \code{len} containing the simulated data,
#' usually interpreted as ``log-returns''.}
#' \item{vol}{vector of length
#' \code{len} containing the simulated instantaneous volatilities.
#' These are \eqn{e^{h_t/2}}{exp(h_t/2)} if \code{nu == Inf}, and they are
#' \eqn{e^{h_t/2} \sqrt{\tau_t}}{exp(h_t/2) * sqrt(tau_t)} for finite \code{nu}.}
#' \item{vol0}{The initial volatility \code{exp(h_0/2)},
#' drawn from the stationary distribution of the latent AR(1) process.}
#' \item{para}{a named list with five elements \code{mu}, \code{phi},
#' \code{sigma}, \code{nu}, and \code{rho}, containing
#' the corresponding arguments.}
#' \item{latent}{vector of the latent state space \eqn{h_t} for \eqn{t > 0}.}
#' \item{latent0}{initial element of the latent state space \eqn{h_0}.}
#' \item{tau}{vector of length \code{len} containing the simulated auxiliary
#' variables for the Student-t residuals when \code{nu} is finite. More precisely,
#' \eqn{\tau_t\sim\text{Gamma}^{-1}(\text{shape}=\nu/2, \text{rate}=\nu/2-1)}{%
#' tau_t ~ InverseGamma(shape = nu / 2, rate = nu / 2 - 1)}.}
#' }
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
#' @example inst/examples/svsim.R
#' @export
svsim <- function(len, mu = -10, phi = 0.98, sigma = 0.2, nu = Inf, rho = 0) {

  name_len <- "length of simulated data set"
  assert_numeric(len, name_len)
  assert_single(len, name_len)
  assert_ge(len, 2, name_len)

  name_mu <- "input parameter mu"
  assert_numeric(mu, name_mu)
  assert_single(mu, name_mu)

  name_phi <- "input parameter phi"
  assert_numeric(phi, name_phi)
  assert_single(phi, name_phi)
  assert_gt(phi, -1, name_phi)
  assert_lt(phi, 1, name_phi)

  name_sigma <- "input parameter sigma"
  assert_numeric(sigma, name_sigma)
  assert_single(sigma, name_sigma)
  assert_positive(sigma, name_sigma)

  name_nu <- "input parameter nu"
  assert_numeric(nu, name_nu)
  assert_single(nu, name_nu)
  assert_gt(nu, 2, name_nu)

  name_rho <- "input parameter rho"
  assert_numeric(rho, name_rho)
  assert_single(rho, name_rho)
  assert_gt(rho, -1, name_rho)
  assert_lt(rho, 1, name_rho)

  len <- as.integer(len)

  h <- rep_len(as.numeric(NA), length.out=len)
  h0 <- rnorm(1, mean=mu, sd=sigma/sqrt(1-phi^2))
  tau <- if (is.finite(nu)) {
    1/rgamma(len, shape=nu/2, rate=nu/2-1)
  } else {
    rep_len(1, length.out=len)
  }
  eta <- rnorm(len)
  eps <- rho*eta + sqrt(1-rho^2)*rnorm(len)

  # simulate w/ simple loop
  h[1] <- mu + phi*(h0-mu) + sigma*rnorm(1)  # same marginal distribution as h0
  for (i in seq_len(len-1)) {
    h[i+1] <- mu + phi*(h[i]-mu) + sigma*eta[i]
  }
  y <- exp(h/2) * sqrt(tau) * eps  # "log-returns"

  ret <- list(y = y,
              vol0 = exp(h0/2),
              vol = exp(h/2) * sqrt(tau),
              para = list(mu = mu,
                          phi = phi,
                          sigma = sigma,
                          rho = rho,
                          nu = nu),
              latent = h,
              latent0 = h0,
              tau = tau)
  class(ret) <- "svsim"
  ret
}

#' @export
print.svsim <- function(x, ...) {
  cat("\nSimulated time series consisting of ", length(x$y), " observations.\n\n",
      "Parameters: level of latent variable                  mu = ", x$para$mu, "\n",
      "            persistence of latent variable           phi = ", x$para$phi, "\n",
      "            standard deviation of latent variable  sigma = ", x$para$sigma, "\n",
      "            degrees of freedom parameter              nu =", x$para$nu, "\n",
      "            leverage effect parameter                rho =", x$para$rho, "\n", sep="")
  cat("\nSimulated initial volatility:", x$vol0*x$correction, "\n")
  cat("\nSimulated volatilities:\n")
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
  ret$y <- summary(100*object$y)
  ret$vol0 <- 100*object$vol0
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
  cat(x$vol0, "\n")
  cat("\nSummary of simulated volatilities (in %):\n")
  print(x$vol)
  cat("\nSummary of simulated data (in %):\n")
  print(x$y)
  invisible(x)
}

