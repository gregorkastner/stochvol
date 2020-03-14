# Functions calculating criteria

# Log-likelihood functions
#
# @param y vector of observations
# @param h vector of log-variance values
# @param theta vector of parameters, order: mu, phi, sigma, and rho or nu
loglik_svl <- function (y, h, theta) {
  mu <- theta[1]
  phi <- theta[2]
  sigma <- theta[3]
  rho <- theta[4]

  rho_const <- sqrt(1-rho^2)
  n <- length(y)
  result <- 0
  result <- result + sum(dnorm(head(y, -1), exp(.5*head(h, -1))*rho*(tail(h, -1)-mu-phi*(head(h, -1)-mu))/sigma, exp(.5*head(h, -1))*rho_const, TRUE))
  result <- result + dnorm(y[n], 0, exp(.5*h[n]), TRUE)
  result
}

loglik_svt <- function (y, h, theta) {
  mu <- theta[1]
  phi <- theta[2]
  sigma <- theta[3]
  nu <- theta[4]

  n <- length(y)
  sum(-.5*h + dt(y*exp(-.5*h), nu, TRUE))
}

loglik_sv <- function (y, h, theta) {
  mu <- theta[1]
  phi <- theta[2]
  sigma <- theta[3]

  n <- length(y)
  sum(dnorm(y, 0, exp(.5*h), TRUE))
}

#' Deviance information criterion
#' @export
dic <- function (x, ...) {
  UseMethod("dic", x)
}

#' Deviance information criterion for stochastic volatility models
#' @method dic svdraws
dic.svdraws <- function (x) {
  if (!inherits(x, "svdraws")) stop("This function expects an 'svdraws' object.")
  fnlik <- if (NCOL(para(x)) == 3) {
    loglik_sv
  } else if ("nu" %in% names(x$priors)) {
    loglik_svt
  } else if ("rho" %in% names(x$priors)) {
    loglik_svl
  } else {
    stop("Object corrupted. Please contact maintainer.")
  }
  
  fndev <- function (...) -2 * fnlik(...)  # deviance function
  # not the most efficient implementation but it probably doesn't matter in the use cases
  y <- x$y
  n <- length(y)
  if (n == 0) stop("Object corrupted, zero length. Please contact maintainer.")
  deviances <- rep.int(NA, n)
  for (i in seq_along(deviances)) {
    deviances[i] <- fndev(y, latent(x)[i, ], para(x)[i, ])
  }
  D_theta_hat <- fndev(y, colMeans(latent(x)), colMeans(para(x)))
  D_avg <- mean(deviances)
  p_D_1 <- D_avg - D_theta_hat
  p_D_2 <- .5 * var(deviances)

  c(DIC1 = D_avg + p_D_1,
    DIC2 = D_avg + p_D_2,
    D.avg = D_avg, p.D.1 = p_D_1, p.D.2 = p_D_2)
}


