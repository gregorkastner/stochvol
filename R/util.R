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

#' @describeIn logret Log returns of vectors
#' @family utilities
#' @export
logret.default <- function (dat, demean = FALSE, standardize = FALSE, ...) {
  logretx <- tail(diff(log(dat)), length(dat) - 1)
  if (all(isTRUE(demean))) logretx <- logretx - mean(logretx)  # TODO 'all' not needed
  if (all(isTRUE(standardize))) logretx <- logretx / sd(logretx)
  logretx
}

#' Specify Prior Distributions for SV Models
#' 
#' This function gives access to a larger set of prior distributions
#' in case the default choice is unsatisfactory.
#' @param mu one of sv_normal and sv_constant
#' @param phi one of sv_beta and sv_constant. If sv_beta, then the specified beta distribution is the prior on (phi+1)/2
#' @param sigma2 one of sv_gamma, sv_inverse_gamma, and sv_constant
#' @param nu one of sv_infinity, sv_exponential, and sv_constant
#' @param rho one of sv_beta and sv_constant
#' @param latent0_variance one of \code{"stationary"} and sv_normal
#' @param beta an sv_multinormal object
#' @family priors
#' @export
specify_priors <- function (mu = sv_normal(mean = 0, sd = 100),
                            phi = sv_beta(alpha = 15, beta = 0.5),
                            sigma2 = sv_gamma(shape = 0.5, rate = 0.5),
                            nu = sv_infinity(),
                            rho = sv_constant(0),
                            latent0_variance = "stationary",
                            beta = sv_multinormal(mean = 0, sd = 100, dim = 1)) {
  # Validation
  ## Check mu, phi, sigma2, nu, rho, and beta
  sv_inherits <- function (x, whatlist) {
    isTRUE(any(sapply(whatlist, function (what, xx) inherits(xx, what), xx = x)))
  }
  enabled_distributions <-
    list(list(x = mu, name = "mu", whatlist = c("sv_constant", "sv_normal", "sv_constant")),
         list(x = phi, name = "phi", whatlist = c("sv_constant", "sv_beta")),
         list(x = sigma2, name = "sigma2", whatlist = c("sv_constant", "sv_gamma", "sv_inverse_gamma")),
         list(x = nu, name = "nu", whatlist = c("sv_constant", "sv_infinity", "sv_exponential")),
         list(x = rho, name = "rho", whatlist = c("sv_constant", "sv_beta")),
         list(x = beta, name = "beta", whatlist = c("sv_multinormal")))
  lapply(enabled_distributions,
         function (x) {
           with(x, if (!sv_inherits(x, whatlist)) {
             stop(name, " should inherit from one of ", prettify(whatlist), "; not ", class(x))
           })
         })
  if (inherits(rho, "sv_constant") && (rho$value <= -1 || rho$value >= 1)) {
    stop("Fixed rho needs to be in range (-1, 1); got rho = ", rho$value)
  }
  ## Check latent0_variance
  if (sv_inherits(latent0_variance, "sv_constant")) {
    assert_positive(latent0_variance$value, "The provided variance for latent0")
  } else if (!isTRUE(latent0_variance == "stationary")) {
    stop("Currently implemented options for 'latent0_variance' are either the string \"stationary\" or an sv_constant object; received ", latent0_variance)
  }
  
  structure(list(mu = mu, phi = phi, sigma2 = sigma2,
                 nu = nu, rho = rho,
                 latent0_variance = latent0_variance, beta = beta),
            class = "sv_priorspec")
}

#' Prior Distributions in \code{stochvol}
#' 
#' The functions below can be supplied to \link{specify_priors}
#' to overwrite the default set of prior distributions in \link{svsample}.
#' The functions have \code{mean}, \code{density} and \code{print} methods.
#' @rdname sv_prior
#' @family priors
#' @export
sv_constant <- function (value) {
  assert_single(value, "sv_constant")
  assert_numeric(value, "sv_constant")

  if (is.finite(value)) {
    structure(list(value = value),
              class = c("sv_constant", "sv_distribution"))
  } else if (isTRUE(value == Inf)) {
    sv_infinity()
  }
}
#' @export
mean.sv_constant <- function (x, ...) {
  x$value
}
#' @export
print.sv_constant <- function (x, ...) {
  cat("Constant(value = ", x$value, ")\n", sep = "")
}
#' @export
density.sv_constant <- function (x, ...) {
  dist <- x
  function (x) {
    as.numeric(abs(x - dist$value) < sqrt(.Machine$double.eps))
  }
}

#' @rdname sv_prior
#' @family priors
#' @export
sv_normal <- function (mean = 0, sd = 1) {
  assert_single(mean, "mean of sv_normal")
  assert_numeric(mean, "mean of sv_normal")

  assert_single(sd, "sd of sv_normal")
  assert_positive(sd, "sd of sv_normal")

  structure(list(mean = mean, sd = sd),
            class = c("sv_normal", "sv_distribution"))
}
#' @export
mean.sv_normal <- function (x, ...) {
  x$mean
}
#' @export
print.sv_normal <- function (x, ...) {
  cat("Normal(mean = ", x$mean, ", sd = ", x$sd, ")\n", sep = "")
}
#' @export
density.sv_normal <- function (x, ...) {
  dist <- x
  function (x) {
    dnorm(x, mean = dist$mean, sd = dist$sd)
  }
}

#' @section Multivariate Normal:
#' Multivariate normal objects can be specified several ways. The most general way is by calling
#' \code{sv_multinormal(mean, precision)}, which allows for arbitrary mean and (valid) precision
#' arguments. Constant mean vectors and constant diagonal precision matrices of dimension \code{D}
#' can be created two ways: either \code{sv_multinormal(mean, sd, dim = D)} or
#' \code{rep(sv_normal(mean, sd), length.out = D)}.
#' @rdname sv_prior
#' @family priors
#' @export
sv_multinormal <- function (mean = 0, precision = NULL, sd = 1, dim = NA) {
  if (!is.null(precision)) {
    assert_numeric(mean, "mean of sv_multinormal")
    assert_numeric(precision, "precision of sv_multinormal")
    
    if (!isTRUE(is.matrix(precision))) {
      stop("precision of sv_multinormal = ", precision, " should be a matrix.")
    }
    if (!isTRUE(dim(precision)[1] == dim(precision)[2])) {
      stop("precision of sv_multinormal = ", precision, " should be a square matrix.")
    }
    if (!isTRUE(length(mean) == dim(precision)[1])) {
      stop("the mean of sv_multinormal = ", mean, " and the precision of sv_multinormal = ", precision, " should be same length/dimension.")
    }
    tryCatch(chol(precision),
             error = function (e) {
               stop("the precision of sv_multinormal = ", precision, " should be positive definite.")
             })
    
    structure(list(mean = mean, precision = precision),
              class = c("sv_multinormal", "sv_distribution"))
  } else if (!is.na(dim)) {
    rep_len(sv_normal(mean, sd), length.out = dim)
  } else {
    stop("Either mean and precision or mean, sd, and dim need to be provided. The first variant takes precedence.")
  }
}
#' @method rep sv_normal
#' @export
rep.sv_normal <- function (x, times = length.out, length.out = times, ...) {
  if (missing(times) && missing(length.out)) {
    stop("Either 'times' or 'length.out' has to be provided")
  }
  if (!identical(times, length.out)) {
    stop("Parameters 'times' and 'length.out' have to be identical when both given")
  }
  rep_len(x, length.out = length.out)
}
#' @method rep.int sv_normal
#' @export
rep.int.sv_normal <- function (x, times) {
  rep_len(x, length.out = times)
}
#' @method rep_len sv_normal
#' @export
rep_len.sv_normal <- function (x, length.out) {
  sv_multinormal(mean = rep_len(mean(x), length.out),
                 precision = diag(rep_len((x$sd)^(-2), length.out),
                                  nrow = length.out, ncol = length.out))
}
#' @export
mean.sv_multinormal <- function (x, ...) {
  x$mean
}
#' @export
print.sv_multinormal <- function (x, ...) {
  cat("MultivariateNormal(...)\n")
}
#' @export
density.sv_multinormal <- function (x, ...) {
  if (!require("mvtnorm")) {
    warning("'density.sv_multinormal' needs the 'mvtnorm' package to be installed")
    function (x) {
      NA_real_
    }
  } else {
    dist <- x
    sigma <- solve(dist$precision)
    function (x) {
      mvtnorm::dmvnorm(x, mean = dist$mean, sigma = sigma)
    }
  }
}

#' @rdname sv_prior
#' @family priors
#' @export
sv_gamma <- function (shape, rate) {
  assert_single(shape, "shape of sv_gamma")
  assert_positive(shape, "shape of sv_gamma")

  assert_single(rate, "rate of sv_gamma")
  assert_positive(rate, "rate of sv_gamma")

  structure(list(shape = shape, rate = rate),
            class = c("sv_gamma", "sv_distribution"))
}
#' @export
mean.sv_gamma <- function (x, ...) {
  x$shape / x$rate
}
#' @export
print.sv_gamma <- function (x, ...) {
  cat("Gamma(shape = ", x$shape, ", rate = ", x$rate, ")\n", sep = "")
}
#' @export
density.sv_gamma <- function (x, ...) {
  dist <- x
  function (x) {
    dgamma(x, shape = dist$shape, rate = dist$rate)
  }
}

#' @rdname sv_prior
#' @family priors
#' @export
sv_inverse_gamma <- function (shape, scale) {
  assert_single(shape, "shape of sv_inverse_gamma")
  assert_positive(shape, "shape of sv_inverse_gamma")

  assert_single(scale, "scale of sv_invgamma")
  assert_gt(scale, 2, "scale of sv_invgamma")

  structure(list(shape = shape, scale = scale),
            class = c("sv_inverse_gamma", "sv_distribution"))
}
#' @export
mean.sv_inverse_gamma <- function (x, ...) {
  x$scale / (x$shape - 1)
}
#' @export
print.sv_inverse_gamma <- function (x, ...) {
  cat("InverseGamma(shape = ", x$shape, ", scale = ", x$scale, ")\n", sep = "")
}
#' @export
density.sv_inverse_gamma <- function (x, ...) {
  dist <- x
  function (x) {
    ifelse(x == 0, 0, dgamma(1/x, shape = dist$shape, rate = dist$scale) * x^{-2})
  }
}

#' @rdname sv_prior
#' @family priors
#' @export
sv_beta <- function (alpha, beta) {
  assert_single(alpha, "shape of sv_beta")
  assert_positive(alpha, "shape of sv_beta")

  assert_single(beta, "rate of sv_beta")
  assert_positive(beta, "rate of sv_beta")

  structure(list(alpha = alpha, beta = beta),
            class = c("sv_beta", "sv_distribution"))
}
#' @export
mean.sv_beta <- function (x, ...) {
  x$alpha / (x$alpha + x$beta)
}
#' @export
print.sv_beta <- function (x, ...) {
  cat("Beta(a = ", x$alpha, ", b = ", x$beta, ")\n", sep = "")
}
#' @export
density.sv_beta <- function (x, ...) {
  dist <- x
  function (x) {
    dbeta(x, shape1 = dist$alpha, shape2 = dist$beta)
  }
}

#' @rdname sv_prior
#' @family priors
#' @export
sv_exponential <- function (rate) {
  assert_single(rate, "rate of sv_exponential")
  assert_positive(rate, "rate of sv_exponential")

  structure(list(rate = rate),
            class = c("sv_exponential", "sv_distribution"))
}
#' @export
mean.sv_exponential <- function (x, ...) {
  1 / x$rate
}
#' @export
print.sv_exponential <- function (x, ...) {
  cat("Exponential(rate = ", x$rate, ")\n", sep = "")
}
#' @export
density.sv_exponential <- function (x, ...) {
  dist <- x
  function (x) {
    dexp(x, rate = dist$rate)
  }
}

#' @rdname sv_prior
#' @family priors
#' @export
sv_infinity <- function () {
  structure(list(),
            class = c("sv_infinity", "sv_distribution"))
}
#' @export
mean.sv_infinity <- function (x, ...) {
  Inf
}
#' @export
print.sv_infinity <- function (x, ...) {
  cat("Infinity\n")
}
#' @export
density.sv_infinity <- function (x, ...) {
  dist <- x
  function (x) {
    ifelse(x == Inf, 1, 0)
  }
}

#' @export
print.sv_priorspec <- function(x, ...) {
  cat("Prior distributions:\n")
  cat("mu        ~ "); print(x$mu)
  cat("(phi+1)/2 ~ "); print(x$phi)
  cat("sigma^2   ~ "); print(x$sigma2)
  cat("nu        ~ "); print(x$nu)
  cat("(rho+1)/2 ~ "); print(x$rho)
  cat("\n")
}

# Find a good initial value for mu
## Posterior mean of the homoskedastic Bayesian linear regression model
## log((y_t - X_t*beta_hat)^2) + 1.27 = mu + epsilon_t
# TODO document and export and use
init_mu <- function (y, priorspec, X = NULL, beta_hat = NULL) {
  moments_prior_mu <- switch(priorspec$mu$distribution,
                             "normal" =
                               c(priorspec$mu$para$mean, priorspec$mu$para$sd^2),
                             stop("nyi"))
  e_prior_mu <- moments_prior_mu[1]; v_prior_mu <- moments_prior_mu[2]

  regression_part <- if (is.null(X)) {
    0
  } else {
    X %*% beta_hat
  }
  left_hand_side <- log((y - regression_part)^2) + 1.27
  len <- length(left_hand_side)
  ols <- mean(left_hand_side)
  (len * ols + e_prior_mu / v_prior_mu) / (len + 1 / v_prior_mu)
}

# Merge lists: match user input to a default, fill in missing parts in the input from the default
apply_default_list <- function (input, default, name_input, name_default) {
  if (is.null(input)) {
    default
  } else if (is.list(default)) {
    elements <- names(default)
    assert_element(names(input), elements, name_input, name_default)

    for (element in elements) {
      default[[element]] <- apply_default_list(input[[element]], default[[element]], paste0(name_input, "$", element), paste0(name_default, "$", element))
    }
    default
  } else {
    input
  }
}

asisprint <- function (x, censtring) {
  if (length(x) %% 2 == 0) {
    toshorten <- rep(as.character(censtring), length(x)/2)
    if (identical(toshorten, as.character(x)))
      return(sprintf("ASISx%d", length(x)/2))
  }
  paste0("(", paste(x, collapse=", "), ")")
}
