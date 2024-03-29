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

#' @describeIn logret Log returns of vectors
#' @family utilities
#' @method logret default
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
#' @param mu one of \code{sv_normal} or \code{sv_constant}
#' @param phi one of \code{sv_beta}, \code{sv_normal}, or \code{sv_constant}. If \code{sv_beta}, then the specified beta distribution is the prior for \code{(phi+1)/2}
#' @param sigma2 one of \code{sv_gamma}, \code{sv_inverse_gamma}, or \code{sv_constant}
#' @param nu one of \code{sv_infinity}, \code{sv_exponential}, or \code{sv_constant}. If \code{sv_exponential}, then the specified exponential distribution is the prior for \code{nu-2}
#' @param rho one of \code{sv_beta} or \code{sv_constant}. If \code{sv_beta}, then the specified beta distribution is the prior for \code{(rho+1)/2}
#' @param latent0_variance either the character string \code{"stationary"} or an \code{sv_constant} object.
#' If \code{"stationary"}, then h0 ~ N(\code{mu}, \code{sigma^2/(1-phi^2)}). If an \code{sv_constant} object with value \code{v}, then h0 ~ N(\code{mu}, \code{sigma^2/v}).
#' Here, N(b, B) stands for mean b and variance B
#' @param beta an \code{sv_multinormal} object
#' @family priors
#' @export
specify_priors <- function (mu = sv_normal(mean = 0, sd = 100),
                            phi = sv_beta(shape1 = 5, shape2 = 1.5),
                            sigma2 = sv_gamma(shape = 0.5, rate = 0.5),
                            nu = sv_infinity(),
                            rho = sv_constant(0),
                            latent0_variance = "stationary",
                            beta = sv_multinormal(mean = 0, sd = 10000, dim = 1)) {
  # Validation
  ## Check mu, phi, sigma2, nu, rho, and beta
  sv_inherits <- function (x, whatlist) {
    isTRUE(any(sapply(whatlist, function (what, xx) inherits(xx, what), xx = x)))
  }
  enabled_distributions <-
    list(list(x = mu, name = "mu", whatlist = c("sv_constant", "sv_normal", "sv_constant")),
         list(x = phi, name = "phi", whatlist = c("sv_constant", "sv_beta", "sv_normal")),
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
  ## Check constant values
  if (sv_inherits(phi, "sv_constant")) {
    assert_gt(phi$value, -1, "The provided constant value for phi")
    assert_lt(phi$value, 1, "The provided constant value for phi")
  }
  if (sv_inherits(sigma2, "sv_constant")) {
    assert_positive(sigma2$value, "The provided constant value for sigma2")
  }
  if (sv_inherits(nu, "sv_constant")) {
    assert_gt(nu$value, 2, "The provided constant value for nu")
  }
  if (sv_inherits(rho, "sv_constant")) {
    assert_gt(rho$value, -1, "The provided constant value for rho")
    assert_lt(rho$value, 1, "The provided constant value for rho")
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
#' The functions below can be supplied to \code{\link{specify_priors}}
#' to overwrite the default set of prior distributions in \code{\link{svsample}}.
#' The functions have \code{mean}, \code{range}, \code{density}, and
#' \code{print} methods.
#' @param value The constant value for the degenerate constant distribution
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
#' @export
range.sv_constant <- function (x, na.rm = FALSE, ...) {
  rep_len(x$value, length.out = 2)
}

#' @param mean Expected value for the univariate normal distribution or mean vector of the multivariate normal distribution
#' @param sd Standard deviation for the univariate normal distribution or constant scale of the multivariate normal distribution
#' @rdname sv_prior
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
#' @export
range.sv_normal <- function (x, na.rm = FALSE, ...) {
  c(-Inf, Inf)
}

#' @section Multivariate Normal:
#' Multivariate normal objects can be specified several ways. The most general way is by calling
#' \code{sv_multinormal(mean, precision)}, which allows for arbitrary mean and (valid) precision
#' arguments. Constant mean vectors and constant diagonal precision matrices of dimension \code{D}
#' can be created two ways: either \code{sv_multinormal(mean, sd, dim = D)} or
#' \code{rep(sv_normal(mean, sd), length.out = D)}.
#' @param precision Precision matrix for the multivariate normal distribution
#' @param dim (optional) Dimension of the multivariate distribution
#' @rdname sv_prior
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
#' @rawNamespace S3method(rep.int,sv_normal,rep_int_sv_normal)
rep_int_sv_normal <- function (x, times) {
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
"[.sv_multinormal" <- function (x, n, drop=TRUE) {
  if (drop) {
    sv_normal(mean = x$mean[n], sd = 1/sqrt(x$precision[n, n]))
  } else {
    sv_multinormal(mean = x$mean[n], precision = x$precision[n, n])
  }
}
#' @export
mean.sv_multinormal <- function (x, ...) {
  x$mean
}
#' @export
print.sv_multinormal <- function (x, short = FALSE, ...) {
  if (length(x$mean) == 1) {
    print(x[1])
  } else {
    if (short) {
      cat("MultivariateNormal(...)\n")
    } else {
      cat("MultivariateNormal with mean vector\n    (", paste(x$mean, collapse = ", "), ")\n    and precision matrix\n")
      print(x$precision)
    }
  }
}
#' @export
density.sv_multinormal <- function (x, ...) {
  if (!requireNamespace("mvtnorm")) {
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
#' @export
range.sv_multinormal <- function (x, na.rm = FALSE, ...) {
  stop("Function 'range' undefined for class 'sv_multinormal'")
}

#' @param shape Shape parameter for the distribution
#' @param rate Rate parameter for the distribution
#' @rdname sv_prior
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
#' @export
range.sv_gamma <- function (x, na.rm = FALSE, ...) {
  c(0, Inf)
}

#' @param scale Scale parameter for the distribution
#' @rdname sv_prior
#' @export
sv_inverse_gamma <- function (shape, scale) {
  assert_single(shape, "shape of sv_inverse_gamma")
  assert_gt(shape, 2, "shape of sv_inverse_gamma")

  assert_single(scale, "scale of sv_invgamma")
  assert_positive(scale, "scale of sv_invgamma")

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
#' @export
range.sv_inverse_gamma <- function (x, na.rm = FALSE, ...) {
  c(0, Inf)
}

#' @param shape1 First shape parameter for the distribution
#' @param shape2 Second shape parameter for the distribution
#' @rdname sv_prior
#' @export
sv_beta <- function (shape1, shape2) {  # rename 2beta_m1
  assert_single(shape1, "shape of sv_beta")
  assert_positive(shape1, "shape of sv_beta")

  assert_single(shape2, "rate of sv_beta")
  assert_positive(shape2, "rate of sv_beta")

  structure(list(shape1 = shape1, shape2 = shape2),
            class = c("sv_beta", "sv_distribution"))
}
#' @export
mean.sv_beta <- function (x, ...) {
  x$shape1 / (x$shape1 + x$shape2)
}
#' @export
print.sv_beta <- function (x, ...) {
  cat("Beta(a = ", x$shape1, ", b = ", x$shape2, ")\n", sep = "")
}
#' @export
density.sv_beta <- function (x, ...) {
  dist <- x
  function (x) {
    dbeta(x, shape1 = dist$shape1, shape2 = dist$shape2)
  }
}
#' @export
range.sv_beta <- function (x, na.rm = FALSE, ...) {
  c(0, 1)
}

#' @param rate Rate parameter for the distribution
#' @rdname sv_prior
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
#' @export
range.sv_exponential <- function (x, na.rm = FALSE, ...) {
  c(0, Inf)
}

#' @rdname sv_prior
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
range.sv_infinity <- function (x, na.rm = FALSE, ...) {
  c(Inf, Inf)
}

#' @export
print.sv_priorspec <- function(x, showbeta = FALSE, ...) {
  cat("Prior distributions:\n")
  cat("mu        ~ "); print(x$mu)
  if (inherits(x$phi, "sv_beta")) {
    cat("(phi+1)/2 ~ ")
  } else {
    cat("phi       ~ ")
  }
  print(x$phi)
  cat("sigma^2   ~ "); print(x$sigma2)
  if (inherits(x$nu, "sv_exponential")) {
    cat("nu-2      ~ "); print(x$nu)
  } else {
    cat("nu        ~ "); print(x$nu)
  }
  if (inherits(x$nu, "sv_beta")) {
    cat("(rho+1)/2 ~ "); print(x$rho)
  } else {
    cat("rho       ~ "); print(x$rho)
  }
  if (showbeta) {
    cat("beta      ~ "); print(x$beta, short = TRUE)
  }
}

# Find good initialization for beta
## Simple OLS
init_beta <- function (y, X) {
  stats::coefficients(stats::lm(y ~ 0 + X))
}

# Find a good initial value for mu
## Posterior mean of the homoskedastic Bayesian linear regression model
## log((y_t - X_t*beta_hat)^2) = mu + epsilon_t
## where epsilon_t ~ N(-1.27, 4.934)
init_mu <- function (y, priorspec, X = NULL, beta_hat = NULL) {
  laplace_approx_mean <- -1.27; laplace_approx_var <- 4.934
  regression_part <- if (is.null(X)) {
    0
  } else {
    X %*% beta_hat
  }
  left_hand_side <- (log((y - regression_part)^2 + 1e-20) - laplace_approx_mean)
  len <- length(left_hand_side)
  ols <- mean(left_hand_side)

  if (inherits(priorspec$mu, "sv_normal")) {
    moments_prior_mu <- c(mean(priorspec$mu), priorspec$mu$sd^2)
    e_prior_mu <- moments_prior_mu[1]; v_prior_mu <- moments_prior_mu[2]
    (v_prior_mu * len * ols + laplace_approx_var * e_prior_mu) / (v_prior_mu * len + laplace_approx_var)
  } else {
    ols
  }
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

