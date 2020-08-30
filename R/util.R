#' @describeIn logret Log returns of vectors
#' @family utilities
#' @export
logret.default <- function (dat, demean = FALSE, standardize = FALSE, ...) {
  logretx <- tail(diff(log(dat)), length(dat) - 1)
  if (all(isTRUE(demean))) logretx <- logretx - mean(logretx)  # TODO 'all' not needed
  if (all(isTRUE(standardize))) logretx <- logretx / sd(logretx)
  logretx
}

# Distributions used as prior in the \code{stochvol} package
# TODO document and export
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

sv_normal <- function (mean = 0, sd = 1) {
  assert_single(mean, "mean of sv_normal")
  assert_numeric(mean, "mean of sv_normal")

  assert_single(sd, "sd of sv_normal")
  assert_positive(sd, "sd of sv_normal")

  structure(list(mean = mean, stdev = sd),
            class = c("sv_normal", "sv_distribution"))
}

#' @export
mean.sv_normal <- function (x, ...) {
  x$mean
}

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

sv_infinity <- function () {
  structure(list(),
            class = c("sv_infinity", "sv_distribution"))
}

#' @export
mean.sv_infinity <- function (x, ...) {
  Inf
}

# TODO document and export
specify_priors <- function (mu = sv_normal(mean = 0, sd = 100),
                            phi = sv_beta(alpha = 15, beta = 0.5),
                            sigma2 = sv_gamma(shape = 0.5, rate = 0.5),
                            nu = sv_infinity(),
                            rho = sv_constant(0),
                            latent0 = "stationary",
                            beta = sv_normal(mean = 0, sd = 100)) {
  # Validation
  ## Check mu, phi, sigma2, nu, rho, and beta
  sv_inherits <- function (x, whatlist) {
    isTRUE(any(sapply(whatlist, function (what, xx) inherits(xx, what), xx = x)))
  }
  enabled_distributions <-
    list(list(x = mu, name = "mu", whatlist = c("sv_normal", "sv_constant")),
         list(x = phi, name = "phi", whatlist = c("sv_beta")),
         list(x = sigma2, name = "sigma2", whatlist = c("sv_gamma", "sv_inverse_gamma")),
         list(x = nu, name = "nu", whatlist = c("sv_infinity", "sv_exponential")),
         list(x = rho, name = "rho", whatlist = c("sv_constant", "sv_beta")),
         list(x = beta, name = "beta", whatlist = c("sv_normal", "sv_constant")))
  lapply(enabled_distributions,
         function (x) {
           with(x, if (!sv_inherits(x, whatlist)) {
             stop(name, " should inherit from one of ", prettify(whatlist), "; not ", class(x))
           })
         })
  ## If constant, check constant value
  if (inherits(mu, "sv_constant") && mu$value != 0) {
    stop("Fixed mu != 0 not yet implemented; got mu = ", mu$value)
  }
  if (inherits(rho, "sv_constant") && rho$value != 0) {
    stop("Fixed rho != 0 not yet implemented; got rho = ", rho$value)
  }
  if (inherits(beta, "sv_constant") && beta$value != 0) {
    stop("Fixed beta != 0 not yet implemented; got beta = ", beta$value)
  }
  ## Check latent0
  if(!(isTRUE(latent0 == "stationary") ||
       (assert_single(latent0, "latent0") &&
        assert_positive(latent0, "latent0")))) {  # TODO sv_inherits(latent0, "sv_normal"))) {
    stop("Currently implemented options for 'latent0' are either the string \"stationary\" or a single positive number; received ", latent0)
  }

  structure(list(mu = mu, phi = phi, sigma2 = sigma2,
                 nu = nu, rho = rho,
                 latent0 = latent0, beta = beta),
            class = "sv_priorspec")
}

# Find a good initial value for mu
## Posterior mean of the homoskedastic Bayesian linear regression model
## log((y_t - X_t*beta_hat)^2) + 1.27 = mu + epsilon_t
# TODO document and export
init_mu <- function (y, priorspec, X = NULL, beta_hat = NULL) {
  moments_prior_mu <- switch(priorspec$mu$distribution,
                             "normal" =
                               c(priorspec$mu$para$mean, priorspec$mu$para$stdev^2),
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
