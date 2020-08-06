#' @describeIn logret Log returns of vectors
#' @family utilities
#' @export
logret.default <- function (dat, demean = FALSE, standardize = FALSE, ...) {
  logretx <- tail(diff(log(dat)), length(dat) - 1)
  if (all(isTRUE(demean))) logretx <- logretx - mean(logretx)  # TODO 'all' not needed
  if (all(isTRUE(standardize))) logretx <- logretx / sd(logretx)
  logretx
}

# Find a good initial value for mu
## Posterior mean of the homoskedastic Bayesian linear regression model
## log((y_t - X_t*beta_hat)^2) + 1.27 = mu + epsilon_t
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

asisprint <- function (x, censtring) {
  if (length(x) %% 2 == 0) {
    toshorten <- rep(as.character(censtring), length(x)/2)
    if (identical(toshorten, as.character(x)))
      return(sprintf("ASISx%d", length(x)/2))
  }
  paste0("(", paste(x, collapse=", "), ")")
}
