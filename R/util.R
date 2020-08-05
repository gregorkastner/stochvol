#' @describeIn logret Log returns of vectors
#' @family utilities
#' @export
logret.default <- function (dat, demean = FALSE, standardize = FALSE, ...) {
  logretx <- tail(diff(log(dat)), length(dat) - 1)
  if (all(isTRUE(demean))) logretx <- logretx - mean(logretx)  # TODO 'all' not needed
  if (all(isTRUE(standardize))) logretx <- logretx / sd(logretx)
  logretx
}

#priorspec <-
#  list(mu = list(distribution = "normal",
#                 para = list(mean = 0, stdev = 100)),
#       phi = list(distribution = "beta",
#                  para = list(alpha = 15, beta = 0.5)),
#       sigma2 = list(distribution = "gamma",
#                     para = list(shape = 0.5, rate = 1)))
#y <- log(stochvol::svsim(20000, mu = -10, phi = 0.95, sigma = 0.8)$y^2)
#init_sv_aux_regression(y, priorspec)
# Heuristic stepwise auxiliary conjugate Bayesian regressions for finding good starting values
## 1. Non-SV regressors (betas)
### We assume time-independent latents
### y_t = X_t * beta + Sigma * epsilon_t
### where Sigma ~ InvGamma(a, b). The hyperparameters 'a', and 'b' are matched to the moments of
### exp(h/2) with h ~ N(mu, sigma^2/(1-phi^2)), and we use the priors of mu, phi, and sigma.
### We take \hat{beta} = E[beta | y, X, hyperparameters(mu, sigma, phi)], a.k.a. the posterior mean.
init_beta_regression <- function (y, X, priorspec) {
  stop("nyi")
}
## 2. Vanilla SV parameters (mu, phi, sigma)
### We take Y_t <- log((y_t - X_t * \hat{beta})^2) + 1.27 for this part
### Bayesian linear regression
### Y_t = mu * (1-phi) + phi * Y_{t-1} + \sqrt{sigma^2 + 4.934*(1+phi^2)} * ksi_t
### Again prior-induced moment matched conjugate priors are given to
### alpha = mu*(1-phi),
### beta = phi,
### Sigma = \sqrt{sigma^2 + 4.934*(1+phi^2)}
### We again take the posterior mean of that regression but w.r.t. the original parameters
### meaning that a density transformation and only then an expectation is needed.
init_sv_aux_regression <- function (y, priorspec) {
  y <- y + 1.27  # -1.27 == mean(normal_approximation_to_chi2)

  # SV parameters prior moments
  ## mu
  moments_prior_mu <- switch(priorspec$mu$distribution,
                             "normal" =
                               c(priorspec$mu$para$mean, priorspec$mu$para$stdev^2),
                             stop("nyi"))
  e_prior_mu <- moments_prior_mu[1]; v_prior_mu <- moments_prior_mu[2]
  ## phi
  moments_prior_phi <- switch(priorspec$phi$distribution,
                              "beta" =
                                (function (a, b)
                                 c((a - b) / (a + b), 4 * a * b / (a + b)^2 / (1 + a + b),
                                   16*a*b * (a^3 - a^2*(b-2) - a*(b^2 + 2*b - 1) + b*(1+b)^2) /
                                     (a+b)^2 / (1+a+b)^2 / (2+a+b) / (3+a+b)))
                                (priorspec$phi$para$alpha, priorspec$phi$para$beta),
                              stop("nyi"))
  e_prior_phi <- moments_prior_phi[1]; v_prior_phi <- moments_prior_phi[2]; v_prior_phi2 <- moments_prior_phi[3]
  ## sigma2
  moments_prior_sigma2 <- switch(priorspec$sigma2$distribution,
                                 "gamma" =
                                   (function (a, b)
                                    c(a / b, a / b^2))
                                   (priorspec$sigma2$para$shape, priorspec$sigma2$para$rate),
                                 "inverse_gamma" =
                                   (function (a, b)
                                    c(b / (a - 1), (b / (a - 1))^2 / (a - 2)))
                                   (priorspec$sigma2$para$shape, priorspec$sigma2$para$scale),
                                 stop("nyi"))
  e_prior_sigma2 <- moments_prior_sigma2[1]; v_prior_sigma2 <- moments_prior_sigma2[2]

  # Auxiliary regression (unconditional) prior moments
  v_normal_approx <- 4.934
  ## alpha & beta
  e_prior_alpha <- e_prior_mu * (1 - e_prior_phi)
  e_prior_beta <- e_prior_phi
  v_prior_alpha <- v_prior_mu * (v_prior_phi + (1 - e_prior_phi)^2) + e_prior_mu^2 * v_prior_phi
  v_prior_beta <- v_prior_phi
  cov_prior_alpha_beta <- -e_prior_mu * v_prior_phi
  ## Sigma^2
  e_prior_Sigma2 <- e_prior_sigma2 + v_normal_approx * (1 + v_prior_phi + e_prior_phi^2)
  v_prior_Sigma2 <- v_prior_sigma2 + v_normal_approx^2 * v_prior_phi2
  a_prior_Sigma2 <- e_prior_Sigma2^2 / v_prior_Sigma2 + 2
  b_prior_Sigma2 <- e_prior_Sigma2 * (e_prior_Sigma2^2 / v_prior_Sigma2 + 1)
  ## Notation of Bayesian linear regressions
  mu_0 <- rbind(e_prior_alpha, e_prior_beta)
  lambda_0 <- solve(matrix(c(v_prior_alpha, cov_prior_alpha_beta, cov_prior_alpha_beta, v_prior_beta), 2, 2) / e_prior_Sigma2)
  a_0 <- a_prior_Sigma2
  b_0 <- b_prior_Sigma2

  # Auxiliary regression posterior
  X <- cbind(1, cbind(head(y, -1)))
  y <- tail(y, -1)
  if (NROW(y) == 1) {
    y <- t(y)
  } else {
    y <- cbind(y)
  }
  lambda_n <- crossprod(X) + lambda_0
  mu_n <- cbind(solve(lambda_n, lambda_0 %*% mu_0 + crossprod(X, y)))
  a_n <- a_0 + 0.5 * length(y)
  b_n <- b_0 + 0.5 * (crossprod(y) + crossprod(mu_0, lambda_0 %*% mu_0) - crossprod(mu_n, lambda_n %*% mu_n))

  # Compute posterior mean of SV parameters
  e_beta <- mu_n[2]
  e_Sigma2 <- c(b_n / (a_n - 1))
  v_alpha_beta <- e_Sigma2 * solve(lambda_n)
  ## mu (with Monte Carlo)
  chol_alpha_beta <- chol(v_alpha_beta)  # upper triangular
  n_variates <- 100000L
  variates_alpha_beta <- t(apply(matrix(rnorm(2 * n_variates), ncol = 2) %*% chol_alpha_beta, 1, function (x, mu) x + mu, mu_n))
  e_mu <- mean(variates_alpha_beta[, 1] / (1 - variates_alpha_beta[, 2]))
  # phi
  e_phi <- e_beta
  # sigma2
  v_beta <- v_alpha_beta[2, 2]
  e_sigma2 <- e_Sigma2 - v_normal_approx * (1 + v_beta + e_beta^2)

  c(mu = e_mu, phi = max(c(-0.999, min(c(0.999, e_phi)))), sigma2 = max(c(0.001, e_sigma2)))
}
## 3. Non-vanilla SV parameters (nu, rho)
### Prior mean should work just fine
init_sv_extra <- function (priorspec) {
  stop("nyi")
}
# TODO document and export
init_parameters <- function (y, X, priorspec, offset = 0) {
  beta_hat <- init_beta_regression(y, X, priorspec)
  vanilla_params <- init_sv_aux_regression(log((y - X %*% beta_hat)^2 + offset), priorspec)
  extra_params <- init_sv_extra(priorspec)
  list(beta = beta_hat, para = c(vanilla_params, extra_params))
}

asisprint <- function (x, censtring) {
  if (length(x) %% 2 == 0) {
    toshorten <- rep(as.character(censtring), length(x)/2)
    if (identical(toshorten, as.character(x)))
      return(sprintf("ASISx%d", length(x)/2))
  }
  paste0("(", paste(x, collapse=", "), ")")
}
