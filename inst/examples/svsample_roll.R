\donttest{
# Simulate from the true model
sim <- svsim(200)

# Perform rolling estimation using the vanilla SV
# model and default priors
roll <- svsample_roll(sim, draws = 5000, burnin = 2000,
                      keep_draws = TRUE,
                      forecast_length = 10,
                      n_ahead = 1, refit_every = 1,
                      refit_window = "moving",
                      calculate_predictive_likelihood = TRUE,
                      calculate_quantile = c(0.01, 0.05))

# Perform rolling estimation by making use
# of two CPU cores, advanced priors, and multiple
# chains with pre-set initial values. Let us combine
# that with an AR(2) specification
prior_beta <- sv_multinormal(c(1,0,-1), rbind(c(1, 0, 0.1),
                                              c(0, 0.3, -0.04),
                                              c(0.1, -0.04, 0.1)))
priorspec <- specify_priors(rho = sv_beta(4, 4),
                            latent0_variance = sv_constant(1),
                            beta = prior_beta,
                            nu = sv_exponential(0.05))
startpara <- list(list(mu = -9, phi = 0.3),
                  list(mu = -11, sigma = 0.1, phi = 0.95),
                  list(phi = 0.99))
roll <- svsample_roll(sim, draws = 5000, burnin = 2000,
                      designmatrix = "ar2",
                      priorspec = priorspec,
                      startpara = startpara,
                      parallel = "snow", n_cpus = 2,
                      n_chains = 3,
                      keep_draws = TRUE,
                      forecast_length = 10,
                      n_ahead = 1, refit_every = 1,
                      refit_window = "expanding",
                      calculate_predictive_likelihood = TRUE,
                      calculate_quantile = c(0.01, 0.05))
}
