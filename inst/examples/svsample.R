###############
# Full examples
###############

# Example 1
## Simulate a short and highly persistent SV process 
sim <- svsim(100, mu = -10, phi = 0.99, sigma = 0.2)

## Obtain 5000 draws from the sampler (that's not a lot)
draws <-
  svsample(sim, draws = 5000, burnin = 100,
           priormu = c(-10, 1), priorphi = c(20, 1.5), priorsigma = 0.2)

## Check out the results
summary(draws)
plot(draws)


\donttest{
# Example 2
## Simulate an asymmetric and conditionally heavy-tailed SV process
sim <- svsim(150, mu = -10, phi = 0.96, sigma = 0.3, nu = 10, rho = -0.3)

## Obtain 10000 draws from the sampler
## Use more advanced prior settings
## Run two parallel MCMC chains
advanced_draws <-
  svsample(sim, draws = 10000, burnin = 5000,
           priorspec = specify_priors(mu = sv_normal(-10, 1),
                                      sigma2 = sv_gamma(0.5, 2),
                                      rho = sv_beta(4, 4),
                                      nu = sv_constant(5)),
           parallel = "snow", n_chains = 2, n_cpus = 2)

## Check out the results
summary(advanced_draws)
plot(advanced_draws)


# Example 3
## AR(1) structure for the mean
data(exrates)
len <- 3000
ahead <- 100
y <- head(exrates$USD, len)

## Fit AR(1)-SVL model to EUR-USD exchange rates
res <- svsample(y, designmatrix = "ar1")

## Use predict.svdraws to obtain predictive distributions
preddraws <- predict(res, steps = ahead)

## Calculate predictive quantiles
predquants <- apply(predy(preddraws), 2, quantile, c(.1, .5, .9))

## Visualize
expost <- tail(head(exrates$USD, len+ahead), ahead)
ts.plot(y, xlim = c(length(y)-4*ahead, length(y)+ahead),
	       ylim = range(c(predquants, expost, tail(y, 4*ahead))))
for (i in 1:3) {
  lines((length(y)+1):(length(y)+ahead), predquants[i,],
        col = 3, lty = c(2, 1, 2)[i])
}
lines((length(y)+1):(length(y)+ahead), expost,
      col = 2)


# Example 4
## Predicting USD based on JPY and GBP in the mean
data(exrates)
len <- 3000
ahead <- 30
## Calculate log-returns
logreturns <- apply(exrates[, c("USD", "JPY", "GBP")], 2,
                    function (x) diff(log(x)))
logretUSD <- logreturns[2:(len+1), "USD"]
regressors <- cbind(1, as.matrix(logreturns[1:len, ]))  # lagged by 1 day

## Fit SV model to EUR-USD exchange rates
res <- svsample(logretUSD, designmatrix = regressors)

## Use predict.svdraws to obtain predictive distributions
predregressors <- cbind(1, as.matrix(logreturns[(len+1):(len+ahead), ]))
preddraws <- predict(res, steps = ahead,
                     newdata = predregressors)
predprice <- exrates[len+2, "USD"] * exp(t(apply(predy(preddraws), 1, cumsum)))

## Calculate predictive quantiles
predquants <- apply(predprice, 2, quantile, c(.1, .5, .9))

## Visualize
priceUSD <- exrates[3:(len+2), "USD"]
expost <- exrates[(len+3):(len+ahead+2), "USD"]
ts.plot(priceUSD, xlim = c(len-4*ahead, len+ahead+1),
	       ylim = range(c(expost, predquants, tail(priceUSD, 4*ahead))))
for (i in 1:3) {
  lines(len:(len+ahead), c(tail(priceUSD, 1), predquants[i,]),
        col = 3, lty = c(2, 1, 2)[i])
}
lines(len:(len+ahead), c(tail(priceUSD, 1), expost),
      col = 2)
}

########################
# Further short examples
########################

y <- svsim(50, nu = 10, rho = -0.1)$y

# Supply initial values
res <- svsample(y,
                startpara = list(mu = -10, sigma = 1))

\donttest{
# Supply initial values for parallel chains
res <- svsample(y,
                startpara = list(list(mu = -10, sigma = 1),
                                 list(mu = -11, sigma = .1, phi = 0.9),
                                 list(mu = -9, sigma = .3, phi = 0.7)),
                parallel = "snow", n_chains = 3, n_cpus = 2)

# Parallel chains with with a pre-defined cluster object
cl <- parallel::makeCluster(2, "PSOCK", outfile = NULL)  # print to console
res <- svsample(y,
                parallel = "snow", n_chains = 3, cl = cl)
parallel::stopCluster(cl)
}

# Turn on correction for model misspecification
## Since the approximate model is fast and it is working very
##   well in practice, this is turned off by default
res <- svsample(y,
                expert = list(correct_model_misspecification = TRUE))
print(res)

\dontrun{
# Parallel multicore chains (not available on Windows)
res <- svsample(y, draws = 30000, burnin = 10000,
                parallel = "multicore", n_chains = 3, n_cpus = 2)

# Plot using a color palette
palette(rainbow(coda::nchain(para(res, "all"))))
plot(res)

# Use functionality from package 'coda'
## E.g. Geweke's convergence diagnostics
coda::geweke.diag(para(res, "all")[, c("mu", "phi", "sigma")])

# Use functionality from package 'bayesplot'
bayesplot::mcmc_pairs(res, pars = c("sigma", "mu", "phi", "h_0", "h_15"))
}

