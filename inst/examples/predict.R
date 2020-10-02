# Example 1
## Simulate a short and highly persistent SV process 
sim <- svsim(100, mu = -10, phi = 0.99, sigma = 0.2)

## Obtain 5000 draws from the sampler (that's not a lot)
draws <- svsample(sim$y, draws = 5000, burnin = 100,
  priormu = c(-10, 1), priorphi = c(20, 1.5), priorsigma = 0.2)

## Predict 10 days ahead
fore <- predict(draws, 10)

## Check out the results
summary(fore$h)
summary(fore$y)
plot(draws, forecast = fore)


# Example 2
## Simulate now an SV process with an AR(1) mean structure
len <- 109L
simar <- svsim(len, phi = 0.93, sigma = 0.15, mu = -9)
for (i in 2:len) {
  simar$y[i] <- 0.1 - 0.7 * simar$y[i-1] + simar$vol[i] * rnorm(1)
}

## Obtain 7000 draws
drawsar <- svsample(simar$y, draws = 7000, burnin = 300,
  designmatrix = "ar1", priormu = c(-10, 1), priorphi = c(20, 1.5),
  priorsigma = 0.2)

## Predict 7 days ahead (using AR(1) mean for the returns)
forear <- predict(drawsar, 7)

## Check out the results
plot(forear)
plot(drawsar, forecast = forear)

\dontrun{
# Example 3
## Simulate now an SV process with leverage and with non-zero mean
len <- 96L
regressors <- cbind(rep_len(1, len), rgamma(len, 0.5, 0.25))
betas <- rbind(-1.1, 2)
simreg <- svsim(len, rho = -0.42)
simreg$y <- simreg$y + as.numeric(regressors %*% betas)

## Obtain 12000 draws
drawsreg <- svsample(simreg$y, draws = 12000, burnin = 3000,
  designmatrix = regressors, priormu = c(-10, 1), priorphi = c(20, 1.5),
  priorsigma = 0.2, priorrho = c(4, 4))

## Predict 5 days ahead using new regressors
predlen <- 5L
predregressors <- cbind(rep_len(1, predlen), rgamma(predlen, 0.5, 0.25))
forereg <- predict(drawsreg, predlen, predregressors)

## Check out the results
summary(forereg$h)
summary(forereg$y)
plot(forereg)
plot(drawsreg, forecast = forereg)
}

