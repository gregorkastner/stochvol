## Simulate a highly persistent SV process 
sim <- svsim(500, mu = -10, phi = 0.99, sigma = 0.2)
 
## Obtain 4000 draws from the sampler (that's too few!)
draws <- svsample(sim$y, draws = 4000, burnin = 100, priormu = c(-10, 1),
                 priorphi = c(20, 1.2), priorsigma = 0.2)
 
## Predict 20 days ahead
fore <- predict(draws, 20)
 
## plot the results
plot(draws, forecast = fore)
 
\dontrun{
## Simulate an SV process with leverage
sim <- svsim(500, mu = -10, phi = 0.95, sigma = 0.2, rho=-0.5)
 
## Obtain 8000 draws from the sampler (that's too little!)
draws <- svsample(sim$y, draws = 4000, burnin = 3000, priormu = c(-10, 1),
                  priorphi = c(20, 1.2), priorsigma = 0.2,
                  priorrho = c(1, 1))
 
## Predict 20 days ahead
fore <- predict(draws, 20)
 
## plot the results
plot(draws, forecast = fore)
}