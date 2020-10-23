# Simulate data
sim <- svsim(150)

# Draw from vanilla SV
draws <- svsample(sim, draws = 2000)

## Summarize all parameter draws as a merged mcmc object
summary(para(draws))
## Extract the draws as an mcmc.list object
para(draws, chain = "all")

\donttest{
## Further short examples
summary(latent0(draws))
summary(latent(draws))
summary(vola(draws))
sampled_parameters(draws)
priors(draws)

# Draw 3 independent chains from heavy-tailed and asymmetric SV with AR(2) structure
draws <- svsample(sim, draws = 20000, burnin = 3000,
                  designmatrix = "ar2",
                  priornu = 0.1, priorrho = c(4, 4),
                  n_chains = 3)

## Extract beta draws from the second chain
svbeta(draws, chain = 2)
## ... tau draws from all chains merged/concatenated together
svtau(draws)
## Create a new svdraws object from the first and third chain
second_chain_excluded <- draws[c(1, 3)]

# Draw from the predictive distribution
pred <- predict(draws, steps = 2)

## Extract the predicted observations as an mcmc.list object
predicted_y <- predy(pred, chain = "all")
## ... the predicted standard deviations from the second chain
predicted_sd <- predvola(pred, chain = 2)
## Create a new svpredict object from the first and third chain
second_chain_excluded <- pred[c(1, 3)]
}

