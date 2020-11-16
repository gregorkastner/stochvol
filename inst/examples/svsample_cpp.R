# Draw one sample using fast SV and general SV
y <- svsim(40)$y
params <- list(mu = -10, phi = 0.9, sigma = 0.1,
               nu = Inf, rho = 0, beta = NA,
               latent0 = -10)
res_fast <- svsample_fast_cpp(y,
  startpara = params, startlatent = rep(-10, 40))
res_gen <- svsample_general_cpp(y,
  startpara = params, startlatent = rep(-10, 40))

# Embed SV in another sampling scheme
## vanilla SV
len <- 40L
draws <- 1000L
burnin <- 200L
param_store <- matrix(NA, draws, 3,
                      dimnames = list(NULL,
                                      c("mu", "phi", "sigma")))
startpara <- list(mu = 0, phi = 0.9, sigma = 0.1,
                  nu = Inf, rho = 0, beta = NA,
                  latent0 = 0)
startlatent <- rep(0, len)
for (i in seq_len(burnin+draws)) {
  # draw the data in the bigger sampling scheme
  # now we simulate y from vanilla SV
  y <- svsim(len, mu = 0, phi = 0.9, sigma = 0.1)$y
  # call SV sampler
  res <- svsample_fast_cpp(y, startpara = startpara,
                           startlatent = startlatent)
  # administrate values
  startpara[c("mu","phi","sigma")] <-
    as.list(res$para[, c("mu", "phi", "sigma")])
  startlatent <- drop(res$latent)
  # store draws after the burnin
  if (i > burnin) {
    param_store[i-burnin, ] <-
      res$para[, c("mu", "phi", "sigma")]
  }
}
### quick look at the traceplots
ts.plot(param_store, col = 1:3)

