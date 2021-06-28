# Simulate data
n <- 50L
dat <- data.frame(x = runif(n, 3, 4),
                  z = runif(n, -1, -0.5))
designmatrix <- matrix(c(dat$x, dat$x^2, log10(dat$x),
                         dat$z), ncol = 4)
betas <- matrix(c(-1, 1, 2, 0), ncol = 1)
y <- designmatrix %*% betas + svsim(n)$y
dat$y <- y
# Formula interface
res <- svlm(y ~ 0 + x + I(x^2) + log10(x) + z, data = dat)
# Prediction
predn <- 10L
preddat <- data.frame(x = runif(predn, 3, 4),
                      z = runif(predn, -1, -0.5))
pred <- predict(res, newdata = preddat, steps = predn)
