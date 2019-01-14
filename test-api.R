mydevlib <- "~/R/under_development"
library(stochvol, lib.loc = mydevlib)

set.seed(421)

# simulation
len <- 100L
dat <- list()
## svsim
### without extras
dat <- c(dat, list(svsim(len)))
dat <- c(dat, list(svsim(len, phi=0.9, mu=-10, sigma=0.1)))
### with nu
dat <- c(dat, list(svsim(len, phi=0.9, mu=-10, sigma=0.1, nu=4)))
### with rho
dat <- c(dat, list(svsim(len, phi=0.9, mu=-10, sigma=0.1, rho=-0.4)))
### with rho and nu (not allowed)
iserr <- try(svsim(len, phi=0.9, mu=-10, sigma=0.1, rho=-0.4, nu=4), silent=TRUE)
if (!inherits(iserr, "try-error")) warning("should've gotten an error")
### with mean
dmat <- cbind(rep_len(1, len), rgamma(len, 0.5, 0.25))
datnew <- svsim(len, phi=0.9, mu=-10, sigma=0.1, rho=-0.4)
datnew$y <- datnew$y + as.numeric(dmat %*% rbind(-1, 2))
dat <- c(dat, list(datnew))
### with AR(1) structure
datnew <- svsim(len, phi=0.9, mu=-10, sigma=0.1, nu=4)
datnew$y <- cumsum(datnew$y)
dat <- c(dat, list(datnew))

# sampling
res <- list()
res2 <- list()
## svsample
### normal error
#### vanilla
res <- c(res, lapply(dat, function (x) {svsample(x$y, draws = 200, quiet = TRUE)}))
res <- c(res, lapply(dat, function (x) {svsample(x$y, draws = 200, priorlatent0 = 2, quiet = TRUE)}))
res <- c(res, lapply(dat, function (x) {svsample(x$y, draws = 200, expert = list(gammaprior=FALSE), quiet = TRUE)}))
#### designmatrix
res <- c(res, lapply(dat, function (x, dm) {svsample(x$y, draws = 180, burnin = 1000, designmatrix = dm,
                     priormu = c(0, 1000), priorphi = c(20, 1.5), priorsigma = 0.2, quiet = TRUE,
                     priorbeta = c(0, 100), thinpara = 1, thinlatent = 2)}, dmat))
res <- c(res, lapply(dat, function (x) {svsample(x$y, draws = 180, burnin = 1000, designmatrix = "ar1",
                     priormu = c(0, 1000), priorphi = c(20, 1.5), priorsigma = 0.2, quiet = TRUE,
                     priorbeta = c(0, 100), thinpara = 1, thinlatent = 1)}))
### t error
#### vanilla
res <- c(res, lapply(dat, function (x) {svsample(x$y, draws = 200, priornu = c(2, 8), quiet = TRUE, thinpara = 4)}))
#### designmatrix
res <- c(res, lapply(dat, function (x, dm) {svsample(x$y, draws = 180, burnin = 1000, designmatrix = dm, quiet = TRUE,
                     priormu = c(0, 1000), priorphi = c(20, 1.5), priorsigma = 0.2,
                     priornu = c(2, 8), priorbeta = c(0, 100),
                     thinpara = 2, thinlatent = 2)}, dmat))
res <- c(res, lapply(dat, function (x) {svsample(x$y, draws = 280, burnin = 700, designmatrix = "ar1", quiet = TRUE,
                     priormu = c(0, 1000), priorphi = c(20, 3.5), priorsigma = 5,
                     priornu = c(2, 8), priorbeta = c(0, 100),
                     thinpara = 3, thinlatent = 1)}))
## svsample2
res2 <- c(res2, lapply(dat, function (x) {svsample2(x$y, draws = 160, burnin = 100,
                     priormu = c(0, 1000), priorphi = c(20, 1.5), priorsigma = 0.2, quiet = TRUE,
                     startlatent = rep_len(-9, length(x$y)), startpara = list(phi=-0.1, mu=2, sigma=5))}))
## svlsample
#### vanilla
res <- c(res, lapply(dat, function (x) {svlsample(x$y, draws = 200, quiet = TRUE)}))
res <- c(res, lapply(dat, function (x) {svlsample(x$y, draws = 200, expert = list(gammaprior=FALSE), quiet = TRUE)}))
#### designmatrix
res <- c(res, lapply(dat, function (x, dm) {svlsample(x$y, draws = 180, burnin = 1000, designmatrix = dm,
                     priormu = c(0, 1000), priorphi = c(20, 1.5), priorsigma = 0.2,
                     priorrho = c(3, 5), quiet = TRUE,
                     priorbeta = c(0, 100), thinpara = 1, thinlatent = 2)}, dmat))
res <- c(res, lapply(dat, function (x) {svlsample(x$y, draws = 180, burnin = 1000, designmatrix = "ar1",
                     priormu = c(0, 1000), priorphi = c(20, 1.5), priorsigma = 0.2,
                     priorrho = c(3, 9), quiet = TRUE,
                     priorbeta = c(0, 100), thinpara = 1, thinlatent = 1)}))
## svlsample2
res2 <- c(res2, lapply(dat, function (x) {svlsample2(x$y, draws = 160, burnin = 100,
                     priormu = c(0, 1000), priorphi = c(20, 1.5), priorsigma = 0.2, quiet = TRUE,
                     startlatent = rep_len(-9, length(x$y)), startpara = list(phi=-0.1, mu=2, sigma=5, rho=-0.4))}))
# prediction
pred <- lapply(res, function (x) predict(x, steps = 10L))
# plotting
png("todelete%03d.png", width = 240, height = 240, pointsize = 6)
lapply(res, plot)
dev.off()
# summary
summaries <- lapply(res, summary)
# residuals
resids <- lapply(res, residuals)
# print function
trash <- capture.output(lapply(dat, print))
#trash <- capture.output(lapply(res, print))  # really long time
# helpers
