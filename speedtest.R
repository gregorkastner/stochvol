zz <- file("speedtest_output.txt", open = "wt")
sink(zz)
sink(zz, type = "message")
cat("#####\n", date(), "\n")
sessionInfo()

cat("### Some benchmark tests ###\nSource: https://mac.r-project.org/benchmarks/bench.R\n\nBEGIN\n")

# BEGIN benchmark tests
hilbert<-function(n) 1/(outer(seq(n),seq(n),"+")-1)
print("hilbert n=500")
print(system.time(eigen(hilbert(500))))
print(system.time(eigen(hilbert(500))))
print(system.time(eigen(hilbert(500))))
print("hilbert n=1000")
print(system.time(eigen(hilbert(1000))))
print(system.time(eigen(hilbert(1000))))
print(system.time(eigen(hilbert(1000))))
print("sort n=6")
print(system.time(sort(rnorm(10^6))))
print(system.time(sort(rnorm(10^6))))
print(system.time(sort(rnorm(10^6))))
print("sort n=7")
print(system.time(sort(rnorm(10^7))))
print(system.time(sort(rnorm(10^7))))
print(system.time(sort(rnorm(10^7))))
# loess
loess.me<-function(n) {
print(paste("loess n=",as.character(n),sep=""))
for (i in 1:5) {
    x<-rnorm(10^n); y<-rnorm(10^n); z<-rnorm(10^n)
    print(system.time(loess(z~x+y)))
    }
}
loess.me(3)
loess.me(4)
# END benchmark tests

cat("END\n\n#####\n\nSTART stochvol speedtest\n\n")

devtools::load_all(".")
set.seed(19891109)
dat <- list(
  normsmall = rnorm(20),
  normbig = rnorm(2000),
  svlsmall = svsim(len = 20, mu = -9, phi = 0.95, sigma = 0.1, rho = 0, nu = Inf)$y,
  svlbig = svsim(len = 2000, mu = -9, phi = 0.95, sigma = 0.1, rho = -0.3, nu = Inf)$y
)

for (dataset in names(dat)) {
  tmp <- replicate(2, svsample(dat[[dataset]], draws = 10000, burnin = 1000), simplify = FALSE)
  tmp <- replicate(2, svtsample(dat[[dataset]], draws = 10000, burnin = 1000), simplify = FALSE)
  tmp <- replicate(2, svlsample(dat[[dataset]], draws = 10000, burnin = 1000,
                                priormu = c(-9, 1),  # stability
                                priorphi = c(10, 10), priorsigma = 0.01, priorrho = c(10, 10),
                                expert = list(#parameterization = c("centered", "noncentered"),
                                              mhcontrol = list(use.mala = FALSE),
                                              correct.latent.draws = FALSE)),
                   simplify = FALSE)
}

cat("END stochvol speedtest\n\n")
sink(type = "message")
sink()


# some code for testing (AWOL produced NaNs for Darjus, numerically unstable)
if (FALSE) {
  dataset <- "normsmall"
  set.seed(3)
  tmp <- svlsample(dat[[dataset]], draws = 182, burnin = 1000,
                   #priormu = c(-9, 1),
                                expert = list(#parameterization = c("centered", "noncentered"),
                                              mhcontrol = list(use.mala = FALSE),
                                              correct.latent.draws = FALSE))
  plot(tmp)
  plot(tmp$adaptation_centered[, 1], type = "l")
  plot(tmp$adaptation_centered[, 2], type = "l")
  plot(tmp$adaptation_centered[, 3], type = "l")
  plot(tmp$adaptation_noncentered[, 1], type = "l")
  plot(tmp$adaptation_noncentered[, 2], type = "l")
  plot(tmp$adaptation_noncentered[, 3], type = "l")
}

