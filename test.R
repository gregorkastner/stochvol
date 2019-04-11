library(devtools)

load_all()

set.seed(5556)
sim <- svsim(500, mu = -10, phi = 0.9, sigma = 0.3, rho = -0.6)
draws1 <- svlsample(sim, burnin = 1000, draws = 50000,
                    priormu = c(-10, 0.0001), dontupdatemu = FALSE,
                    startpara = list(phi = 0.9, mu = -10, sigma = 0.3, rho = -0.6))
draws2 <- svlsample(sim, burnin = 1000, draws = 50000,
                    dontupdatemu = TRUE,
                    startpara = list(phi = 0.9, mu = -10, sigma = 0.3, rho = -0.6))

png("svlsample-no-mu%03d.png", width = 1024, height = 1024)
plot(draws1, simobj = sim)
plot(draws2, simobj = sim)
qqplot(draws1$para[, "phi"], draws2$para[, "phi"]); abline(a = 0, b = 1)
qqplot(draws1$para[, "sigma"], draws2$para[, "sigma"]); abline(a = 0, b = 1)
qqplot(draws1$para[, "rho"], draws2$para[, "rho"]); abline(a = 0, b = 1)
qqplot(draws1$latent[, 30], draws2$latent[, 30]); abline(a = 0, b = 1)
dev.off()

