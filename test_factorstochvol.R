# > example("fsvsample", package="factorstochvol")
library(factorstochvol)

# Load exchange rate data (ships with stochvol):
data(exrates, package = "stochvol")
exrates$date <- NULL

# Compute the de-meaned percentage log returns:
dat <- 100 * logret(exrates, demean = TRUE)

# We are going to fit a one-factor model so the ordering is irrelevant
# NOTE that these are very few draws, you probably want more...
set.seed(45)
run <- system.time(
  res <- fsvsample(dat, factors = 2, draws = 4, burnin = 0, runningstore = 6)
)
print(400/run["user.self"])
# 26.5 for original stochvol + factorstochvol
# 25.7 afterwards
# 27.2 after turning everything into arma
#saveRDS(res, "./test_factorstochvol_new.RDS")
saveRDS(res, "./test_factorstochvol_old.RDS")
res_new <- readRDS("./test_factorstochvol_new.RDS")
isTRUE(all.equal(res$f, res_new$f)) &&
  isTRUE(all.equal(res$para, res_new$para)) &&
  isTRUE(all.equal(res$h, res_new$h)) &&
  isTRUE(all.equal(res$facload, res_new$facload)) &&
  isTRUE(all.equal(res$h0, res_new$h0))

voltimeplot(res)

corimageplot(res, nrow(dat), plotCI = 'circle')

oldpar <- par(ask = TRUE)
plot(res)
par(oldpar)
