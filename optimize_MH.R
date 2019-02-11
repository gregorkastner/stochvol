options(warn = 1)

library(xts)
library(stochvol, lib.loc = paste0(Sys.getenv("HOME"), "/R/under_development"))

dat <- readRDS("dat.RDS")

version <- "1"
iTaskID <- as.integer(Sys.getenv("SGE_TASK_ID", 2204))
cat("Task ID: ", iTaskID, "\n")

vDataIdx <- seq_along(dat)
vMhgrid <- exp(seq(log(0.001), log(1), length.out=10L))
vRgrid <- exp(seq(log(0.001), log(2), length.out=10L))
vCentering <- c("centered", "noncentered")
dfGrid <- expand.grid(iDataIdx = vDataIdx,
                      dMhgrid = vMhgrid,
                      dRgrid = vRgrid,
                      sCentering = vCentering,
                      stringsAsFactors = FALSE)
iRow <- iTaskID
lParams <- as.list(dfGrid[iRow, ])

iDataIdx <- lParams$iDataIdx
dMhgrid <- lParams$dMhgrid
dRgrid <- lParams$dRgrid
sCentering <- lParams$sCentering

# get data
dData <- as.numeric(tail(diff(log(dat[[iDataIdx]])), -1))
dData <- dData - mean(dData)
sDataName <- names(dat)[iDataIdx]

# set parameters
dSigma2Prior <- 1
dMuMean <- -10
dMuSigma <- 10
dPhiA <- 20
dPhiB <- 1.5
dRhoA <- 4
dRhoB <- 4
iNsim <- 20000L
iBurnin <- 20000L
iSvburnin <- 10000L

# get random seed
rng <- file("/dev/urandom","rb") # open connection
iSeed <- readBin(rng,what="integer",n=1) # read some 8-byte integers 
close(rng) # close the connection
iSeed <- abs(as.integer(Sys.time()) %% 1000003L + iSeed %% 1000003L) %% 1000003L
print(iSeed)

# run sampler
lResult <- svlsample(y = dData, draws = iNsim, burnin = iBurnin,
                     designmatrix = NA, priormu = c(dMuMean, dMuSigma),
                     priorphi = c(dPhiA, dPhiB), priorsigma = dSigma2Prior,
                     priorrho = c(dRhoA, dRhoB),
                     expert = list(mhcontrol = dMhgrid,
                                   init.with.svsample = iSvburnin,
                                   parameterization = sCentering,
                                   var.rho = dRgrid))

lToSave <- list(data.ts = dData,
                data.name = sDataName,
                post.summary = summary(lResult),
                runtime = lResult$runtime,
                seed = iSeed,
                accept.rate.param = lResult$accept.phi/iNsim,
                accept.rate.latent = lResult$accept.h/iNsim)

new.folder <- paste0("results", version, "/", iTaskID)
dir.create(new.folder, recursive = TRUE, showWarnings = FALSE)
setwd(new.folder)
saveRDS(lToSave, "result.RDS")

# plots
cat("Plots\n", file=stderr())

pdf("posterior.pdf", width=8, height=8)
plot(lResult)
dev.off()

pdf("data.pdf", width=8, height=5)
ts.plot(dData)
dev.off()

