#' Main function
#' 
#' @export
svlsample_r <- function (n=10000, y, thin=1, burnin=1000, init=NULL, prior=NULL,
                       param.sampler=c("rwMH", "constant", "auxiliary"),  # "MALA", "HMC"
                       latent.sampler=c("auxiliary", "auxiliaryMH", "constant"),  # "auxiliaryMH", "particle.filter"
                       sampler.params=NULL,
                       seed=NULL, strategy=c("centered", "non-centered")) {
  
  cParamSampler <- match.arg(param.sampler)
  cLatentSampler <- match.arg(latent.sampler)
  cStrategy <- match.arg(strategy, c("centered", "non-centered"), several.ok = TRUE)
  cSamplerParams <- sampler.params
  if (is.numeric(seed)) set.seed(seed)
  
  vY0 <- y
  
  if (!is.list(init)) init <- list()
  dPhi0 <- if (is.null(init$phi)) 0.9 else init$phi
  dRho0 <- if (is.null(init$rho)) -0.45 else init$rho
  dSigma20 <- if (is.null(init$sigma2)) 0.05 else init$sigma2
  dMu0 <- if (is.null(init$mu)) -10 else init$mu
  vVol0 <- if (is.null(init$vol)) rep_len(exp(dMu0/2), length(vY0)) else init$vol
  
  if (!is.list(prior)) prior <- list()
  cPriorPhi <- if (is.null(prior$phi)) c(20, 1.5) else prior$phi
  cPriorRho <- if (is.null(prior$rho)) c(3, 5) else prior$rho
  cPriorSigma2 <- if (is.null(prior$sigma2)) c(0.5, 0.5) else prior$sigma2
  cPriorMu <- if (is.null(prior$mu)) c(-10, 10) else prior$mu
  
  vY <- vY0 - mean(vY0)
  vYStar <- log(vY^2+1e-7)
  vD <- 2*(vY > 0) - 1
  lData <- list(y = vY, y.star = vYStar, d = vD)
  vH0 <- 2*log(vVol0)
  
  vVolIndices <- seq_along(vY)
  mParams <- matrix(NA_real_, floor(n/thin), 4)
  colnames(mParams) <- c("phi", "rho", "sigma2", "mu")
  mVols <- matrix(NA_real_, floor(n/thin), length(vVolIndices))
  .GlobalEnv$a <- rep_len(NA_real_, length(vY))
  .GlobalEnv$b <- rep_len(NA_real_, length(vY))
  .GlobalEnv$m <- rep_len(NA_real_, length(vY))
  .GlobalEnv$p <- rep_len(NA_real_, length(vY))
  .GlobalEnv$v <- rep_len(NA_real_, length(vY))
  .GlobalEnv$mspd <- matrix(NA_real_, length(vY), 10)  # mixture state posterior distribution
  
  bRedrawHWhenThinning <- TRUE  # has only effect for thinning > 0. if TRUE, then we always redraw H, otherwise we only redraw theta in the thinning rounds
  dPhi <- dPhi0; dRho <- dRho0; dSigma2 <- dSigma20; dMu <- dMu0
  lLatent <- list(h = vH0, ht = (vH0-dMu0)/sqrt(dSigma20), s = NULL)
  vH <- vH0
  vHt <- (vH0-dMu0)/sqrt(dSigma20)
  oRuntime <- system.time({
    for (i in seq(-burnin+1, n)) {
      bThinningRound <- (thin > 1) && (i %% thin != 0)  # is this a thinning round?
      ipburnin <- i+burnin
      npburnin <- n+burnin
      if (ipburnin %% round(npburnin/10) == 0) {
        if (ipburnin <= burnin) {
          cat("Burnin:   ")
        } else {
          cat("Sampling: ")
        }
        cat(round(ipburnin/npburnin*100), "%\n", sep="")
      }
      
      if (!bThinningRound || bRedrawHWhenThinning) {
        lLatent <- fnDrawLatent(lData, lLatent, dPhi, dRho, dSigma2, dMu, cStrategy[1], cLatentSampler)
      }
      
      for (cs in cStrategy) {
        vTheta <- fnDrawParameters(dPhi, dRho, dSigma2, dMu, lData, lLatent,
                                   cPriorPhi, cPriorRho, cPriorSigma2, cPriorMu,
                                   cs, cParamSampler, cSamplerParams)
        dPhi <- vTheta[1]
        dRho <- vTheta[2]
        dSigma2 <- vTheta[3]
        dMu <- vTheta[4]
        if (cs == "centered") {
          lLatent$ht <- (lLatent$h-dMu)/sqrt(dSigma2)
        } else if (cs == "non-centered") {
          lLatent$h <- sqrt(dSigma2)*lLatent$ht + dMu
        } else {
          stop("Invalid strategy")
        }
      }
      
      if ((i >= 1) && !bThinningRound) {
        mParams[i/thin, ] <- vTheta
        mVols[i/thin, ] <- exp(lLatent$h[vVolIndices]/2)
      }
    }
  }, gcFirst = FALSE)
  invisible(list(param = mParams, vol = mVols, runtime = oRuntime))
}

fnDrawLatent <- function (lData, lLatent, dPhi, dRho, dSigma2, dMu, cStrategy, cLatentSampler) {
  if (cLatentSampler == "auxiliary") {
    if (cStrategy == "centered") {
      vS <- draw_s_auxiliary(lData$y.star, lData$d, lLatent$h, dPhi, dRho, dSigma2, dMu, cStrategy, dfModelConstants)
      vH <- draw_h_auxiliary(lData$y.star, lData$d, vS, dPhi, dRho, dSigma2, dMu, cStrategy, dfModelConstants)
      vHt <- (vH-dMu)/sqrt(dSigma2)
    } else if (cStrategy == "non-centered") {
      vS <- draw_s_auxiliary(lData$y.star, lData$d, lLatent$ht, dPhi, dRho, dSigma2, dMu, cStrategy, dfModelConstants)
      vHt <- draw_h_auxiliary(lData$y.star, lData$d, vS, dPhi, dRho, dSigma2, dMu, cStrategy, dfModelConstants)
      vH <- sqrt(dSigma2)*vHt + dMu
    } else {
      stop("Invalid strategy")
    }
    list(h = vH, ht = vHt, s = vS)
  } else if (cLatentSampler == "auxiliaryMH") {
    vS <- NULL
    if (cStrategy == "centered") {
      vH <- draw_latent_auxiliaryMH(lData$y, lData$y.star, lData$d, lLatent$h, dPhi, dRho, dSigma2, dMu, dfModelConstants)
      vHt <- (vH-dMu)/sqrt(dSigma2)
    } else if (cStrategy == "non-centered") {
      vH <- draw_latent_auxiliaryMH(lData$y, lData$y.star, lData$d, lLatent$h, dPhi, dRho, dSigma2, dMu, dfModelConstants)
      vHt <- (vH-dMu)/sqrt(dSigma2)
      #vHt <- draw_latent_auxiliaryMH(lData$y, lData$y.star, lData$d, lLatent$ht, vS, dPhi, dRho, dSigma2, dMu, dfModelConstants)
      #vH <- sqrt(dSigma2)*vHt + dMu
    } else {
      stop("Invalid strategy")
    }
    list(h = vH, ht = vHt, s = vS)
  } else if (cLatentSampler == "constant") {  # mixture states are still drawn
    if (cStrategy == "centered") {
      vS <- draw_s_auxiliary(lData$y.star, lData$d, lLatent$h, dPhi, dRho, dSigma2, dMu, cStrategy, dfModelConstants)
    } else if (cStrategy == "non-centered") {
      vS <- draw_s_auxiliary(lData$y.star, lData$d, lLatent$ht, dPhi, dRho, dSigma2, dMu, cStrategy, dfModelConstants)
    } else {
      stop("Invalid strategy")
    }
    list(h = lLatent$h, ht = lLatent$ht, s = vS)
  } else {
    stop("Invalid latent sampler")
  }
}

fnDrawParameters <- function (dPhi, dRho, dSigma2, dMu, lData, lLatent,
                              cPriorPhi, cPriorRho, cPriorSigma2, cPriorMu,
                              cStrategy, cParamSampler, cSamplerParams) {
  if (cParamSampler == "rwMH") {
    if (cStrategy == "centered") {
      vTheta <- draw_theta_rwMH(dPhi, dRho, dSigma2, dMu, lData$y, lLatent$h,
                                cPriorPhi, cPriorRho, cPriorSigma2, cPriorMu,
                                cStrategy, cSamplerParams$stdev)
    } else if (cStrategy == "non-centered") {
      vTheta <- draw_theta_rwMH(dPhi, dRho, dSigma2, dMu, lData$y, lLatent$ht,
                                cPriorPhi, cPriorRho, cPriorSigma2, cPriorMu,
                                cStrategy, cSamplerParams$stdev)
    } else {
      stop("Invalid strategy")
    }
  } else if (cParamSampler == "constant") {
    c(dPhi, dRho, dSigma2, dMu)
  } else if (cParamSampler == "auxiliary") {
    vTheta <- draw_theta_auxiliary(dPhi, dRho, dSigma2, dMu, lData$y.star, lData$d, lLatent$s,
                                   cPriorPhi, cPriorRho, cPriorSigma2, cPriorMu, dfModelConstants)
  } else {
    stop("Invalid parameter sampler")
  }
}
