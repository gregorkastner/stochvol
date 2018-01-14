#' Samplers

#' Main MCMC sampler
#' 
#' Combines the samplers and features burn-in.
#' @author hdarjus \email{hdarjus@gmail.com}
#' @param zooY zoo object of observations
#' @param iNsim number of simulation steps
#' @param lPriors list containing the prior parameters, elements: \code{mu.mean, mu.var,
#'   phi.a, phi.b, sigma2.shape, sigma2.rate, rho.a, rho.b}
#' @param lInit list of initial values, elements: \code{phi, sigma2, rho, mu, h}
#' @param iBurnin (optional) number of steps ignore in the beginning
#' @return List of three elements: a dataframe of samples with columns \code{phi, sigma2, rho, mu},
#'   an \code{iNsim x length(zooY)} matrix of log volatility samples, and a vector of weights for reweighting.
#'   The returned samples are not reweighted.
#' @details The function uses the mixture constants in the background. Some recommended
#'   initial parameter values: phi = .95, sigma2 = .01, mu = 2*log(.01), rho = -.4, and also
#'   h_i = mu, s_i = 5.
#' @examples 
#'   require("zoo")
#'   zooData <- fnGenLogSV(200, -0.5, 0.1, 0.9, -9, iSeed = 42)
#'   lResult <- fnMCMCSampler(zooData[, "y"],
#'                            iNsim = 10,
#'                            lPriors = list(mu.mean = -10, mu.var = 10, phi.a = 20, phi.b = 1.5, sigma2.shape = 2.5, sigma2.rate = 0.025, rho.a = 1, rho.b = 1),
#'                            lInit = list(phi = 0.9, sigma2 = 0.03, rho = -0.9, mu = -9, h = zoo::coredata(zooData[, "h"])),
#'                            iBurnin = 0,
#'                            sSigma2.Prior = "inv.gamma",
#'                            sCentering = "centered")
#' @export
fnMCMCSampler <- function (zooY, iNsim, lPriors, lInit, iBurnin = 1000,
                           sSigma2.Prior = c("inv.gamma", "gamma"),
                           sCentering = c("centered", "non-centered", "GIS_C", "GIS_NC"),
                           bReturnWeights = FALSE) {
  # match input
  sSigma2.Prior <- match.arg(sSigma2.Prior)
  sCentering <- match.arg(sCentering)
  
  # data
  vY <- zoo::coredata(zooY)
  vY <- vY - mean(vY)  # demean
  iN <- length(vY)
  vYStar <- log(vY^2 + 1e-7)
  vD <- 2*(vY >= 0) - 1
  
  # set current values to initials
  dCurPhi <- lInit$phi
  dCurSigma2 <- lInit$sigma2
  dCurMu <- lInit$mu
  dCurRho <- lInit$rho
  # Nakajima uses initial values m_i = -1.2704, v_i = sqrt(4.93), a_i = 1, b_i = .5, I use s_i = 5. I only save s_i directly.
  vCurH <- if (is.null(lInit$h)) {
    rep(if (sCentering %in% c("centered", "GIS_C")) dCurMu else 0, iN)
  } else {
    lInit$h
  }
  vCurS <- NULL
  
  # container for the samples
  dfSamples <- data.frame(phi = rep(NA_real_, iNsim), sigma2 = NA_real_, rho = NA_real_, mu = NA_real_)
  mSampleVol <- matrix(NA_real_, nrow = iNsim, ncol = iN)
  
  # there's nothing to reuse in the first loop
  lReuseDiffOptim <- NULL
  lReuseFrequency <- 1

  # setup ASIS
  vCenteringOrder <- if (sCentering == "GIS_C") {
    c("centered", "non-centered")
  } else if (sCentering == "GIS_NC") {
    c("non-centered", "centered")
  } else {
    sCentering
  }
  
  # sampling loop
  for (iI in seq(-iBurnin+1, iNsim)) {
    if (iI %% 500 == 0) {
      cat(iI, "\n  phi   ", dCurPhi, "\n  sigma2 ", dCurSigma2, "\n  rho   ", dCurRho, "\n  mu    ", dCurMu, "\n")
      print(summary(vCurH))
      print(summary(vCurS))
    }
    
    # sample mixture state
    vCurS <- fnSampleMixtureState(vYStar, vCurH, vD, dCurMu, dCurSigma2, dCurRho, dCurPhi, vCenteringOrder[1])
    
    # sample theta = (phi, sigma2, rho)
    for (sC in rev(vCenteringOrder)) {
      lThetaResult <- fnSampleTheta(dCurPhi, dCurSigma2, dCurRho,
                                    vYStar, vD,
                                    vCurS,
                                    c(lPriors$mu.mean, lPriors$mu.var),
                                    c(lPriors$phi.a, lPriors$phi.b),
                                    c(lPriors$sigma2.shape, lPriors$sigma2.rate),
                                    c(lPriors$rho.a, lPriors$rho.b),
                                    lReuseDiffOptim = if (iI %% lReuseFrequency == 0) NULL else lReuseDiffOptim,
                                    sSigma2.Prior = sSigma2.Prior, sCentering = sC)
    
      dCurPhi <- lThetaResult$phi
      dCurSigma2 <- lThetaResult$sigma2
      dCurRho <- lThetaResult$rho
      lFilterResults <- lThetaResult$filter.results
      lReuseDiffOptim <- lThetaResult$reuse.diff.optim
    }
    
    # sample mu
    dCurMu <- fnSampleMu(dq = lFilterResults$q, dQ = lFilterResults$Q)
    
    # sample volatility
    vCurH <- fnSampleH(dCurPhi, dCurSigma2, dCurRho, dCurMu, vD, vCurS, lFilterResults, vCenteringOrder[1])
    
    # save samples
    if (iI >= 1) {
      dfSamples$phi[iI] <- dCurPhi
      dfSamples$sigma2[iI] <- dCurSigma2
      dfSamples$rho[iI] <- dCurRho
      dfSamples$mu[iI] <- dCurMu
      mSampleVol[iI, ] <- if (vCenteringOrder[1] == "centering") exp(vCurH/2) else exp(sqrt(dCurSigma2) * vCurH/2 + dCurMu)
    }
  }
  
  vWeights <- if (bReturnWeights) {
    fnCalcWeights(vYStar, vD, mSampleVol, dfSamples, vCenteringOrder[1])
  } else {
    NULL
  }
  
  return(list(samples = dfSamples,
              vol = mSampleVol,
              weights = vWeights))
}

#' Sampler for \eqn{\theta}
#' 
#' This is step 2a, to sample \eqn{\theta=(\rho,\sigma^2,\phi)} conditional
#' on \eqn{(y*, d, s)}.
#' @author hdarjus \email{hdarjus@gmail.com}
#' @param dPhiOld current parameter \eqn{\phi}
#' @param dSigma2Old current parameter \eqn{\sigma} squared
#' @param dRhoOld current parameter \eqn{\rho}
#' @param vYStar vector of transformed measurements
#' @param vD vector of signs
#' @param vS vector of current mixture states
#' @param vMuPriorParams vector of length 2, mean and variance of the normal prior of \eqn{\mu}
#' @param vPhiPriorParams vector of length 2, two parameters of the beta prior of \eqn{\phi}
#' @param vSigma2PriorParams vector of length 2, shape and rate of the gamma prior of inverse square \eqn{\sigma}
#' @param vRhoPriorParams vector of length 2, two parameters of the beta prior of \eqn{\rho}
#' @param lReuseDiffOptim list of cached values, elements: theta.star, theta.hat, log.post.hat, gradient, hess, sigma.star
#' @return List of 4 elements, \code{phi, sigma2, rho, filter.results}, the a posteriori sampled \eqn{\theta} and a list of the augmented Kalman filter by-products.
#' @seealso \code{fnAugKalmanFilter, fnAugKalmanFilterCpp}
fnSampleTheta <- function (dPhiOld, dSigma2Old, dRhoOld, vYStar, vD, vS, vMuPriorParams,
                           vPhiPriorParams, vSigma2PriorParams, vRhoPriorParams,
                           lReuseDiffOptim = NULL, sSigma2.Prior, sCentering) {
  va <- dfModelConstants$a[vS]
  vb <- dfModelConstants$b[vS]
  vm <- dfModelConstants$m[vS]
  vv <- dfModelConstants$v[vS]
  dFallbackVariance <- 0.25
  mThetaOld <- rbind(dPhiOld, dSigma2Old, dRhoOld)
  
  if (is.null(lReuseDiffOptim)) {
    # vX = c(log((1+dPhiOld)/(1-dPhiOld)), log(dSigma2Old), log((1+dRhoOld)/(1-dRhoOld)))
    fnLogPostMax <- function (vX, va, vb, vm, vv, vD, vYStar,
                              vMuPriorParams, vPhiPriorParams, vSigma2PriorParams, vRhoPriorParams,
                              sSigma2.Prior, sCentering) {
      #cat(".")
      fnThetaLogPostDistTransInput(vX[1], vX[2], vX[3],
                                   va, vb, vm, vv, vD, vYStar,
                                   vMuPriorParams, vPhiPriorParams, vSigma2PriorParams, vRhoPriorParams,
                                   sSigma2.Prior, sCentering)$result
    }
    # vX = c(dPhi, dSigma2, dRho)
    fnLogPostDiff <- function (vX, va, vb, vm, vv, vD, vYStar,
                               vMuPriorParams, vPhiPriorParams, vSigma2PriorParams, vRhoPriorParams,
                               sSigma2.Prior, sCentering) {
      #cat(";")
      fnThetaLogPostDist(vX[1], vX[2], vX[3],
                         va, vb, vm, vv, vD, vYStar,
                         vMuPriorParams, vPhiPriorParams, vSigma2PriorParams, vRhoPriorParams,
                         sSigma2.Prior, sCentering)$result
    }
    # Find \hat\theta for Eq. (7) in Nakajima
    vInit <- c(log((1+dPhiOld)/(1-dPhiOld)), log(dSigma2Old), log((1+dRhoOld)/(1-dRhoOld)))
    lOptimRes <- tryCatch(
      optim(par = vInit, fn = fnLogPostMax, gr = NULL,
            va = va, vb = vb, vm = vm, vv = vv, vD = vD, vYStar = vYStar,
            vMuPriorParams = vMuPriorParams,
            vPhiPriorParams = vPhiPriorParams,
            vSigma2PriorParams = vSigma2PriorParams,
            vRhoPriorParams = vRhoPriorParams,
            sSigma2.Prior = sSigma2.Prior,
            sCentering = sCentering,
            method = "BFGS", control = list(fnscale = -1, maxit = 6)),
      error = function (e) {  # trying a different optimizer when BFGS doesn't work
        warning(message(e))
        optim(par = vInit, fn = fnLogPostMax, gr = NULL,
              va = va, vb = vb, vm = vm, vv = vv, vD = vD, vYStar = vYStar,
              vMuPriorParams = vMuPriorParams,
              vPhiPriorParams = vPhiPriorParams,
              vSigma2PriorParams = vSigma2PriorParams,
              vRhoPriorParams = vRhoPriorParams,
              sSigma2.Prior = sSigma2.Prior,
              sCentering = sCentering,
              method = "Nelder-Mead", control = list(fnscale = -1, maxit = 100))
      })
    dPhiHat <- 2/(1+exp(-lOptimRes$par[1]))-1
    dSigma2Hat <- exp(lOptimRes$par[2])
    dRhoHat <- 2/(1+exp(-lOptimRes$par[3]))-1
    vThetaHat <- c(dPhiHat, dSigma2Hat, dRhoHat)
    dLogPostHat <- lOptimRes$value
    # Find \theta_\ast and \Sigma_\ast
    dDist = 0.0001
    vThetaHatRestricted <- c(max(min(dPhiHat, 1-dDist), -1+dDist),
                             max(dSigma2Hat, dDist),
                             max(min(dRhoHat, 1-dDist), -1+dDist))
    mGrad <- matrix(numDeriv::grad(fnLogPostDiff, vThetaHatRestricted,
                                   #method = "simple", method.args = list(eps=dDist/2),
                                   method = "Richardson", method.args = list(d=dDist/2, r=3),
                                   va = va, vb = vb, vm = vm, vv = vv, vD = vD, vYStar = vYStar,
                                   vMuPriorParams = vMuPriorParams,
                                   vPhiPriorParams = vPhiPriorParams,
                                   vSigma2PriorParams = vSigma2PriorParams,
                                   vRhoPriorParams = vRhoPriorParams,
                                   sSigma2.Prior = sSigma2.Prior,
                                   sCentering = sCentering),
                    ncol = 1)
    mHessian <- numDeriv::hessian(fnLogPostDiff, vThetaHatRestricted,
                                  method.args = list(d=dDist/2, r=4),
                                  va = va, vb = vb, vm = vm, vv = vv, vD = vD, vYStar = vYStar,
                                  vMuPriorParams = vMuPriorParams,
                                  vPhiPriorParams = vPhiPriorParams,
                                  vSigma2PriorParams = vSigma2PriorParams,
                                  vRhoPriorParams = vRhoPriorParams,
                                  sSigma2.Prior = sSigma2.Prior,
                                  sCentering = sCentering)
    mSigmaStar <- tryCatch(as.matrix(Matrix::nearPD(-solve(mHessian))$m),  # make sure positive definiteness
                           error = function (e) NULL)
    vEigenValues <- sort(-1/eigen(mHessian, symmetric = T, only.values = T)$values)
    bFallbackCovariance <-
      is.null(mSigmaStar) ||
      any(is.infinite(mHessian) | is.nan(mHessian)) ||  # no NaNs or Infs
      (kappa(mHessian, exact = TRUE) <= .Machine$double.eps) ||  # computationally singular
      (head(vEigenValues, 1) <= .Machine$double.eps) ||  # not positive definite
      (tail(vEigenValues, 1) >= Inf)  # too large eigenvalue => possible problems
    if (bFallbackCovariance) {
      warning("Using fallback covariance matrix")
      mSigmaStar <- diag(max(dFallbackVariance, norm(cbind(vThetaHat), "F")^2), nrow = nrow(mGrad))
      mGrad <- 0 * mGrad
    }
    mThetaStar <- vThetaHat + mSigmaStar %*% mGrad
    # cache values for later use
    lReuseNew <- list(
      theta.star = mThetaStar,
      theta.hat = vThetaHat,
      log.post.hat = dLogPostHat,
      gradient = mGrad,
      hess = mHessian,
      sigma.star = mSigmaStar
    )
  } else {
    lReuseNew <- lReuseDiffOptim
    mThetaStar <- lReuseDiffOptim$theta.star
    vThetaHat <- lReuseDiffOptim$theta.hat
    dLogPostHat <- lReuseDiffOptim$log.post.hat
    mGrad <- lReuseDiffOptim$gradient
    mHessian <- lReuseDiffOptim$hess
    mSigmaStar <- lReuseDiffOptim$sigma.star
  }
  # Generate \theta^\ast from the truncated normal
  bFound <- FALSE
  vThetaStar <- drop(mThetaStar)
  dEpsilonBound <- 0.0001
  mTruncNormalBounds <- rbind(lower = c(-1+dEpsilonBound, dEpsilonBound, -1+dEpsilonBound), upper = c(1-dEpsilonBound, Inf, 1-dEpsilonBound))
  while (!bFound) {
    dSupportProb <- mvtnorm::pmvnorm(lower = mTruncNormalBounds["lower", ], upper = mTruncNormalBounds["upper", ], mean = vThetaStar, sigma = mSigmaStar)
    mThetaNew <- if (dSupportProb <= 0.1) {
      tmvtnorm::rtmvnorm(1, mean = vThetaStar, sigma = mSigmaStar, lower = mTruncNormalBounds["lower", ], upper = mTruncNormalBounds["upper", ],
                         algorithm = "gibbs", burn.in.samples = 1000, start.value = c(0, 0.1, 0))
    } else {
      tmvtnorm::rtmvnorm(1, mean = vThetaStar, sigma = mSigmaStar, lower = mTruncNormalBounds["lower", ], upper = mTruncNormalBounds["upper", ],
                         algorithm = "rejection")
    }
    dim(mThetaNew) <- c(3, 1)  # transpose
    bFound <- !any(is.nan(mThetaNew) | is.na(mThetaNew))
    if (!bFound) {
      mSigmaStar <- diag(max(dFallbackVariance, norm(mThetaStar, "F")^2), 3)
    }
  }
  dPhiNew <- mThetaNew[1]
  dSigma2New <- mThetaNew[2]
  dRhoNew <- mThetaNew[3]
  # Intermediate results
  lIntResOld <- fnThetaLogPostDist(dPhiOld, dSigma2Old, dRhoOld,
                                   va, vb, vm, vv, vD, vYStar,
                                   vMuPriorParams, vPhiPriorParams, vSigma2PriorParams, vRhoPriorParams,
                                   sSigma2.Prior, sCentering)
  lIntResNew <- fnThetaLogPostDist(dPhiNew, dSigma2New, dRhoNew,
                                   va, vb, vm, vv, vD, vYStar,
                                   vMuPriorParams, vPhiPriorParams, vSigma2PriorParams, vRhoPriorParams,
                                   sSigma2.Prior, sCentering)
  # Metropolis-Hastings step
  # here Nakajima takes the second order Taylor expansion as the log density, which should be equivalent to this;
  # constants of the distributions are omitted since they cancel out;
  dLogPiOld <- lIntResOld$result
  dLogPiNew <- lIntResNew$result
  fnLogProposal <- function (mTheta) {
    return(-0.5 * t(mTheta - mThetaStar) %*% solve(mSigmaStar) %*% (mTheta - mThetaStar))
  }
  dLogProposalOld <- fnLogProposal(mThetaOld)
  dLogProposalNew <- fnLogProposal(mThetaNew)
  
  dAcceptance <- (dLogPiNew - dLogPiOld) - (dLogProposalNew - dLogProposalOld)
  if (dAcceptance >= log(runif(1, 0, 1))) {  # accept
    mThetaResult <- mThetaNew
    lFilterResults <- lIntResNew$filter.results
  } else {  # reject
    mThetaResult <- mThetaOld
    lFilterResults <- lIntResOld$filter.results
  }
  
  return(list(
    phi = mThetaResult[1],
    sigma2 = mThetaResult[2],
    rho = mThetaResult[3],
    filter.results = lFilterResults,
    reuse.diff.optim = lReuseNew
  ))
}

#' Hidden state sampler
#' 
#' Sampling the volatility, \eqn{h}. It calls the Gaussian simulation smoother
#' to sample the \eqn{\eta} error terms and then the volatility is constructed
#' using its autoregressive structure.
#' @author hdarjus \email{hdarjus@gmail.com}
#' @param dPhi newest \eqn{\phi} value
#' @param dSigma2 newest \eqn{\sigma^2} value
#' @param dRho newest \eqn{\rho} value
#' @param dMu newest \eqn{\mu} value
#' @param vD vector of signs
#' @param vS vector of current mixture states
#' @param lFilterResults results of \code{fnAugKalmanFilter}
#' @return Vector of volatility values.
#' @seealso \code{fnAugKalmanFilter, fnSimulationSmoother}
fnSampleH <- function (dPhi, dSigma2, dRho, dMu, vD, vS, lFilterResults, sCentering) {
  va <- dfModelConstants$a[vS]
  vm <- dfModelConstants$m[vS]
  
  lSimSmoothResult <- fnSimulationSmoother(dMu, lFilterResults, sCentering)
  vEth <- lSimSmoothResult$eta
  dEth0 <- lSimSmoothResult$eta0
  
  iPeriods <- length(lFilterResults$D)
  vH <- rep(0, iPeriods)
  if (sCentering == "centered") {
    vH[1] <- dMu + dEth0
    vDt <- dMu * (1-dPhi) + dRho * sqrt(dSigma2) * vD * va * exp(vm/2)
  } else if (sCentering == "non-centered") {
    vH[1] <- dEth0
    vDt <- dRho * vD * va * exp(vm/2)
  } else {
    stop("Invalid centering!")
  }
  
  for (i in 1:(iPeriods-1)) {
    vH[i+1] <- vDt[i] + dPhi * vH[i] + vEth[i]
  }
  
  return(vH)
}

#' Mixture state sampler
#' 
#' This implements the first step of Omori, the mixture states sampling.
#' @author hdarjus \email{hdarjus@gmail.com}
#' @param vYStar vector of transformed measurments
#' @param vH vector of hidden states
#' @param vD vector of \eqn{d_t}
#' @param dMu parameter \eqn{\mu}
#' @param dSigma2 parameter \eqn{\sigma^2}
#' @param dRho parameter \eqn{\rho}
#' @param dPhi parameter \eqn{\phi}
#' @return Vector of samples of the mixture state based on its posterior distribution.
#' @seealso \code{\link{fnStatePostDist}, \link{fnDrawDirichlet}}
fnSampleMixtureState <- function (vYStar, vH, vD, dMu, dSigma2, dRho, dPhi, sCentering) {
  if (sCentering == "centered") {
    vEpsStar <- vYStar - vH
    vEta <- (vH[-1] - dMu) - dPhi*(vH[-length(vH)] - dMu)
  } else if (sCentering == "non-centered") {
    vEpsStar <- vYStar - dMu - sqrt(dSigma2)*vH
    vEta <- vH[-1] - dPhi*vH[-length(vH)]
  } else {
    stop("Invalid centering!")
  }
  mPostDist <- fnMixtureStatePostDist(vEpsStar, vEta, vD, dMu, dSigma2, dRho, sCentering)
  mNewStates <- fnDrawDirichlet(1, mPostDist)
  return(as.numeric(mNewStates))
}

#' Sampler for \eqn{\mu}
#' 
#' This is the first part of step 2b, sampling \eqn{\mu} using
#' results from the augmented Kalman filter acquired in step 2a.
#' @author hdarjus \email{hdarjus@gmail.com}
#' @param dq output of augmented Kalman filter from step 2a
#' @param dQ output of augmented Kalman filter from step 2a
#' @return Column vector of newly sampled \eqn{\mu} values.
#' @seealso \code{fnSampleTheta, fnAugKalmanFilter, fnSimulSmoother}
fnSampleMu <- function (dq, dQ) {
  return(rnorm(1, dq/dQ, 1/sqrt(dQ)))
}

#' Correction for misspecification
#' 
#' Resamples the results coming from \code{fnMCMCSampler}. In principle this is needed
#' since \code{fnMCMCSampler} samples the approximate posterior density.
#' @param iN size of the new sample
#' @param lMCMCResult list coming from the \code{fnMCMCSampler function}
#' @return Data frame of corrected samples. Rows correspond to samples,
#'   columns \code{phi, rho, sigma2, mu} to the parameters.
#' @seealso \code{fnMCMCSampler}
#' @export
fnReweight <- function (iN, lMCMCResult) {
  vProbs <- lMCMCResult$weights
  vChosenIndices <- sample.int(n = length(vProbs), size = iN, replace = T, prob = vProbs)
  return(lMCMCResult$samples[vChosenIndices, ])
}
