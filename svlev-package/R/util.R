#' Utility functions

#' Augmented Kalman filter
#' 
#' Implementation of the augmented Kalman filter for the special
#' case of the current model. The augmented Kalman filter is a
#' modified version of the standard Kalman filter that
#' considers the the error terms' correlation and the constant
#' \eqn{\mu} in the state equation.
#' It is used to calculate the likelihood of \eqn{\theta=(\phi,\sigma,\rho)},
#' and to calculate the posterior of \eqn{\mu}.
#' For a specification, see Appendix B of Nakajima, 2009.
#' @author hdarjus \email{hdarjus@gmail.com}
#' @param dPhi parameter \eqn{\phi}
#' @param dSigma2 square of parameter \eqn{\sigma}
#' @param dRho parameter \eqn{\rho}
#' @param va numeric vector of mixture states constants
#' @param vb numeric vector of mixture states constants
#' @param vm numeric vector of mixture states constants
#' @param vv numeric vector of mixture states constants
#' @param vD vector of signs
#' @param vYStar vector of y*
#' @param dMuMu prior mean of \eqn{\mu}
#' @param dSigma2Mu prior variance of \eqn{\mu}
#' @return List with elements \code{sigma2, D, J1, L, f, F, hts, v, Q, q, jt22, h1var}, from which
#'   \code{D, J1, L, f, F, h1var, v} are numeric vectors of length T,
#'   and \code{sigma2, Q, q, jt22, h1var} are numbers. Except for \code{sigma2}, all of them are
#'   partial results helping later calculations.
#' @references Nakajima, Jouchi, and Omori, Yasuhiro.
#'   "Leverage, heavy-tails and correlated jumps in stochastic volatility models."
#'   Computational Statistics & Data Analysis 53.6 (2009): 2335-2353.
#' @export
fnAugKalmanFilter <- function (dPhi, dSigma2, dRho, va, vb, vm, vv, vD, vYStar, dMuMu, dSigma2Mu, sCentering) {
  lFilterResult <- switch(sCentering,
    centered = fnAugKalmanFilterCentered(dPhi, dSigma2, dRho, va, vb, vm, vv, vD, vYStar, dMuMu, dSigma2Mu),
    "non-centered" = fnAugKalmanFilterNonCentered(dPhi, dSigma2, dRho, va, vb, vm, vv, vD, vYStar, dMuMu, dSigma2Mu),
    stop("Invalid centering input"))
  return(lFilterResult)
}

fnAugKalmanFilterCentered <- function (dPhi, dSigma2, dRho, va, vb, vm, vv, vD, vYStar, dMuMu, dSigma2Mu) {
  iN <- length(vYStar)
  
  dSigma <- sqrt(dSigma2)
  vgammats <- vD * dRho * dSigma * exp(vm/2)
  vhts <- vb * vv * vgammats
  djt22 <- dSigma2 * (1 - dRho^2)  # square of second value of J
  dH1var <- dSigma2/(1-dPhi^2)  # unstable
  
  da <- 0
  dA <- -1
  dP <- dH1var
  
  # Initialize returned values
  vDret <- rep(0, iN)
  vJ <- rep(0, iN)  # only the first values of J
  vL <- rep(0, iN)
  vf <- rep(0, iN)
  vF <- rep(0, iN)
  dQ <- 1/dSigma2Mu
  dq <- dMuMu * dQ
  
  # Calculate everything inside a foor loop
  for (iI in seq_len(iN)) {
    vDret[iI] <- dP + vv[iI]^2
    dK <- (dPhi * dP + vhts[iI] * vv[iI])/vDret[iI]
    vL[iI] <- dPhi - dK
    vJ[iI] <- vhts[iI] - dK * vv[iI]
    dP <- dPhi * dP * vL[iI] + vhts[iI] * vJ[iI] + djt22
    vf[iI] <- vYStar[iI] - vm[iI] - da
    vF[iI] <- -dA
    da <- va[iI] * vgammats[iI] + dPhi * da + dK * vf[iI]
    dA <- dPhi - 1 + dPhi * dA + dK * vF[iI]
    
    dq <- dq + vF[iI] * vf[iI] / vDret[iI]
    dQ <- dQ + vF[iI]^2 / vDret[iI]
  }
  
  return(list(
    sigma2 = dSigma2,
    D = vDret,
    J1 = vJ,  # only the first values of J
    L = vL,
    f = vf,
    "F" = vF,
    hts = vhts,
    v = vv,
    Q = dQ,
    q = dq,
    jt22 = djt22,
    h1var = dH1var
  ))
}

fnAugKalmanFilterNonCentered <- function (dPhi, dSigma2, dRho, va, vb, vm, vv, vD, vYStar, dMuMu, dSigma2Mu) {
  iN <- length(vYStar)
  
  dSigma <- sqrt(dSigma2)
  vgammats <- vD * dRho * exp(vm/2)
  vhts <- vb * vv * vgammats
  djt22 <- (1 - dRho^2)  # square of second value of J
  dH1var <- 1/(1-dPhi^2)  # unstable
  
  da <- 0
  dA <- 0
  dP <- dH1var
  
  # Initialize returned values
  vDret <- rep(0, iN)
  vJ <- rep(0, iN)  # only the first values of J
  vL <- rep(0, iN)
  vf <- rep(0, iN)
  vF <- rep(0, iN)
  dQ <- 1/dSigma2Mu
  dq <- dMuMu * dQ
  
  # Calculate everything inside a foor loop
  for (iI in seq_len(iN)) {
    vDret[iI] <- dSigma2 * dP + vv[iI]^2
    dK <- (dSigma * dPhi * dP + vhts[iI] * vv[iI])/vDret[iI]
    vL[iI] <- dPhi - dK * dSigma
    vJ[iI] <- vhts[iI] - dK * vv[iI]
    dP <- dPhi * dP * vL[iI] + vhts[iI] * vJ[iI] + djt22
    vf[iI] <- vYStar[iI] - vm[iI] - da * dSigma
    vF[iI] <- 1 - dA * dSigma
    da <- va[iI] * vgammats[iI] + dPhi * da + dK * vf[iI]
    dA <- dPhi * dA + dK * vF[iI]
    
    dq <- dq + vF[iI] * vf[iI] / vDret[iI]
    dQ <- dQ + vF[iI]^2 / vDret[iI]
  }
  
  return(list(
    sigma2 = dSigma2,
    D = vDret,
    J1 = vJ,  # only the first values of J
    L = vL,
    f = vf,
    "F" = vF,
    hts = vhts,
    v = vv,
    Q = dQ,
    q = dq,
    jt22 = djt22,
    h1var = dH1var
  ))
}

#' Gaussian simulation smoother
#' 
#' The Gaussian simulation smoother is used to sample the hidden state in
#' Gaussian hidden state models. It samples the error term, \eqn{\eta}, we
#' can calculate the hidden state using \eqn{\eta}. Uses values acquired
#' in the \eqn{\theta} sampling step.
#' @author hdarjus \email{hdarjus@gmail.com}
#' @param dMu parameter \eqn{\mu}
#' @param lFilterResults results of \code{fnAugKalmanFilter}
#' @return List of 2 elements: \code{eth, eth0}, where \code{eth} is a vector of error terms
#'   1 to T, and \code{eth0} is the 0th \eqn{\eta}, a number.
#' @references De Jong, Piet, and Shephard, Neil.
#'   "The simulation smoother for time series models."
#'   Biometrika (1995): 339-350.
#' @export
fnSimulationSmoother <- function (dMu, lFilterResults, sCentering) {
  lSmoothingResult <- switch(sCentering,
    centered = fnSimulationSmootherCentered(dMu, lFilterResults),
    "non-centered" = fnSimulationSmootherNonCentered(dMu, lFilterResults),
    stop("Invalid centering input"))
  return(lSmoothingResult)
}

fnSimulationSmootherCentered <- function (dMu, lFilterResults) {
  vDret <- lFilterResults$D
  vJ <- lFilterResults$J1
  vL <- lFilterResults$L
  vhts <- lFilterResults$hts
  vv <- lFilterResults$v
  djt22 <- lFilterResults$jt22
  dH1var <- lFilterResults$h1var
  ve <- lFilterResults$f - lFilterResults[["F"]] * dMu
  
  # Results to return
  vEta <- rep(0, length(vDret))
  dEta0 <- NULL
  
  # Init values, then loop
  dr <- 0
  dU <- 0
  for (iI in rev(seq_along(vDret))) {
    dhjpj <- vhts[iI] * vJ[iI] + djt22
    dDinv <- 1/vDret[iI]
    
    dC <- (vhts[iI]^2 + djt22) - (vhts[iI] * vv[iI])^2 * dDinv - dU * (vhts[iI] * vJ[iI] + djt22)^2
    dEps <- rnorm(1, 0, sqrt(dC))
    dV <- vhts[iI] * vv[iI] * dDinv + dU * vL[iI] * dhjpj
    vEta[iI] <- vhts[iI] * vv[iI] * ve[iI] * dDinv + dhjpj * dr + dEps
    
    dr <- ve[iI] * dDinv + vL[iI] * dr - dV * dEps / dC
    dU <- dDinv + dU * vL[iI]^2 + dV^2 / dC
  }
  
  # Case i = 0
  dC <- dH1var * (1 - dU * dH1var)
  dEps <- rnorm(1, 0, sqrt(dC))
  dEta0 <- dH1var * dr + dEps
  
  return(list(eta = vEta, eta0 = dEta0))
}

fnSimulationSmootherNonCentered <- function (dMu, lFilterResults) {
  dSigma2 <- lFilterResults$sigma2
  vDret <- lFilterResults$D
  vJ <- lFilterResults$J1
  vL <- lFilterResults$L
  vhts <- lFilterResults$hts
  vv <- lFilterResults$v
  djt22 <- lFilterResults$jt22
  dH1var <- lFilterResults$h1var
  ve <- lFilterResults$f - lFilterResults[["F"]] * dMu
  
  dSigma <- sqrt(dSigma2)
  
  # Results to return
  vEta <- rep(0, length(vDret))
  dEta0 <- NULL
  
  # Init values, then loop
  dr <- 0
  dU <- 0
  for (iI in rev(seq_along(vDret))) {
    dhjpj <- vhts[iI] * vJ[iI] + djt22
    dDinv <- 1/vDret[iI]
    
    dC <- (vhts[iI]^2 + djt22) - (vhts[iI] * vv[iI])^2 * dDinv - dU * (vhts[iI] * vJ[iI] + djt22)^2
    dEps <- rnorm(1, 0, sqrt(dC))
    dV <- dSigma * vhts[iI] * vv[iI] * dDinv + dU * vL[iI] * dhjpj
    vEta[iI] <- vhts[iI] * vv[iI] * ve[iI] * dDinv + dhjpj * dr + dEps
    
    dr <- dSigma * ve[iI] * dDinv + vL[iI] * dr - dV * dEps / dC
    dU <- dSigma2 * dDinv + dU * vL[iI]^2 + dV^2 / dC
  }
  
  # Case i = 0
  dC <- dH1var * (1 - dU * dH1var)
  dEps <- rnorm(1, 0, sqrt(dC))
  dEta0 <- dH1var * dr + dEps
  
  return(list(eta = vEta, eta0 = dEta0))
}

#' Calculates the reweighting weights
#' 
#' Calculates weights for the correction of the posterior distribution.
#' @author hdarjus \email{hdarjus@gmail.com}
#' @param vYStar vector of transformed measurements
#' @param vD vector of signs
#' @param mSampleH matrix of all samples from the posterior of \eqn{h}; rows correspond to draws, columns to timepoints
#' @param dfSamples all samples from the posterior of \eqn{\mu,\sigma,\rho,\phi}; rows correspond to draws, columns \code{phi, sigma2, rho, mu} to parameters
#' @return Vector of weights.
fnCalcWeights <- function (vYStar, vD, mSampleH, dfSamples, sCentering) {
  cSamp <- nrow(mSampleH)
  cTime <- ncol(mSampleH)
  mYStar <- outer(rep(1, cSamp), vYStar)
  if (sCentering == "centered") {
    mMu <- outer(dfSamples$mu, rep(1, cTime-1))
    mPhi <- outer(dfSamples$phi, rep(1, cTime-1))
    mEpsStar <- mYStar - mSampleH
    mEta <- (mSampleH[, -1] - mMu) - mPhi * (mSampleH[, -cTime] - mMu)
    mEta <- cbind(mEta, 0)
  } else if (sCentering == "non-centered") {
    mMu <- outer(dfSamples$mu, rep(1, cTime))
    mSigma <- outer(sqrt(dfSamples$sigma2), rep(1, cTime))
    mPhi <- outer(dfSamples$phi, rep(1, cTime-1))
    mEpsStar <- mYStar - mMu - mSigma * mSampleH
    mEta <- mSampleH[, -1] - mPhi * mSampleH[, -cTime]
    mEta <- cbind(mEta, 0)
  } else {
    stop("Invalid centering method")
  }
  vWeights <- rep(0, cSamp)
  for (iSamp in seq_len(cSamp)) {
    vWeights[iSamp] <- prod(fnErrorTrueDensity(mEpsStar[iSamp, ], mEta[iSamp, ], vD, sqrt(dfSamples[iSamp, "sigma2"]), dfSamples[iSamp, "rho"], sCentering) /
                              fnErrorApproxDensity(mEpsStar[iSamp, ], mEta[iSamp, ], vD, sqrt(dfSamples[iSamp, "sigma2"]), dfSamples[iSamp, "rho"], sCentering))
  }
  vWeights <- vWeights/sum(vWeights)
  return(vWeights)
}
