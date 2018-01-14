#' Posterior distributions/likelihoods

#' Mixture state posteriors
#' 
#' Calculates the posterior distribution of the mixture components.
#' @author hdarjus \email{hdarjus@gmail.com}
#' @param vEpsStar vector of \eqn{\epsilon*_t}
#' @param vEta vector of \eqn{\eta_t}
#' @param vD vector of \eqn{d_t}
#' @param dMu parameter \eqn{\mu}
#' @param dSigma2 parameter \eqn{\sigma^2}
#' @param dRho parameter \eqn{\rho}
#' @return Matrix of posterior probabilities of the mixture component states. Each row
#'   is a time period, each column contains the probabilities of a mixture component state the periods.
#' @details The derivations of the posterior distribution can be found under
#'   \url{https://www.sharelatex.com/project/58a89475bf3be6ae4fe5e6a0}.
fnMixtureStatePostDist <- function (vEpsStar, vEta, vD, dMu, dSigma2, dRho, sCentering) {
  # Parameter validation
  if (length(vEpsStar) != length(vD))
    stop("vEpsStar and vD have to be of same length!")
  if (length(vEta)+1 != length(vEpsStar))
    stop("vEta should have one less element than vEpsStar!")
  if (length(c(dMu, dSigma2, dRho)) != 3)
    stop("dMu, dSigma2 and dRho should be single numbers!")
  if (dSigma2 <= 0) 
    stop("dSigma2 has to be positive!")
  
  iN <- length(vEpsStar)
  iMixCount <- length(dfModelConstants$v)
  mPrior <- rep(1, iN) %o% dfModelConstants$p  # state prior probabilities
  mEpsStarLik <- outer(vEpsStar,  # likelihood of epsilon*
                       seq_len(iMixCount),
                       FUN = function (dEpsStar, iJ, vM, vV, vV2)
                         exp((-.5)/vV2[iJ] * (dEpsStar-vM[iJ])^2) / (sqrt(2*pi)*vV[iJ]),
                       dfModelConstants$m,
                       dfModelConstants$v,
                       dfModelConstants$v2)
  dSigma2Used <- if (sCentering == "centered") dSigma2 else 1
  mEtaLik <- outer(seq_len(iN-1),  # likelihood of eta
                   seq_len(iMixCount),
                   FUN = function (iT, iJ, vEta, vEpsStar, vD, vA, vB, vM)
                     exp((-.5)/(dSigma2Used*(1-dRho^2)) * (vEta[iT] - dRho*sqrt(dSigma2Used)*vD[iT]*exp(vM[iJ]/2)*(vA[iJ]+vB[iJ]*(vEpsStar[iT]-vM[iJ])))^2) / (sqrt(2*pi*(1-dRho^2)*dSigma2Used)),
                   vEta,
                   vEpsStar[-iN],  # to make it clear
                   vD[-iN],  # to make it clear
                   dfModelConstants$a,
                   dfModelConstants$b,
                   dfModelConstants$m)
  mEtaLik <- rbind(mEtaLik, 1)
  mResult <- mPrior * mEpsStarLik * mEtaLik
  vRowSums <- rowSums(mResult)
  vNiceRows <- (vRowSums > 0) & !is.na(vRowSums) & is.finite(vRowSums)
  mResult[vNiceRows, ] <- mResult[vNiceRows, ]/vRowSums[vNiceRows]  # using that matrices are represented by column
  if (!all(vNiceRows)) {
    mResult[!vNiceRows, ] <- 1/10
    warning("Mixture state posterior contained ugly rows!")
  }
  return(mResult)
}

#' Log posterior of theta/gamma
#' 
#' Calculates the log posterior distribution of \eqn{\gamma} = \eqn{(\phi, \sigma, \rho)},
#' used in the maximization of step 2a. It is called \eqn{\theta} before the maximization.
#' Constants are ignored.
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
#' @param vMuPriorParams vector of length 2, mean and variance of the normal prior of \eqn{\mu}
#' @param vPhiPriorParams vector of length 2, two parameters of the beta prior of \eqn{\phi}
#' @param vSigma2PriorParams vector of length 2, shape and rate of the gamma prior of inverse square \eqn{\sigma}
#' @param vRhoPriorParams vector of length 2, two parameters of the beta prior of \eqn{\rho}
#' @return List with 2 elements: \code{result} is the log of posterior of the joint variable \eqn{(\phi, \sigma, \rho)}
#'   evaluated at \code{(dPhi, dSigma2, dRho)}, conditional on the remaining arguments; \code{filter.results} is
#'   the result of the augmented Kalman filter. The mixture constants
#'   defined in \code{\link{dfModelConstants}} are also used with the mixture states.
#' @seealso \code{\link{fnGammaLogPostDistTransInput}, \link{dfModelConstants}, \link{fnAugKalmanFilter}}
fnThetaLogPostDist <- function (dPhi, dSigma2, dRho, va, vb, vm, vv, vD, vYStar, vMuPriorParams,
                                vPhiPriorParams, vSigma2PriorParams, vRhoPriorParams, sSigma2.Prior, sCentering) {
  lFilterResults <- fnAugKalmanFilterCpp(dPhi, dSigma2, dRho, va, vb, vm, vv, vD, vYStar, vMuPriorParams[1], vMuPriorParams[2], sCentering)
  # Likelihood
  dLik <- sum(log(abs(lFilterResults$D)))
  dLik <- dLik + log(abs(lFilterResults$Q))
  dLik <- dLik + sum(lFilterResults$f^2 / lFilterResults$D)
  dLik <- dLik - lFilterResults$q^2 / lFilterResults$Q
  dLik <- -0.5 * dLik
  # Prior
  dPri <- log(0.5) + dbeta((dPhi+1)/2, vPhiPriorParams[1], vPhiPriorParams[2], log = T)
  dPri <- dPri + switch(sSigma2.Prior,
                        inv.gamma = -2*log(dSigma2) + dgamma(1/dSigma2,
                                                             shape = vSigma2PriorParams[1],
                                                             rate = vSigma2PriorParams[2],
                                                             log = T),
                        gamma = dgamma(dSigma2,
                                       shape = vSigma2PriorParams[1],
                                       rate = vSigma2PriorParams[2],
                                       log = T))
  dPri <- dPri + log(0.5) + dbeta((dRho+1)/2, vRhoPriorParams[1], vRhoPriorParams[2], log = T)
  
  dRes <- dLik + dPri
  
  return(list(result = dRes, filter.results = lFilterResults))
}

#' Log posterior of theta/gamma with transformed inputs
#' 
#' Calculates the log posterior distribution of \eqn{\gamma} = \eqn{(\phi, \sigma, \rho)},
#' used in the maximization of step 2a. It is called \eqn{\theta} before the maximization.
#' This function does the same as \code{fnGammaLogPostDist}, just uses transformed inputs.
#' @author hdarjus \email{hdarjus@gmail.com}
#' @param dTransPhi \eqn{\log(1+\phi)/\log(1-\phi)}
#' @param dTransSigma2 \eqn{\log(\sigma^2)}
#' @param dTransRho \eqn{\log(1+\rho)/\log(1-\rho)}
#' @param va numeric vector of mixture states constants
#' @param vb numeric vector of mixture states constants
#' @param vm numeric vector of mixture states constants
#' @param vv numeric vector of mixture states constants
#' @param vD vector of signs
#' @param vYStar vector of y*
#' @param vMuPriorParams vector of length 2, mean and variance of the normal prior of \eqn{\mu}
#' @param vPhiPriorParams vector of length 2, two parameters of the beta prior of \eqn{\phi}
#' @param vSigmaPriorParams vector of length 2, shape and rate of the gamma prior of inverse square \eqn{\sigma}
#' @param vRhoPriorParams vector of length 2, two parameters of the beta prior of \eqn{\rho}
#' @return List with 2 elements: \code{result} is the log of posterior of the joint variable \eqn{(\phi, \sigma, \rho)}
#'   evaluated at \code{(2/(1+exp(-dTransPhi))-1, exp(dTransSigma2), 2/(1+exp(-dTransRho))-1)},
#'   conditional on the remaining arguments; \code{filter.results} is
#'   the result of the augmented Kalman filter. The mixture constants
#'   defined in \code{\link{dfModelConstants}} are also used with the mixture states.
#' @seealso \code{\link{fnGammaLogPostDist}, \link{dfModelConstants}}
fnThetaLogPostDistTransInput <- function (dTransPhi, dTransSigma2, dTransRho, va, vb, vm, vv, vD, vYStar, vMuPriorParams,
                                          vPhiPriorParams, vSigmaPriorParams, vRhoPriorParams, sSigma2.Prior, sCentering) {
  dPhi <- 2/(1+exp(-dTransPhi))-1
  dSigma2 <- min(c(exp(dTransSigma2), 1e200))
  dRho <- 2/(1+exp(-dTransRho))-1
  if (exp(dTransSigma2) > 1e200) {
    warning("dTransSigma2 was too large!")
  }
  dRes <- fnThetaLogPostDist(dPhi, dSigma2, dRho, va, vb, vm, vv, vD, vYStar, vMuPriorParams,
                             vPhiPriorParams, vSigmaPriorParams, vRhoPriorParams, sSigma2.Prior, sCentering)
  return(dRes)
}
