#' @export
fnGenLogSV <- function (iN, dRho, dSigma, dPhi, dMu, iSeed=NA) {
  if (is.numeric(iSeed))
    set.seed(iSeed)
  # generate residuals, h and y
  vEps <- rnorm(iN)
  vEta <- dSigma * (dRho * vEps + sqrt(1-dRho^2) * rnorm(iN))
  vH <- rep(0, n=iN)
  vH[1] <- rnorm(1, dMu, sd = dSigma/sqrt(1-dPhi^2))
  for (iI in 1:(iN-1)) {
    vH[iI+1] <- dMu + dPhi*(vH[iI] - dMu) + vEta[iI]
  }
  vY <- vEps * exp(vH/2)
  list(params=c(phi=dPhi, mu=dMu, sigma2=dSigma2, rho=dRho),
       y=vY, h=vH, eps=vEps, eta=vEta)
}

# Generate centered and non-centered way
gen.rho.cnc <- function () {
  iT <- 10000
  Z <- matrix(rnorm(2*iT), iT, 2)
  mTrans <- matrix(c(dRho, sqrt(1-dRho^2), -sqrt(1-dRho^2), dRho), 2, 2, byrow = TRUE)
  mTransChange <- diag(c(1, -1))
  V1 <- tcrossprod(Z, mTrans)
  
  vY <- matrix(NA, iT, 2)
  vH <- rep(NA, iT)
  vH[1] <- dMu
  for (it in seq_len(iT)) {
    vY[it, 1] <- exp(vH[it]/2) * Z[it, 1]  # Z[it, 1] == cbind(dRho, -sqrt(1-dRho^2)) %*% V1[it, ]
    vY[it, 2] <- exp(vH[it]/2) * cbind(dRho, sqrt(1-dRho^2)) %*% V1[it, ]
    if (it < iT) {
      vH[it+1] <- dMu + dPhi*(vH[it]-dMu) + sqrt(dSigma2)*V1[it, 1]
    }
  }
  vHt <- (vH-dMu)/sqrt(dSigma2)
  cbind(y=vY, h=vH, ht=vHt)
}
