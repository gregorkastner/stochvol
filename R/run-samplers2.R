#' Main function
#' 
#' @export
svlsample <- function (n=10000, y, thin=1, burnin=1000, init=NULL, prior=NULL,
                        stdev=NA,
                        seed=NULL, strategy=c("centered", "non-centered")) {
  
  cStrategy <- match.arg(strategy, c("centered", "non-centered"), several.ok = TRUE)
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
  
  dPhi <- dPhi0; dRho <- dRho0; dSigma2 <- dSigma20; dMu <- dMu0
  lLatent <- list(h = vH0, ht = (vH0-dMu0)/sqrt(dSigma20), s = NULL)
  vH <- vH0
  vHt <- (vH0-dMu0)/sqrt(dSigma20)
  oRuntime <- system.time({
    post_sample <- svlsample_cpp(n, vY, vYStar, vD, thin, burnin, dPhi0, dRho0, dSigma20, dMu0,
                                 vH0, cPriorPhi[1], cPriorPhi[2], cPriorRho[1], cPriorRho[2],
                                 cPriorSigma2[1], cPriorSigma2[2], cPriorMu[1], cPriorMu[2],
                                 stdev, cStrategy, dfModelConstants)
  }, gcFirst = FALSE)
  params <- post_sample[, 1:4]
  colnames(params) <- c("phi", "rho", "sigma2", "mu")
  invisible(list(param = params, vol = post_sample[, -4:-1], runtime = oRuntime))
}

