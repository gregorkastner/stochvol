#' Main function
#' 
#' @export
svlsample <- function (y, draws = 10000, burnin = 1000, designmatrix = NA,
                       priormu = c(0, 100), priorphi = c(5, 1.5), priorsigma = 1,
                       priorrho = c(3, 5), priorbeta = c(0, 10000),
                       thinpara = 1, thinlatent = 1, thintime = 1,
                       quiet = FALSE, startpara, startlatent, expert, ...) {
  # TODO quiet designmatrix priorbeta gammaprior(expert)
  # Some error checking for y
  if (inherits(y, "svsim")) {
    y <- y[["y"]]
    warning("Extracted data vector from 'svsim'-object.")
  }
  if (!is.numeric(y)) stop("Argument 'y' (data vector) must be numeric.")
  if (length(y) < 2) stop("Argument 'y' (data vector) must contain at least two elements.")

  # Some error checking for draws
  if (!is.numeric(draws) || length(draws) != 1 || draws < 1) {
    stop("Argument 'draws' (number of MCMC iterations after burn-in) must be a single number >= 1.")
  } else {
    draws <- as.integer(draws)
  }

  # Some error checking for burnin
  if (!is.numeric(burnin) || length(burnin) != 1 || burnin < 0) {
    stop("Argument 'burnin' (burn-in period) must be a single number >= 0.")
  } else {
    burnin <- as.integer(burnin)
  }

  # Some error checking for designmatrix
  if (any(is.na(designmatrix))) {
    designmatrix <- matrix(NA)
  } else {
    if (any(grep("ar[0-9]+$", as.character(designmatrix)[1]))) {
      order <- as.integer(gsub("ar", "", as.character(designmatrix)))
      if (length(y) <= order + 1) stop("Time series 'y' is to short for this AR process.")
      designmatrix <- matrix(rep(1, length(y) - order), ncol = 1)
      colnames(designmatrix) <- c("const")
      if (order >= 1) {
        for (i in 1:order) {
          oldnames <- colnames(designmatrix)
          designmatrix <- cbind(designmatrix, y[(order-i+1):(length(y)-i)])
          colnames(designmatrix) <- c(oldnames, paste0("ar", i))
        }
        y <- y[-(1:order)]
      }
    }
    if (!is.numeric(designmatrix)) stop("Argument 'designmatrix' must be a numeric matrix or an AR-specification.")
    if (!is.matrix(designmatrix)) {
      designmatrix <- matrix(designmatrix, ncol = 1)
    }
    if (!nrow(designmatrix) == length(y)) stop("Number of columns of argument 'designmatrix' must be equal to length(y).")
  }

  # Some error checking for the prior parameters 
  if (!is.numeric(priormu) || length(priormu) != 2) {
    stop("Argument 'priormu' (mean and variance for the Gaussian prior for mu) must be numeric and of length 2.")
  }

  if (!is.numeric(priorphi) || length(priorphi) != 2) {
    stop("Argument 'priorphi' (shape1 and shape2 parameters for the Beta prior for (phi+1)/2) must be numeric and of length 2.")
  }

  if (!is.numeric(priorsigma) || length(priorsigma) != 1 || priorsigma <= 0) {
    stop("Argument 'priorsigma' (scaling of the chi-squared(df = 1) prior for sigma^2) must be a single number > 0.")
  }

  if (!is.numeric(priorrho) || length(priorrho) != 2) {
    stop("Argument 'priorrho' (shape1 and shape2 parameters for the Beta prior for (rho+1)/2) must be numeric and of length 2.")
  }

  if (!is.numeric(priorbeta) || length(priorbeta) != 2) {
    stop("Argument 'priorbeta' (means and sds for the independent Gaussian priors for beta) must be numeric and of length 2.")
  }

  # Some error checking for thinpara
  if (!is.numeric(thinpara) || length(thinpara) != 1 || thinpara < 1) {
    stop("Argument 'thinpara' (thinning parameter for mu, phi, sigma, and rho) must be a single number >= 1.")
  } else {
    thinpara <- as.integer(thinpara)
  }

  # Some error checking for thinlatent
  if (!is.numeric(thinlatent) || length(thinlatent) != 1 || thinlatent < 1) {
    stop("Argument 'thinlatent' (thinning parameter for the latent log-volatilities) must be a single number >= 1.")
  } else {
    thinlatent <- as.integer(thinlatent)
  }

  # Some error checking for thintime
  if (!is.numeric(thintime) || length(thintime) != 1 || thintime < 1) {
    if (thintime == 'firstlast') {
      thintime <- length(y) - 1L
    } else {
      stop("Argument 'thintime' (thinning parameter for time) must be a single number >= 1 or 'firstlast'.")
    }
  } else {
    thintime <- as.integer(thintime)
  }

  # Some input checking for startpara
  startparadefault <- list(mu = priormu[1],
                           phi = 2 * (priorphi[1] / sum(priorphi)) - 1,
                           sigma = priorsigma,
                           rho = 2 * (priorrho[1] / sum(priorrho)) - 1)
  if (missing(startpara)) {
    startpara <- startparadefault
  } else {
    if (!is.list(startpara))
      stop("Argument 'startpara' must be a list. Its elements must be named 'mu', 'phi', 'sigma', and 'rho'.")

    if (!is.numeric(startpara[["mu"]]))
      stop('Argument \'startpara[["mu"]]\' must exist and be numeric.')

    if (!is.numeric(startpara[["phi"]]))
      stop('Argument \'startpara[["phi"]]\' must exist and be numeric.')

    if (abs(startpara[["phi"]]) >= 1)
      stop('Argument \'startpara[["phi"]]\' must be between -1 and 1.')

    if (!is.numeric(startpara[["sigma"]]))
      stop('Argument \'startpara[["sigma"]]\' must exist and be numeric.')

    if (startpara[["sigma"]] <= 0)
      stop('Argument \'startpara[["sigma"]]\' must be positive.')

    if (!is.numeric(startpara[["rho"]]))
      stop('Argument \'startpara[["rho"]]\' must exist and be numeric.')

    if (abs(startpara[["rho"]]) >= 1)
      stop('Argument \'startpara[["rho"]]\' must be between -1 and 1.')
  }

  # Some input checking for startlatent
  if (missing(startlatent)) {
    startlatent <- rep(-10, length(y))
  } else {
    if (!is.numeric(startlatent) | length(startlatent) != length(y))
      stop("Argument 'startlatent' must be numeric and of the same length as the data 'y'.")
  }

  if (!quiet) {  # TODO
    cat(paste("\nCalling ", parameterization, " MCMC sampler with ", draws+burnin, " iter. Series length is ", length(y), ".\n",sep=""), file=stderr())
    flush.console()
  }

  if (.Platform$OS.type != "unix") myquiet <- TRUE else myquiet <- quiet  # Hack to prevent console flushing problems with Windows

  # Some error checking for expert
  strategies <- c("centered", "non-centered")
  expertdefault <- list(parameterization = rep(strategies, 5),  # default: ASISx5
                        mhcontrol = 0.1,
                        gammaprior = TRUE)  # TODO
  if (missing(expert)) {
    parameterization <- expertdefault$parameterization
    mhcontrol <- expertdefault$mhcontrol
    gammaprior <- expertdefault$gammaprior
  } else {
    expertnames <- names(expert)
    if (!is.list(expert) || is.null(expertnames) || any(expertnames == ""))
      stop("Argument 'expert' must be a named list with nonempty names.")
    if (length(unique(expertnames)) != length(expertnames))
      stop("No duplicate elements allowed in argument 'expert'.")
    allowednames <- c("parameterization", "mhcontrol", "gammaprior")
    exist <- pmatch(expertnames, allowednames)
    if (any(is.na(exist)))
      stop(paste("Illegal element '", paste(expertnames[is.na(exist)], collapse="' and '"), "' in argument 'expert'.", sep=''))

    expertenv <- list2env(expert) 

    if (exists("parameterization", expertenv)) {
      if (!is.character(expert[["parameterization"]]))
        stop("Argument 'parameterization' must be either a vector of 'centered', 'non-centered' values or a character string of form 'asis#' with # a positive integer.")
      nmatches <- grep("^asis[1-9][0-9]*$", expert[["parameterization"]])
      if (length(nmatches) == 0) {
        parameterization <- match.arg(expert[["parameterization"]], strategies, several.ok = TRUE)
      } else if (length(nmatches) == 1) {
        parameterization <- rep(strategies, nmatches)
      } else {
        parameterization <- NA
      }
      if (!all(parameterization %in% strategies)) {
        stop("Argument 'parameterization' must be either a vector of 'centered', 'non-centered' values or a character string of form 'asis#' with # a positive integer.")
      }
    } else {
      parameterization <- expertdefault$parameterization
    }

    # Remark: mhcontrol > 0 controls stepsize of log-random-walk proposal
    if (exists("mhcontrol", expertenv)) {
      mhcontrol <- expert[["mhcontrol"]]
      if (!is.numeric(mhcontrol) || length(mhcontrol) != 1 || mhcontrol <= 0)
        stop("Argument 'mhcontrol' must be a single positive number.")
    } else {
      mhcontrol <- expertdefault$mhcontrol
    }

    # use a Gamma prior for sigma^2 in C?
    if (exists("gammaprior", expertenv)) {
      gammaprior <- expert[["gammaprior"]]
      if (!is.logical(gammaprior)) stop("Argument 'gammaprior' must be TRUE or FALSE.")
    } else {
      gammaprior <- expertdefault$gammaprior
    }
  }
  
  myoffset <- if (any(y^2 == 0)) sd(y)/10000 else 0
  ystar <- log(y^2+myoffset)
  if (dOffest > 0) {
    warning(paste("Argument 'y' (data vector) contains zeros. I am adding an offset constant of size ", myoffset, " to do the auxiliary mixture sampling. If you want to avoid this, you might consider de-meaning the returns before calling this function.", sep=""))
  }
  d <- 2*(y > 0) - 1
  h <- startlatent
  
  phi <- startpara$phi; rho <- startpara$rho; sigma2 <- startpara$sigma^2; mu <- startpara$mu

  runtime <- system.time({
    res <- svlsample_cpp(draws, y, ystar, d, burnin, thinpara, thinlatent, thintime,
                                 phi, rho, sigma2, mu, h,
                                 priorphi[1], priorphi[2], priorrho[1], priorrho[2],
                                 priorsigma[1], priorsigma[2], priormu[1], priormu[2],
                                 mhcontrol, parameterization, gammaprior, dfModelConstants)
  })
  
  res$para <- res$para[, c(4,1,3,2)]
  res$para[, 3] <- sqrt(res$para[, 3])
  colnames(res$para) <- c("mu", "phi", "sigma", "rho")
  colnames(res$latent) <- paste0('h_', seq(1, length(y), by=thintime))
  # create svldraws class
  res$runtime <- runtime
  res$y <- y
  res$para <- coda::mcmc(res$para, burnin+thin, burnin+draws, thinpara)
  res$latent <- coda::mcmc(res$latent, burnin+thin, burnin+draws, thinlatent)
  res$thinning <- list(para = thinpara, latent = thinlatent, time = thintime)
  res$priors <- list(mu = priormu, phi = priorphi, sigma = priorsigma, rho = priorrho)
  if (!any(is.na(designmatrix))) {
    res$beta <- mcmc(res$beta[seq(burnin+thinpara+1, burnin+draws+1, thinpara),,drop=FALSE], burnin+thinpara, burnin+draws, thinpara)  # TODO
    colnames(res$beta) <- paste("b", 0:(NCOL(designmatrix)-1), sep = "_")
    res$priors <- c(res$priors, "beta" = list(priorbeta), "designmatrix" = list(designmatrix))
  }
  class(res) <- c("svldraws", "svdraws")

  if (!quiet) {
    cat("Done!\n", file=stderr())
    cat("Summarizing posterior draws... ", file=stderr())
  }
  res <- updatesummary(res, ...)

  if (!quiet) cat("Done!\n\n", file=stderr())
  res
}

