# Functions that assert some properties of the inputs

sv_stop <- function (what, name, should, and_not) {
  stop(name, " = ", prettify(what), " should ", should, "; not ", and_not)
}

prettify <- function (x) {
  if (is.character(x)) {
    x <- paste0("\"", x, "\"")
  }
  out <- paste0("c(", paste(x, collapse = ", "), ")")
  if (nchar(out) > 40L) {
    out <- paste(substr(out, start = 1, stop = 40L), "...)")
  }
  out
}

assert_numeric <- function (x, name) {
  if (!(isTRUE(is.numeric(x)) &&
        isTRUE(all(!is.na(x))))) {
    sv_stop(what = x, name = name, should = "be numeric", and_not = typeof(x))
  }
  TRUE
}

assert_logical <- function (x, name) {
  if (!(isTRUE(is.logical(x)) &&
        isTRUE(all(!is.na(x))))) {
    sv_stop(what = x, name = name, should = "be logical", and_not = typeof(x))
  }
  TRUE
}

assert_single <- function (x, name) {
  assert_length(x, 1L, name)
}

assert_length <- function (x, len, name) {
  if (!isTRUE(length(x) == len)) {
    sv_stop(what = x, name = name, should = paste("have length", len), and_not = length(x))
  }
  TRUE
}

assert_infinite <- function (x, name) {
  assert_numeric(x, name)

  if (!isTRUE(all(is.infinite(x)))) {
    sv_stop(what = x, name = name, should = "be infinite", and_not = x)
  }
  TRUE
}

assert_ge <- function (x, value, name) {
  assert_numeric(x, name)
  assert_numeric(value, "right hand side")

  if (!isTRUE(all(x >= value))) {
    sv_stop(what = x, name = name, should = paste("be greater-equal than", prettify(value)), and_not = "smaller")
  }
  TRUE
}

assert_gt <- function (x, value, name) {
  assert_numeric(x, name)
  assert_numeric(value, "right hand side")

  if (!isTRUE(all(x > value))) {
    sv_stop(what = x, name = name, should = paste("be greater than", prettify(value)), and_not = "smaller or equal")
  }
  TRUE
}

assert_lt <- function (x, value, name) {
  assert_numeric(x, name)
  assert_numeric(value, "right hand side")

  if (!isTRUE(all(x < value))) {
    sv_stop(what = x, name = name, should = paste("be less than", prettify(value)), and_not = "greater or equal")
  }
  TRUE
}

assert_positive <- function (x, name) {
  assert_gt(x, 0, name)
}

assert_nonnegative <- function (x, name) {
  assert_ge(x, 0, name)
}

assert_element <- function (x, v, name_x, name_v) {
  if (!((isTRUE(typeof(x) == typeof(v)) ||
         isTRUE(is.numeric(x) && is.numeric(v))) &&
        isTRUE(all(x %in% v)))) {
    sv_stop(what = x, name = name_x, should = paste("be a subset of", name_v, "=", prettify(v)), and_not = "different")
  }
  TRUE
}

validate_sv_priors <- function (priormu, priorphi, priorsigma, priornu, priorrho, priorbeta, priorlatent0) {
  name_priormu <- "Argument 'priormu' (mean and sd for the Gaussian prior for mu)"
  assert_numeric(priormu, name_priormu)
  assert_length(priormu, 2, name_priormu)
  assert_positive(priormu[2], "The second element in argument 'priormu' (sd for the Gaussian prior for mu)")

  name_priorphi <- "Argument 'priorphi' (shape1 and shape2 parameters for the Beta prior for (phi + 1) / 2)"
  assert_positive(priorphi, name_priorphi)
  assert_length(priorphi, 2, name_priorphi)

  name_priorsigma2 <- "Argument 'priorsigma' (scaling of the chi-squared(df = 1) prior for sigma^2)"
  assert_positive(priorsigma, name_priorsigma2)
  assert_single(priorsigma, name_priorsigma2)

  name_priornu <- "Argument 'priornu' (rate parameter for the exponential prior for the df)"
  assert_nonnegative(priornu, name_priornu)
  assert_single(priornu, name_priornu)

  name_priorrho <- "Argument 'priorrho' (shape1 and shape2 parameters for the Beta prior for (rho + 1) / 2)"
  if (isTRUE(is.numeric(priorrho))) {
    assert_positive(priorrho, name_priorrho)
    assert_length(priorrho, 2, name_priorrho)
  } else if (!isTRUE(is.na(priorrho))) {
    stop(name_priorrho, " should be either NA or a positive vector of length 2")
  }

  name_priorbeta <- "Argument 'priorbeta' (means and sds for the independent Gaussian priors for beta)"
  assert_numeric(priorbeta, name_priorbeta)
  assert_length(priorbeta, 2, name_priorbeta)
  assert_positive(priorbeta[2], "The second element in argument 'priorbeta' sds for the independent Gaussian priors for beta)")

  name_priorlatent0 <- "Argument 'priorlatent0'"
  if (isTRUE(is.numeric(priorlatent0))) {
    assert_positive(priorlatent0, name_priorlatent0)
    assert_single(priorlatent0, name_priorlatent0)
  } else if (!isTRUE(priorlatent0 == "stationary")) {
    stop(name_priorlatent0, " should be either the character string \"stationary\" or a single positive number.")
  }
}

validate_thinning <- function (thinpara, thinlatent, keeptime) {
  name_thinpara <- "Argument 'thinpara' (thinning parameter for mu, phi, sigma, nu, and rho)"
  assert_single(thinpara, name_thinpara)
  assert_ge(thinpara, 1L, name_thinpara)

  name_thinlatent <- "Argument 'thinlatent' (thinning parameter for the latent log-volatilities)"
  assert_single(thinlatent, name_thinlatent)
  assert_ge(thinlatent, 1L, name_thinlatent)

  name_keeptime <- "Argument 'keeptime'"
  assert_single(keeptime, name_keeptime)
  assert_element(keeptime, c("all", "last"))
}

validate_initial_values <- function (startpara, startlatent, y, x) {
  assert_element(names(startpara), c("mu", "phi", "sigma", "nu", "rho", "beta", "latent0"),
                 "Elements of argument 'startpara'", "parameters of the model")

  assert_single(startpara$mu, "Provided initial value for parameter mu in 'startpara'")
  assert_numeric(startpara$mu, "Provided initial value for parameter mu in 'startpara'")

  assert_single(startpara$phi, "Provided initial value for parameter phi in 'startpara'")
  assert_gt(startpara$phi, -1, "Provided initial value for parameter phi in 'startpara'")
  assert_lt(startpara$phi, 1, "Provided initial value for parameter phi in 'startpara'")

  assert_single(startpara$sigma, "Provided initial value for parameter sigma in 'startpara'")
  assert_positive(startpara$sigma, "Provided initial value for parameter sigma in 'startpara'")

  assert_single(startpara$nu, "Provided initial value for parameter nu in 'startpara'")
  assert_gt(startpara$nu, 2, "Provided initial value for parameter nu in 'startpara'")

  assert_single(startpara$rho, "Provided initial value for parameter rho in 'startpara'")
  assert_gt(startpara$rho, -1, "Provided initial value for parameter rho in 'startpara'")
  assert_lt(startpara$rho, 1, "Provided initial value for parameter rho in 'startpara'")

  assert_length(startpara$beta, NCOL(x), "Provided initial values for the regressors beta in 'startpara'")
  assert_numeric(startpara$beta, "Provided initial values for the regressors beta in 'startpara'")

  assert_single(startpara$latent0, "Provided initial value for parameter latent0 in 'startpara'")
  assert_numeric(startpara$latent0, "Provided initial value for parameter latent0 in 'startpara'")

  assert_length(startlatent, length(y), "Argument 'startlatent'")
  assert_numeric(startlatent, "Argument 'startlatent'")
}

validate_expert <- function (expert) {
  ### joint arguments
  assert_logical(expert$correct_model_misspecification, "Expert argument 'correct_model_misspecification'")
  assert_single(expert$correct_model_misspecification, "Expert argument 'correct_model_misspecification'")

  assert_logical(expert$interweave, "Expert argument 'interweave'")
  assert_single(expert$interweave, "Expert argument 'interweave'")

  ### fast SV arguments UNDOCUMENTED (expert use!)
  assert_element(expert$fast_sv$baseline_parameterization, c("centered", "noncentered"),
                 "Fast SV expert argument 'baseline_parameterization'",
                 "the allowed values")
  assert_single(expert$fast_sv$baseline_parameterization,
                "Fast SV expert argument 'baseline_parameterization'")

  assert_element(expert$fast_sv$proposal_phi, c("immediate acceptance-rejection", "repeated acceptance-rejection"),
                 "Fast SV expert argument 'proposal_phi'",
                 "the allowed values")
  assert_single(expert$fast_sv$proposal_phi,
                "Fast SV expert argument 'proposal_phi'")

  assert_element(expert$fast_sv$proposal_sigma2, c("independence", "log random walk"),
                 "Fast SV expert argument 'proposal_sigma2'",
                 "the allowed values")
  assert_single(expert$fast_sv$proposal_sigma2,
                "Fast SV expert argument 'proposal_sigma2'")

  assert_positive(expert$fast_sv$proposal_intercept_var,
                  "Fast SV expert argument 'proposal_intercept_var'")
  assert_single(expert$fast_sv$proposal_intercept_var,
                "Fast SV expert argument 'proposal_intercept_var'")

  assert_positive(expert$fast_sv$proposal_phi_var,
                  "Fast SV expert argument 'proposal_phi_var'")
  assert_single(expert$fast_sv$proposal_phi_var,
                "Fast SV expert argument 'proposal_phi_var'")

  assert_positive(expert$fast_sv$proposal_sigma2_rw_scale,
                  "Fast SV expert argument 'proposal_sigma2_rw_scale'")
  assert_single(expert$fast_sv$proposal_sigma2_rw_scale,
                "Fast SV expert argument 'proposal_sigma2_rw_scale'")

  assert_element(expert$fast_sv$mh_blocking_steps, 1:3,
                 "Fast SV expert argument 'mh_blocking_steps'",
                 "the allowed values")
  assert_single(expert$fast_sv$mh_blocking_steps,
                "Fast SV expert argument 'mh_blocking_steps'")

  assert_logical(expert$fast_sv$store_indicators,
                 "Fast SV expert argument 'store_indicators'")
  assert_single(expert$fast_sv$store_indicators,
                "Fast SV expert argument 'store_indicators'")

  assert_element(expert$fast_sv$init_indicators, 1:10,
                 "Fast SV expert argument 'init_indicators'",
                 "the allowed values")

  ### general SV arguments UNDOCUMENTED (expert use!)
  assert_positive(expert$general_sv$multi_asis,
                  "General SV expert argument 'multi_asis'")
  assert_single(expert$general_sv$multi_asis,
                "General SV expert argument 'multi_asis'")

  assert_element(expert$general_sv$starting_parameterization, c("centered", "noncentered"),
                 "General SV expert argument 'starting_parameterization'",
                 "the allowed values")
  assert_single(expert$general_sv$starting_parameterization,
                "General SV expert argument 'starting_parameterization'")

  assert_element(expert$general_sv$proposal_para, c("random walk", "metropolis-adjusted langevin algorithm"),
                 "General SV expert argument 'proposal_para'",
                 "among the allowed values")
  assert_single(expert$general_sv$proposal_para,
                "General SV expert argument 'proposal_para'")

  if (isTRUE(is.list(expert$general_sv$proposal_diffusion_ken))) {
    assert_length(expert$general_sv$proposal_diffusion_ken, 2,
                   "General SV expert argument 'proposal_diffusion_ken' (the second moment of the random proposal)")
    assert_element(names(expert$general_sv$proposal_diffusion_ken), c("scale", "covariance"),
                   "General SV expert argument 'proposal_diffusion_ken' (the second moment of the random proposal)",
                   "among the allowed values")

    assert_single(expert$general_sv$proposal_diffusion_ken$scale,
                  "General SV expert argument 'proposal_diffusion_ken$scale' (the scaling for the proposal covariance matrix for the model parameters)")
    assert_positive(expert$general_sv$proposal_diffusion_ken$scale,
                    "General SV expert argument 'proposal_diffusion_ken$scale' (the scaling for the proposal covariance matrix for the model parameters)")

    if (!isTRUE(is.matrix(expert$general_sv$proposal_diffusion_ken$covariance))) {
      stop("General SV expert argument 'proposal_diffusion_ken$covariance' (the unscaled proposal covariance matrix for the model parameters) should be a matrix.")
    }
    covdims <- dim(expert$general_sv$proposal_diffusion_ken$covariance)
    if (any(covdims != 4)) {
      stop("General SV expert argument 'proposal_diffusion_ken$covariance' (the unscaled proposal covariance matrix for the model parameters) should be a 4x4 matrix; got dimensions ", covdims, ".")
    }
    tryCatch({
      chol(expert$general_sv$proposal_diffusion_ken$covariance)
    }, error = function (e) {
      stop("General SV expert argument 'proposal_diffusion_ken$covariance' (the unscaled proposal covariance matrix for the model parameters) should be a covariance matrix; cholesky factorization failed.")
    })
  } else if (!isTRUE(is.null(expert$general_sv$proposal_diffusion_ken))) {
    stop("General SV expert argument 'proposal_diffusion_ken' (the second moment of the random proposal) should be either NULL or a list with elements 'scale' and 'covariance'; received type ", typeof(expert$general_sv$proposal_diffusion_ken))
  }
}

