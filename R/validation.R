#  #####################################################################################
#  R package stochvol by
#     Gregor Kastner Copyright (C) 2013-2018
#     Gregor Kastner and Darjus Hosszejni Copyright (C) 2019-
#  
#  This file is part of the R package stochvol: Efficient Bayesian
#  Inference for Stochastic Volatility Models.
#  
#  The R package stochvol is free software: you can redistribute it
#  and/or modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation, either version 2 or
#  any later version of the License.
#  
#  The R package stochvol is distributed in the hope that it will be
#  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#  General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with the R package stochvol. If that is not the case, please
#  refer to <http://www.gnu.org/licenses/>.
#  #####################################################################################

# Functions that assert some properties of the inputs

sv_stop <- function (what, name, should, and_not) {
  stop(name, " = ", prettify(what), " should ", should, "; not ", and_not)
}

prettify <- function (x) {
  if (is.character(x)) {
    x <- paste0("\"", x, "\"")
  }
  out <- paste0("c(", paste(x, collapse = ", "), ")")
  max_length <- 80L
  if (nchar(out) > max_length) {
    out <- paste(substr(out, start = 1, stop = max_length), "...)")
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
  assert_element(keeptime, c("all", "last"), name_keeptime, "'all' or 'last'")
}

validate_initial_values <- function (startpara, startlatent, y, x) {
  if (!isTRUE(is.list(startpara))) stop("Failed validation: startpara should be a list.")
  if (!isTRUE(all(sapply(startpara, is.list)))) stop("Failed validation: startpara should be a list of lists.")

  if (!isTRUE(is.list(startlatent))) stop("Failed validation: startlatent should be a list.")
  if (!isTRUE(all(sapply(startlatent, length) == length(startlatent[[1]])))) stop("Failed validation: elements of startlatent should be of same length.")

  lapply(startpara, function (startpara_i, designmatrix) {
    assert_element(names(startpara_i), c("mu", "phi", "sigma", "nu", "rho", "beta", "latent0"),
                   "Elements of argument 'startpara'", "parameters of the model")

    assert_single(startpara_i$mu, "Provided initial value for parameter mu in 'startpara'")
    assert_numeric(startpara_i$mu, "Provided initial value for parameter mu in 'startpara'")

    assert_single(startpara_i$phi, "Provided initial value for parameter phi in 'startpara'")
    assert_gt(startpara_i$phi, -1, "Provided initial value for parameter phi in 'startpara'")
    assert_lt(startpara_i$phi, 1, "Provided initial value for parameter phi in 'startpara'")

    assert_single(startpara_i$sigma, "Provided initial value for parameter sigma in 'startpara'")
    assert_positive(startpara_i$sigma, "Provided initial value for parameter sigma in 'startpara'")

    assert_single(startpara_i$nu, "Provided initial value for parameter nu in 'startpara'")
    assert_gt(startpara_i$nu, 2, "Provided initial value for parameter nu in 'startpara'")

    assert_single(startpara_i$rho, "Provided initial value for parameter rho in 'startpara'")
    assert_gt(startpara_i$rho, -1, "Provided initial value for parameter rho in 'startpara'")
    assert_lt(startpara_i$rho, 1, "Provided initial value for parameter rho in 'startpara'")

    assert_length(startpara_i$beta, NCOL(designmatrix), "Provided initial values for the regressors beta in 'startpara'")
    assert_numeric(startpara_i$beta, "Provided initial values for the regressors beta in 'startpara'")

    assert_single(startpara_i$latent0, "Provided initial value for parameter latent0 in 'startpara'")
    assert_numeric(startpara_i$latent0, "Provided initial value for parameter latent0 in 'startpara'")
  }, designmatrix = x)

  lapply(startlatent, function (startlatent_i, observations) {
    assert_length(startlatent_i, length(observations), "Argument 'startlatent'")
    assert_numeric(startlatent_i, "Argument 'startlatent'")
  }, observations = y)

  invisible(NULL)
}

validate_expert <- function (expert) {
  ### joint arguments
  assert_logical(expert$correct_model_misspecification, "Expert argument 'correct_model_misspecification'")
  assert_single(expert$correct_model_misspecification, "Expert argument 'correct_model_misspecification'")

  assert_logical(expert$interweave, "Expert argument 'interweave'")
  assert_single(expert$interweave, "Expert argument 'interweave'")

  ### fast SV arguments
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

  assert_logical(expert$fast_sv$update$parameters,
                 "Fast SV expert argument 'update$parameters'")
  assert_single(expert$fast_sv$update$parameters,
                "Fast SV expert argument 'update$parameters'")

  assert_logical(expert$fast_sv$update$latent_vector,
                 "Fast SV expert argument 'update$latent_vector'")
  assert_single(expert$fast_sv$update$parameters,
                "Fast SV expert argument 'update$latent_vector'")

  assert_logical(expert$fast_sv$update$mixture_indicators,
                 "Fast SV expert argument 'update$mixture_indicators'")
  assert_single(expert$fast_sv$update$mixture_indicators,
                "Fast SV expert argument 'update$mixture_indicators'")

  assert_element(expert$fast_sv$init_indicators, 1:10,
                 "Fast SV expert argument 'init_indicators'",
                 "the allowed values")

  assert_positive(expert$fast_sv$init_tau,
                 "Fast SV expert argument 'init_tau'")

  ### general SV arguments
  assert_positive(expert$general_sv$multi_asis,
                  "General SV expert argument 'multi_asis'")
  assert_single(expert$general_sv$multi_asis,
                "General SV expert argument 'multi_asis'")

  if (!is.null(expert$general_sv$theta_asis_setup)) {
    assert_length(expert$general_sv$theta_asis_setup, 3,
                  "General SV expert argument 'theta_asis_setup'")
  }

  assert_length(expert$general_sv$nu_asis_setup, 3,
                "General SV expert argument 'nu_asis_setup'")

  assert_element(expert$general_sv$starting_parameterization, c("centered", "noncentered"),
                 "General SV expert argument 'starting_parameterization'",
                 "the allowed values")
  assert_single(expert$general_sv$starting_parameterization,
                "General SV expert argument 'starting_parameterization'")

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
  } else if (!isTRUE(is.logical(expert$general_sv$proposal_diffusion_ken)) && !expert$general_sv$proposal_diffusion_ken) {
    stop("General SV expert argument 'proposal_diffusion_ken' (the second moment of the random proposal) should be either FALSE or a list with elements 'scale' and 'covariance'; received type ", typeof(expert$general_sv$proposal_diffusion_ken))
  }

  assert_positive(expert$general_sv$init_tau,
                 "General SV expert argument 'init_tau'")

  invisible(NULL)
}

#' Validate and Process Argument 'expert'
#' 
#' A helper function that validates the input and extends it with
#' default values if there are missing parts for argument 'expert'.
#' @param expert list, the input values for expert.
#' @param priorspec a \code{priorspec} object created by
#' \code{\link{specify_priors}}
#' @return A list that is the input extended by default values. If
#' the input is invalid, an error is thrown.
#' @family validation
#' @seealso \code{\link{specify_priors}}
#' @export
validate_and_process_expert <- function (expert = NULL, priorspec = specify_priors()) {
  expertdefault <-
    list(correct_model_misspecification = FALSE,  # online correction for general_sv and post-correction for fast_sv
         interweave = TRUE,
         fast_sv =  # UNDOCUMENTED; very expert settings of the fast_sv sampler
           get_default_fast_sv(),
         general_sv =  # UNDOCUMENTED; very expert settings of the general_sv sampler
           get_default_general_sv(priorspec))
  expert <- apply_default_list(expert, expertdefault, "Names in expert", "allowed names in expert")
  validate_expert(expert)
  expert
}

