/*
 * R package stochvol by
 *     Gregor Kastner Copyright (C) 2013-2020
 *     Darjus Hosszejni Copyright (C) 2019-2020
 *  
 *  This file is part of the R package stochvol: Efficient Bayesian
 *  Inference for Stochastic Volatility Models.
 *  
 *  The R package stochvol is free software: you can redistribute it
 *  and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation, either version 2 or
 *  any later version of the License.
 *  
 *  The R package stochvol is distributed in the hope that it will be
 *  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 *  General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with the R package stochvol. If that is not the case, please
 *  refer to <http://www.gnu.org/licenses/>.
 */

/*
 * utils.h
 * 
 * Utility functions that don't fit anywhere else.
 */

#ifndef _STOCHVOL_UTILS_H_
#define _STOCHVOL_UTILS_H_

#include <R.h>
#include <type_definitions.h>

namespace stochvol {

// Our own assert function that throws R exceptions
inline
void R_assert(const bool assert_statement, const char* message) {
#ifndef NDEBUG
  if (!assert_statement) {
    Rf_error(message);
  }
#endif
  return;
}

// Determine the inverse variance of h0 in the
// noncentered parameterization.
inline
double determine_Bh0inv(
    const double phi,
    const PriorSpec& prior_spec) {
  switch (prior_spec.latent0.variance) {
    case PriorSpec::Latent0::STATIONARY:
      return 1. - std::pow(phi, 2);
    case PriorSpec::Latent0::CONSTANT:
      return 1. / prior_spec.latent0.constant.value;
  }
}

// Determine parameter 'thintime'
inline
int determine_thintime(
    const int T,
    const Rcpp::CharacterVector& keeptime_in) {
  const std::string keeptime = Rcpp::as<std::string>(keeptime_in);
  if (keeptime == "all") {
    return 1;
  } else if (keeptime == "last") {
    return T;
  } else {
    Rf_error("Unknown value for 'keeptime'; got \"%s\"", keeptime.c_str());
  }
}

// Transform the latent vector from centered to
// noncentered.
template<typename T>
T centered_to_noncentered(
    const double mu,
    const double sigma,
    const T& h) {
  return (h - mu) / sigma;
}

// Transform the latent vector from noncentered to
// centered.
template<typename T>
T noncentered_to_centered(
    const double mu,
    const double sigma,
    const T& ht) {
  return mu + sigma * ht;
}

}

#endif

