#ifndef _STOCHVOL_UTILS_H_
#define _STOCHVOL_UTILS_H_

#include <R.h>

namespace stochvol {

inline
void R_assert(const bool assert_statement, const char* message) {
#ifndef NDEBUG
  if (!assert_statement) {
    Rf_error(message);
  }
#endif
  return;
}

template<typename T>
T centered_to_noncentered(
    const double mu,
    const double sigma,
    const T& h) {
  return (h - mu) / sigma;
}

template<typename T>
T noncentered_to_centered(
    const double mu,
    const double sigma,
    const T& ht) {
  return mu + sigma * ht;
}

}

#endif

