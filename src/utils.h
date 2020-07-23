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

}

#endif

