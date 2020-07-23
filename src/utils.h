#ifndef _STOCHVOL_UTILS_H_
#define _STOCHVOL_UTILS_H_

#include <R.h>

inline
void R_assert(const bool error_statement, const char* message) {
#ifndef NDEBUG
  if (error_statement) {
    Rf_error(message);
  }
#endif
  return;
}

#endif

