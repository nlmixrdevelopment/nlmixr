#include <sys/stat.h> 
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>   /* dj: import intptr_t */
#include <errno.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("RxODE", String)
/* replace pkg as appropriate */
#else
#define _(String) (String)
#endif

int _setSilentErr=0;
extern void setSilentErr(int silent){
  _setSilentErr = silent;
}

SEXP _nlmixr_setSilentErr(SEXP in) {
  SEXP ret = PROTECT(Rf_allocVector(LGLSXP, 1));
  int t = TYPEOF(in);
  if (Rf_length(in) > 0) {
    if (t == INTSXP) {
      if (INTEGER(in)[0] > 0) {
	_setSilentErr = 1;
	INTEGER(ret)[0] = 1;
	UNPROTECT(1);
	return ret;
      } else {
	_setSilentErr = 0;
	INTEGER(ret)[0] = 0;
	UNPROTECT(1);
	return ret;
      }
    } else if (t == LGLSXP) {
      if (INTEGER(in)[0] > 0) {
	_setSilentErr = 1;
	INTEGER(ret)[0] = 1;
	UNPROTECT(1);
	return ret;
      } else {
	_setSilentErr = 0;
	INTEGER(ret)[0] = 0;
	UNPROTECT(1);
	return ret;
      }
    } else if (t == REALSXP) {
      if (REAL(in)[0] > 0) {
	_setSilentErr = 1;
	INTEGER(ret)[0] = 1;
	UNPROTECT(1);
	return ret;
      } else {
	_setSilentErr = 0;
	INTEGER(ret)[0] = 0;
	UNPROTECT(1);
	return ret;
	return R_NilValue;
      }
    }
  }
  _setSilentErr = 0;
  INTEGER(ret)[0] = 0;
  UNPROTECT(1);
  return ret;
}

extern void RSprintf(const char *format, ...) {
  if (_setSilentErr == 0) {
    va_list args;
    va_start(args, format);
    Rvprintf(format, args);
    va_end(args);
  }
}
