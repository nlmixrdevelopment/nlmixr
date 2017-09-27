#include <R.h>
#include <Rinternals.h>
#include <Rmath.h> //Rmath includes math.
#include <R_ext/Rdynload.h>

SEXP _nlmixr_chkSolvedInf(SEXP evid,SEXP solved){
  double *ev = REAL(evid);
  int n = length(evid);
  int i;
  int isSolved = INTEGER(solved)[0];
  int isInf = -1;
  SEXP out = PROTECT(allocVector(LGLSXP,1));
  // For a proper solved object ALL EVIDs must be 0,
  for (i = 0; i < n; i++){
    if (ev[i] == 0){
    } else if (fabs(ev[i]) < 101){
      error("Non-compatible EVID found(EVID=%d); Please check dataset!",ev[i]);
    } else if (isSolved){
      // All EVIDs need to be 101 or 10101
      if (isInf == -1){
      } else if (isInf == 0){
	if (ev[i] != 101){
	  error("With EVID=101, all EVIDs have to be 101 in a solved system.");
	}
      } else if (isInf == 0){
	error("With EVID=10101, all EVIDs have to be 10101 in a solved system.");
      }
    }
  }
  if (isSolved){
    LOGICAL(out)[0] = isInf;
  } else {
    LOGICAL(out)[0] = NA_LOGICAL;
  }
  UNPROTECT(1);
  return out;
}
