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
  int cid;
  for (i = 0; i < n; i++){
    cid = (int)ev[i];
    if (cid == 0){
    } else if (abs(cid) < 101){
      error("Non-compatible EVID found(EVID=%d); Please check dataset!",cid);
    } else if (isSolved){
      // All EVIDs need to be 101 or 10101
      if (isInf == -1){
	if (cid == 101){
	  isInf = 0;
	} else if (cid == 10101) {
	  isInf = 1;
	} else {
	  error("Solved systems can only have EVID=101 or EVID=10101; (found EVID=%d)",cid);
	}
      } else if (isInf == 0){
	if (cid != 101){
	  error("With EVID=101, all EVIDs have to be 101 in a solved system.");
	}
      } else if (isInf == 0){
	error("With EVID=10101, all EVIDs have to be 10101 in a solved system.");
      }
    }
  }
  if (isSolved){
    LOGICAL(out)[0] = (isInf==1);
  } else {
    LOGICAL(out)[0] = NA_LOGICAL;
  }
  UNPROTECT(1);
  return out;
}


SEXP _nlmixr_chkSortIDTime(SEXP _id,SEXP _time){
  double *time = REAL(_time);
  int *id = INTEGER(_id);
  int len = length(_time), i = 0, lastID;
  double lastTime;
  SEXP out = PROTECT(allocVector(INTSXP,1));
  INTEGER(out)[0] = 1;
  if (len != length(_id)){
    error("TIME and ID need the same length.");
  } else{
    for (i=0; i < len; i++){
      if (i == 0){
	lastID = id[i];
	lastTime = time[i];
      } else if (lastID==id[i]){
	if (lastTime > time[i]){
	  INTEGER(out)[0] = 0;
	  UNPROTECT(1);
	  return out;
	} else {
	  lastTime = time[i];
	}
      } else if (lastID < id[i]){
	INTEGER(out)[0] = 2;
	lastID = id[i];
	lastTime = time[i];
      } else {
	INTEGER(out)[0] = 0;
	UNPROTECT(1);
        return out;
      }
    }
  }
  UNPROTECT(1);
  return out;
}
