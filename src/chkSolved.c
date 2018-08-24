#include <R.h>
#include <Rinternals.h>
#include <Rmath.h> //Rmath includes math.
#include <R_ext/Rdynload.h>

SEXP _nlmixr_allDose(SEXP evid, SEXP ids){
  SEXP outs = PROTECT(allocVector(LGLSXP,length(evid)));
  int *ev = INTEGER(evid);
  int *id = INTEGER(ids);
  int lastId = id[length(ids)-1]+1;
  int found = 0;
  int lastI = 0;
  for (unsigned int i = length(evid); i--;){
    if (lastId != id[i]){
      lastId = id[i];
      if (found == 0){
	for (unsigned int j = i+1; j < lastI+1; j++){
	  LOGICAL(outs)[j] = 0;
	}
      } else {
        for (unsigned int j = i+1; j < lastI+1; j++){
          LOGICAL(outs)[j] = 1;
	}
      }
      lastI = i;
      found = 0;
      if (ev[i] == 0) found=1;
    } else if (ev[i] == 0) {
      found=1;
    }
  }
  if (found == 0){
    for (unsigned int j = 0; j < lastI+1; j++){
      LOGICAL(outs)[j] = 0;
    }
  } else {
    for (unsigned int j = 0; j < lastI+1; j++){
      LOGICAL(outs)[j] = 1;
    }
  }
  UNPROTECT(1);
  return outs;
}

SEXP _nlmixr_convertEvid(SEXP evid, SEXP cmt){
  int *ev = INTEGER(evid);
  int *amt = INTEGER(cmt);
  SEXP outs = PROTECT(allocVector(INTSXP,length(evid)));
  int *out = INTEGER(outs);
  for (unsigned int i = length(evid); i--;){
    if (ev[i]){
      if (amt[i] > 99){
	int amt100 = amt[i]/100;
	int amt99  = amt[i]-amt100*100;
        out[i] = amt100*1e5+amt99*100+1;
      } else {
	out[i] = 100*amt[i]+1;
      }
    } else {
      out[i]=0;
    }
  }
  UNPROTECT(1);
  return(outs);
}

SEXP _nlmixr_convertEvidRate(SEXP evid, SEXP cmt){
  int *ev = INTEGER(evid);
  int *amt = INTEGER(cmt);
  SEXP outs = PROTECT(allocVector(INTSXP,length(evid)));
  int *out = INTEGER(outs);
  for (unsigned int i = length(evid); i--;){
    if (ev[i]){
      if (amt[i] > 99){
        int amt100 = amt[i]/100;
        int amt99  = amt[i]-amt100*100;
        out[i] = amt100*1e5+amt99*100+10001;
      } else {
        out[i] = 100*amt[i]+10001;
      }
    } else {
      out[i]=0;
    }
  }
  UNPROTECT(1);
  return(outs);
}

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


SEXP _nlmixr_chkSortIDTime(SEXP _id,SEXP _time, SEXP _evid){
  double *time = REAL(_time);
  int *id = INTEGER(_id);
  int *evid = INTEGER(_evid);
  int len = length(_time), i = 0, lastID=-1;
  double lastTime=-1;
  SEXP out = PROTECT(allocVector(INTSXP,1));
  INTEGER(out)[0] = 1;
  if (len != length(_id)){
    error("TIME and ID need the same length.");
  } else{
    for (i=0; i < len; i++){
      if (evid[i] != 0 && abs(evid[i]) < 101){
	INTEGER(out)[0] = 3;
        UNPROTECT(1);
	return out;
      }
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
