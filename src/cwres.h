#ifndef __CWRES_H__
#define __CWRES_H__
#if defined(__cplusplus)
#include "armahead.h"
#include "censResid.h"
#include "shrink.h"
#include "res.h"
#include "utilc.h"

extern "C" {
#endif
  SEXP _nlmixr_cwresCalc(SEXP ipredPredListSEXP, SEXP omegaMatSEXP,
			 SEXP etasDfSEXP, SEXP dvIn, SEXP evidIn, SEXP censIn, SEXP limitIn,
			 SEXP relevantLHSSEXP, SEXP stateSXP, SEXP covSXP, SEXP cwresOpt);
  
#if defined(__cplusplus)
}
#endif
#endif
