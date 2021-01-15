#ifndef __RES_H__
#define __RES_H__

#if defined(__cplusplus)
#include "armahead.h"
#include "censResid.h"
#include "shrink.h"
void calculateDfFull(arma::ivec& ID, arma::mat &etas,
		     List &etasDfFull, int &nid, unsigned int &neta);

extern "C" {
#endif
  SEXP _nlmixr_resCalc(SEXP ipredPredListSEXP, SEXP omegaMatSEXP,
		       SEXP etasDfSEXP, SEXP dvIn, SEXP evidIn, SEXP censIn, SEXP limitIn,
		       SEXP resOpt);
#if defined(__cplusplus)
}
#endif
#endif
