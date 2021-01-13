#ifndef __RES_H__
#define __RES_H__
#include "armahead.h"
#include "censResid.h"

void calculateDfFull(arma::ivec& ID, arma::mat &etas,
		     List &etasDfFull, int &nid, unsigned int &neta);

extern "C" SEXP _nlmixr_cwresCalc(SEXP ipredPredListSEXP, SEXP omegaMatSEXP,
				  SEXP etasDfSEXP, SEXP dvIn, SEXP evidIn, SEXP censIn, SEXP limitIn,
				  SEXP cwresOpt);

#endif
