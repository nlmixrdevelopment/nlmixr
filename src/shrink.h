#ifndef __SHRINK_H__
#define __SHRINK_H__

#if defined(__cplusplus)
#include "armahead.h"

void calcShrinkFinalize(arma::mat &omegaMat, int &nid, List& etaLst, arma::vec &iwres, arma::ivec &evid,
			CharacterVector &etaNames, int doIwres);

#endif

SEXP _nlmixr_calcShrinkOnly(SEXP omegaMatSEXP, SEXP etaLstSEXP);

#endif
