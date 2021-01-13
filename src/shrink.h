#ifndef __SHRINK_H__
#define __SHRINK_H__
#include "armahead.h"

#if defined(__cplusplus)

void calcShrinkFinalize(arma::mat &omegaMat, unsigned int &nid, List& etaLst, arma::vec &iwres, arma::ivec &evid, CharacterVector &etaNames);

#endif

#endif
