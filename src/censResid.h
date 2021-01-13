#ifndef __CENSRESID_H__
#define __CENSRESID_H__
#include "armahead.h"
#define CENS_OMIT 1 
#define CENS_CDF 2
#define CENS_TNORM 3

#if defined(__cplusplus)

static inline double truncnorm(double mean, double sd, double low, double hi){
  NumericMatrix sigma(1,1);
  sigma(0,0)=sd;
  SEXP ret =RxODE::rxRmvnSEXP(wrap(IntegerVector::create(1)),
			      wrap(NumericVector::create(mean)),
			      wrap(sigma),
			      wrap(NumericVector::create(low)),
			      wrap(NumericVector::create(hi)),
			      wrap(IntegerVector::create(1)),
			      wrap(LogicalVector::create(false)),
			      wrap(LogicalVector::create(false)),
			      wrap(NumericVector::create(0.4)),
			      wrap(NumericVector::create(2.05)),
			      wrap(NumericVector::create(1e-10)),
			      wrap(IntegerVector::create(100)));
  return REAL(ret)[0];
}

bool censTruncatedMvnReturnInterestingLimits(arma::vec& dv, arma::vec& dvt,
					     arma::vec& ipred, arma::vec &ipredt, 
					     arma::ivec &cens, arma::vec &limit,
					     arma::vec& lambda, arma::vec &yj, arma::vec& low, arma::vec& hi, 
					     arma::vec &lowerLim, arma::vec &upperLim, arma::vec &ri,
					     bool &doSim);

#endif

#endif
