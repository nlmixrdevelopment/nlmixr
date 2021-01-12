#include "censResid.h"

bool censTruncatedMvnReturnInterestingLimits(arma::vec& dv, arma::vec& dvt,
					     arma::vec& ipred, arma::vec &ipredt, 
					     arma::ivec &cens, arma::vec &limit,
					     arma::vec& lambda, arma::vec &yj, arma::vec& low, arma::vec& hi, 
					     arma::vec &lowerLim, arma::vec &upperLim, arma::vec &ri) {
  bool interestingLim = false;
  for (int i = dv.size(); i--;) {
    // ipredt is the transformation information
    ipred[i]= _powerDi(ipredt[i], lambda[i], (int)yj[i], low[i], hi[i]);
    switch (cens[i]){
    case 1:
      // (limit, dv) ; limit could be -inf
      if (R_FINITE(limit[i])){
	double lim0TBS = _powerD(limit[i], lambda[i], (int)yj[i], low[i], hi[i]);
	double lim1TBS = _powerD(dv[i], lambda[i], (int) yj[i], low[i], hi[i]);
	double sd = sqrt(ri[i]);
	interestingLim=true;
	lowerLim[i] = limit[i];
	upperLim[i] = dv[i];
	dvt[i] = truncnorm(ipredt[i], sd, lim0TBS, lim1TBS);
	dv[i] =_powerDi(dvt[i], lambda[i], (int) yj[i], low[i], hi[i]);
      } else {
	// (-Inf, dv)
	double lim1TBS = _powerD(dv[i], lambda[i], (int) yj[i], low[i], hi[i]);
	double sd = sqrt(ri[i]);
	interestingLim=true;
	lowerLim[i] = R_NegInf;
	upperLim[i] = dv[i];
	dvt[i] = truncnorm(ipredt[i], sd, R_NegInf, lim1TBS);
	dv[i] = _powerDi(dvt[i], lambda[i], (int) yj[i], low[i], hi[i]);
      }
      break;
    case -1:
      // (dv, limit); limit could be +inf
      if (R_FINITE(limit[i])){
	//(dv, limit)
	double lim1TBS = _powerD(limit[i], lambda[i], (int)yj[i], low[i], hi[i]);
	double lim0TBS = _powerD(dv[i], lambda[i], (int) yj[i], low[i], hi[i]);
	double sd = sqrt(ri[i]);
	interestingLim=true;
	lowerLim[i] = dv[i];
	upperLim[i] = limit[i];
	dvt[i] = truncnorm(ipredt[i], sd, lim0TBS, lim1TBS);
	dv[i] = _powerDi(dvt[i], lambda[i], (int) yj[i], low[i], hi[i]);
      } else {
	// (dv, Inf)
	double lim1TBS = _powerD(dv[i], lambda[i], (int) yj[i], low[i], hi[i]);
	double sd = sqrt(ri[i]);
	interestingLim=true;
	lowerLim[i] = dv[i];
	upperLim[i] = R_PosInf;
	dvt[i] = truncnorm(ipredt[i], sd, lim1TBS, R_PosInf);
	dv[i] = _powerDi(dvt[i], lambda[i], (int) yj[i], low[i], hi[i]);
      }
      break;
    case 0:
      lowerLim[i] = NA_REAL;
      upperLim[i] = NA_REAL;
      dvt[i] = _powerD(dv[i], lambda[i], (int)yj[i], low[i], hi[i]);
      break;
    }
  }
  return interestingLim;
}
