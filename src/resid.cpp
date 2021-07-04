#define STRICT_R_HEADER
#include "armahead.h"
#include "censResid.h"
using namespace R;

//[[Rcpp::export]]
RObject augPredTrans(NumericVector& pred, NumericVector& ipred, NumericVector& lambda,
		     RObject& yjIn, NumericVector& low, NumericVector& hi){
  IntegerVector yj = as<IntegerVector>(yjIn);
  for (int i = pred.size(); i--;){
    pred[i] = _powerDi(pred[i], lambda[i], yj[i], low[i], hi[i]);
    ipred[i] = _powerDi(ipred[i], lambda[i], yj[i], low[i] ,hi[i]);
  }
  return R_NilValue;
}
