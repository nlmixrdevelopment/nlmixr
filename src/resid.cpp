#include "armahead.h"
#include "censResid.h"
using namespace R;

// [[Rcpp::export]]
List nlmixrShrink(NumericMatrix &omegaMat,DataFrame etasDf, List etaLst){
  unsigned int nid=etasDf.nrows();
  double tc = sqrt((double)nid);
  double om;
  int j;
  unsigned int neta = omegaMat.nrow();
  for (j = neta;j--;){
    NumericVector cur = etaLst[j];
    // Finalize by adding shrinkage.
    om =omegaMat(j,j);
    cur[5]= (1.0-cur[1]/om)*100.0;
    cur[6]= (1.0-cur[2]/sqrt(om))*100.0;
    cur[7] = cur[0]*tc/(cur[2]);
    cur[8] = 2*Rf_pt(-fabs(cur[7]),(double)(nid-1),1,0);
  }
  CharacterVector dimN= (as<List>(omegaMat.attr("dimnames")))[1];
  etaLst.attr("names") = dimN;
  etaLst.attr("row.names") = CharacterVector::create("mean","var","sd","skewness", "kurtosis","var shrinkage (%)", "sd shrinkage (%)","t statistic","p-value");
  etaLst.attr("class") = "data.frame";
  return etaLst;
}

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
