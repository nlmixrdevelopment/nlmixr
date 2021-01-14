#include "res.h"

void calculateDfFull(arma::ivec& ID, arma::mat &etas,
		     List &etasDfFull, int &nid, unsigned int &neta) {
  int lastId = ID[ID.size()-1], lastCol = nid-1, lastIndex=ID.size()-1;
  int etaFulli = nid-1;
  double curEta=0.0;
  for (unsigned int j = ID.size(); j--; ){
    if (lastId != ID[j]){
      // Fill in full eta data frame
      for (unsigned int i = neta; i--;){
	curEta = etas(etaFulli, i);//(as<NumericVector>(etasDf1[i]))[etaFulli];
	NumericVector cur = etasDfFull[i];
	std::fill_n(cur.begin()+j+1,lastIndex-j,curEta);
      }
      etaFulli--;
      lastId=ID[j];
      lastIndex=j;
      lastCol--;
      if (lastCol == 0){
	// Finalize ETA
	for (unsigned int i = neta; i--;){
	  curEta = etas(0, i);//(as<NumericVector>(etasDf1[i]))[0];
	  NumericVector cur = etasDfFull[i];
	  std::fill_n(cur.begin(),lastIndex+1,curEta);
	}
	break;
      }
    }
  }
}

extern "C" SEXP _nlmixr_resCalc(SEXP ipredPredListSEXP, SEXP omegaMatSEXP,
				SEXP etasDfSEXP, SEXP dvIn, SEXP evidIn, SEXP censIn, SEXP limitIn,
				SEXP resOpt) {
  List ipredPredList = as<List>(ipredPredListSEXP);
  if (ipredPredList.size() !=2) stop("malformed cwres calc");
  List ipredL = ipredPredList[0];
  List predL = ipredPredList[1];
  int ncalc = Rf_length(ipredL[0]);
  List etasDf = as<List>(etasDfSEXP);
  int nid = Rf_length(etasDf[0]);

  arma::vec ipredt(REAL(ipredL[2]), ncalc, false, true);
  arma::vec ipred(ipredt.size());

  arma::vec predt(REAL(predL[2]), ncalc, false, true);

  arma::vec dv(REAL(dvIn), ncalc, false, true);
  arma::vec dvt(ncalc);

  arma::vec rpv(REAL(predL[3]), ncalc, false, true);
  arma::vec riv(REAL(ipredL[3]), ncalc, false, true);


  arma::ivec cens;
  if (Rf_isNull(censIn)) {
    cens = arma::ivec(ncalc, fill::zeros);
  } else {
    cens = as<arma::ivec>(censIn);
  }
  arma::ivec evid;
  if (Rf_isNull(evidIn)) {
    evid = arma::ivec(ncalc, fill::zeros);
  } else {
    evid = as<arma::ivec>(evidIn);
  }
  arma::vec limit;
  if (Rf_isNull(limitIn)) {
    limit = arma::vec(ncalc);
  } else {
    limit = as<arma::vec>(limitIn);
  }

  arma::vec     hi(REAL(ipredL[ipredL.size()-1]), ncalc, false, true);
  arma::vec    low(REAL(ipredL[ipredL.size()-2]), ncalc, false, true);
  arma::vec     yj(REAL(ipredL[ipredL.size()-3]), ncalc, false, true);
  arma::vec lambda(REAL(ipredL[ipredL.size()-4]), ncalc, false, true);
  arma::vec lowerLim(ncalc);
  arma::vec upperLim(ncalc);
  
  arma::mat omegaMat = as<arma::mat>(omegaMatSEXP);
  unsigned int neta = omegaMat.n_rows;

  bool doSim = true;
  List opt = as<List>(resOpt);
  if (opt.containsElementNamed("doSim")) {
    RObject tmp = opt["doSim"];
    if (TYPEOF(tmp) == LGLSXP) {
      doSim = as<bool>(tmp);
    }
  }
  int censMethod = CENS_TNORM;
  if (opt.containsElementNamed("censMethod")) {
    RObject tmp = opt["censMethod"];
    if (TYPEOF(tmp) == INTSXP) {
      censMethod = as<int>(opt["censMethod"]);
    }
  }

  bool interestingLimits = censTruncatedMvnReturnInterestingLimits(dv, dvt, ipred, ipredt, cens, limit,
  								   lambda, yj, low, hi, lowerLim, upperLim,
  								   riv, doSim);
    

  arma::ivec ID(INTEGER(predL[0]), ncalc, false, true);

  arma::mat etas(nid, neta);
  List etasDfFull(neta);
  //etasDfFull.names()=etasDf1.names();
  CharacterVector etaN1 = etasDf.names();
  CharacterVector etaN2(neta);
  for (unsigned int j = neta; j--;) {
    etas.col(j) = arma::vec(REAL(etasDf[j+1]), nid, false, true);
    etaN2[j] = etaN1[j+1];
    etasDfFull[j] = NumericVector(ncalc);
  }
  etasDfFull.names() = etaN2;
  etasDfFull.attr("row.names")=IntegerVector::create(NA_INTEGER,-ncalc);
  etasDfFull.attr("class") = "data.frame";
  
  calculateDfFull(ID, etas, etasDfFull, nid, neta);

  arma::vec pred(predt.size());
  for (unsigned int i = pred.size(); i--;){
    pred[i] = _powerDi(predt[i], lambda[i], (int)yj[i], low[i], hi[i]);
  }
  
  arma::vec res = dv - pred;
  arma::vec iwres=(dvt-ipredt)/sqrt(riv);
  arma::vec ires = dv - ipred;
  
  for (unsigned int j = ires.size(); j--; ) {
    if (censMethod == CENS_OMIT && cens[j] != 0) {
      dv[j]	= NA_REAL;
      pred[j]	= NA_REAL;
      res[j]	= NA_REAL;
      ipred[j]	= NA_REAL;
      ires[j]	= NA_REAL;
      iwres[j]	= NA_REAL;
    } else if (evid[j] != 0) {
      dv[j]	= NA_REAL;
      res[j]	= NA_REAL;
      ires[j]	= NA_REAL;
      iwres[j]	= NA_REAL;
    }
  }
  DataFrame retDF;

  if (interestingLimits) {
    retDF = DataFrame::create(_["PRED"]=wrap(pred),
			      _["RES"]=wrap(res),
			      _["IPRED"]=wrap(ipred),
			      _["IRES"]=wrap(ires),
			      _["IWRES"]=wrap(iwres),
			      _["lowerLim"] = wrap(lowerLim),
			      _["upperLim"] = wrap(upperLim));
  } else {
    retDF = DataFrame::create(_["PRED"]=wrap(pred),
			      _["RES"]=wrap(res),
			      _["IPRED"]=wrap(ipred),
			      _["IRES"]=wrap(ires),
			      _["IWRES"]=wrap(iwres));
  }
  List ret(3);
  ret[0] = retDF;
  ret[1] = dv;
  ret[2] = etasDfFull;
  // ret[1] = etaLst;
  // ret[2] = etasDfFull;
  // ret[3] = dv;
  return wrap(ret);
}
