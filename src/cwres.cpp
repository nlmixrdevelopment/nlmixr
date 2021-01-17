#include "cwres.h"
static inline void calculateCwresDerr(arma::mat& fppm, arma::mat& fpim,
				      arma::ivec& ID, arma::mat &etas,
				      arma::vec &dErr_dEta_i, arma::vec &dErr_dEta_p,
				      List &etasDfFull, int &nid, unsigned int &neta) {
  int lastId = ID[ID.size()-1], lastCol = nid-1, lastIndex=ID.size()-1;
  int etaFulli = nid-1;
  double curEta=0.0;
  for (unsigned int j = fppm.n_rows; j--; ){
    if (lastId != ID[j]){
      // Fill in full eta data frame
      for (unsigned int i = neta; i--;){
	curEta = etas(etaFulli, i);//(as<NumericVector>(etasDf1[i]))[etaFulli];
	NumericVector cur = etasDfFull[i];
	std::fill_n(cur.begin()+j+1,lastIndex-j,curEta);
      }
      etaFulli--;
      // FIXME do it without copy?
      arma::vec tmp = fppm.rows(j+1, lastIndex) * trans(etas.row(lastCol));
      std::copy(tmp.begin(),tmp.end(),dErr_dEta_p.begin()+j+1);
      tmp = fpim.rows(j+1, lastIndex) * trans(etas.row(lastCol));
      std::copy(tmp.begin(),tmp.end(),dErr_dEta_i.begin()+j+1);
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
	// Finalize dErr_dEta
	arma::vec tmp = fppm.rows(0, lastIndex) * trans(etas.row(lastCol));
	std::copy(tmp.begin(),tmp.end(),dErr_dEta_p.begin());
	tmp = fpim.rows(0, lastIndex) * trans(etas.row(lastCol));
	std::copy(tmp.begin(),tmp.end(),dErr_dEta_i.begin());
	break;
      }
    }
  }
}

extern "C" SEXP _nlmixr_cwresCalc(SEXP ipredPredListSEXP, SEXP omegaMatSEXP,
				  SEXP etasDfSEXP, SEXP dvIn, SEXP evidIn, SEXP censIn, SEXP limitIn,
				  SEXP relevantLHSSEXP, SEXP stateSXP, SEXP cwresOpt) {
BEGIN_RCPP
  List ipredPredList = as<List>(ipredPredListSEXP);
  if (ipredPredList.size() != 4) return R_NilValue; //Rcpp::stop("malformed cwres calc");
  List ipredL = ipredPredList[0];
  List predL = ipredPredList[1];
  List etaLst = ipredPredList[2];
  List ebeL = ipredPredList[3];
  int ncalc = Rf_length(ipredL[0]);
  List etasDf = as<List>(etasDfSEXP);
  int nid = Rf_length(etasDf[0]);
  int npred = getPredIndex(ipredL);

  arma::vec ipredt(REAL(ipredL[npred]), ncalc, false, true);
  arma::vec ipred(ipredt.size());

  arma::vec predt(REAL(predL[npred]), ncalc, false, true);

  arma::vec dv(REAL(dvIn), ncalc, false, true);
  arma::vec dvt(ncalc);

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
  int hasLimit=0;
  getLimitFromInput(limitIn, ncalc, limit, hasLimit);

  arma::vec     hi(REAL(ipredL[ipredL.size()-1]), ncalc, false, true);
  arma::vec    low(REAL(ipredL[ipredL.size()-2]), ncalc, false, true);
  arma::vec     yj(REAL(ipredL[ipredL.size()-3]), ncalc, false, true);
  arma::vec lambda(REAL(ipredL[ipredL.size()-4]), ncalc, false, true);
  arma::vec lowerLim(ncalc);
  arma::vec upperLim(ncalc);

  arma::mat omegaMat = as<arma::mat>(omegaMatSEXP);
  unsigned int neta = omegaMat.n_rows;

  arma::vec rpv(REAL(predL[npred+1+neta]), ncalc, false, true);
  arma::vec riv(REAL(ipredL[npred+1+neta]), ncalc, false, true);

  bool doSim = true;
  List opt = as<List>(cwresOpt);
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
  								   riv, doSim, censMethod);


  arma::ivec ID(INTEGER(predL[0]), ncalc, false, true);

  arma::mat fppm(ncalc,neta);
  arma::mat fpim(ncalc,neta);

  arma::mat etas(nid, neta);
  List etasDfFull(neta);
  //etasDfFull.names()=etasDf1.names();
  CharacterVector etaN1 = etasDf.names();
  CharacterVector etaN2(neta);
  for (unsigned int j = neta; j--;) {
    fppm.col(j) = arma::vec(REAL(predL[j + 1 + npred]), ncalc, false, true);
    fpim.col(j) = arma::vec(REAL(ipredL[j + 1 + npred]), ncalc, false, true);
    etas.col(j) = arma::vec(REAL(etasDf[j+1]), nid, false, true);
    etaN2[j] = etaN1[j+1];
    etasDfFull[j] = NumericVector(ncalc);
  }
  etasDfFull.names() = etaN2;
  etasDfFull.attr("row.names")=IntegerVector::create(NA_INTEGER,-ncalc);
  etasDfFull.attr("class") = "data.frame";

  arma::mat V_fo_p = (fppm * omegaMat * fppm.t()); // From Mentre 2006 p. 352
  arma::mat V_fo_i = (fpim * omegaMat * fpim.t()); // From Mentre 2006 p. 352
  // There seems to be a difference between how NONMEM and R/S types
  // of software calculate WRES.  Mentre 2006 states that the
  // Variance under the FO condition should only be diag(Vfo_full) + Sigma,
  // but Hooker 2007 claims there is a
  // diag(Vfo_full)+diag(dh/deta*Sigma*dh/deta).
  // h = the additional error from the predicted function.
  //
  // In the nlmixr/FOCEi implemented here, the variance of the err
  // term is 1, or Sigma is a 1 by 1 matrix with one element (1)
  //
  // The dh/deta term would be the sd term, or sqrt(r), which means
  // sqrt(r)*sqrt(r)=|r|.  Since r is positive, this would be r.
  //
  // Also according to Hooker, WRES is calculated under the FO
  // assumption, where eta=0, eps=0 for this r term and Vfo term.
  // However, conditional weighted residuals are calculated under
  // the FOCE condition for the Vfo and the FO conditions for
  // dh/deta
  //

  arma::vec Vfop = V_fo_p.diag();
  arma::vec Vfoi = V_fo_i.diag();

  arma::vec dErr_dEta_i(ncalc);
  arma::vec dErr_dEta_p(ncalc);
  calculateCwresDerr(fppm, fpim, ID, etas, dErr_dEta_i, dErr_dEta_p, etasDfFull, nid, neta);

  arma::vec rest = dvt - predt;
  arma::vec wres = rest/sqrt(abs(Vfop+rpv));

  arma::vec cpredt = ipredt - dErr_dEta_i;
  arma::vec crest = dvt - cpredt;

  arma::vec cpred(cpredt.size());
  arma::vec cres(crest.size());
  arma::vec pred(predt.size());
  for (unsigned int i = cres.size(); i--;){
    cpred[i] = _powerDi(cpredt[i], lambda[i], (int)yj[i], low[i], hi[i]);
    cres[i]  = dv[i] - cpred[i];
    pred[i] = _powerDi(predt[i], lambda[i], (int)yj[i], low[i], hi[i]);
  }
  arma::vec res = dv - pred;
  arma::vec cwres = crest/sqrt(Vfoi+riv);
  arma::vec iwres=(dvt-ipredt)/sqrt(riv);
  arma::vec ires = dv - ipred;

  for (unsigned int j = ires.size(); j--; ) {
    if (censMethod == CENS_CPRED && cens[j] != 0) {
      // dvt[j]    = cpredt[j];
      // dv[j]	= cpred[j];
      // ires[j]	= dv[j] - ipred[j];
      // iwres[j]	= (dvt[j] - ipredt[j])/riv[j];
    } else if (censMethod == CENS_PRED && cens[j] != 0) {
      // dvt[j]    = predt[j];
      // dv[j]	= pred[j];
      // res[j]	= 0.0;
      // ires[j]	= dv[j] - ipred[j];
      // iwres[j]	= (dvt[j] - ipredt[j])/riv[j];
    } else if (censMethod == CENS_OMIT && cens[j] != 0) {
      dv[j]	= NA_REAL;
      pred[j]	= NA_REAL;
      res[j]	= NA_REAL;
      wres[j]	= NA_REAL;
      ipred[j]	= NA_REAL;
      ires[j]	= NA_REAL;
      iwres[j]	= NA_REAL;
      cpred[j]	= NA_REAL;
      cres[j]	= NA_REAL;
      cwres[j]	= NA_REAL;
    } else if (evid[j] != 0) {
      dv[j]	= NA_REAL;
      res[j]	= NA_REAL;
      wres[j]	= NA_REAL;
      ires[j]	= NA_REAL;
      iwres[j]	= NA_REAL;
      cres[j]	= NA_REAL;
      cwres[j]	= NA_REAL;
    }
  }
  int ncol = 9;
  if (interestingLimits) {
    ncol += 3 + hasLimit;
  }
  List retDF(ncol);
  CharacterVector nm(ncol);
  int i=0;
  //nm[i] = "DV"; retDF[i++] = wrap(dv);
  nm[i] = "PRED"; retDF[i++] = wrap(pred);
  nm[i] = "RES"; retDF[i++] = wrap(res);
  nm[i] = "WRES"; retDF[i++] = wrap(wres);
  nm[i] = "IPRED"; retDF[i++] = wrap(ipred);
  nm[i] = "IRES"; retDF[i++] = wrap(ires);
  nm[i] = "IWRES"; retDF[i++] = wrap(iwres);
  nm[i] = "CPRED"; retDF[i++] = wrap(cpred);
  nm[i] = "CRES"; retDF[i++] = wrap(cres);
  nm[i] = "CWRES"; retDF[i++] = wrap(cwres);
  if (interestingLimits) {
    nm[i] = "CENS"; retDF[i++] = wrap(cens);
    if (hasLimit){
      nm[i] = "LIMIT"; retDF[i++] = wrap(limit);
    }
    nm[i] = "lowerLim"; retDF[i++] = wrap(lowerLim);
    nm[i] = "upperLim"; retDF[i++] = wrap(upperLim);
  }
  retDF.names() = nm;
  retDF.attr("row.names") = IntegerVector::create(NA_INTEGER,-ncalc);
  retDF.attr("class") = "data.frame";
  calcShrinkFinalize(omegaMat, nid, etaLst, iwres, evid, etaN2, 1);
  List retC = List::create(retDF, etasDfFull, getDfSubsetVars(ipredL, stateSXP),
			   getDfSubsetVars(ebeL, relevantLHSSEXP));
  dfSetStateLhsOps(retC, opt);
  retC = dfCbindList(wrap(retC));
  List ret(4);
  ret[0] = getDfIdentifierCols(ipredL, npred);
  ret[1] = List::create(_["DV"] = wrap(dv));
  ret[2] = retC;
  ret[3] = etaLst;
  return wrap(ret);
END_RCPP
}
