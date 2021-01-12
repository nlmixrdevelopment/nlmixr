#include "armahead.h"

arma::ivec getSimIdLoc(arma::ivec& id, arma::ivec& simId,
		      unsigned int &nid, unsigned int &K) {
  int i = 0;
  int lastId=id[0];
  int lastSimId = simId[0];
  int totNobs = 0;
  nid = 1;
  while (simId[i] == lastSimId) {
    if (lastId != id[i]) {
      nid++;
      lastId = id[i];
    }
    i++; totNobs++;
  }
  K =   id.size() /  totNobs;

  arma::ivec idLoc(nid+1);
  unsigned int j = 0;
  for (i = 0; i < totNobs; i++) {
    if (lastId != id[i]){
      idLoc[j++] = i;
      lastId = id[i];
    }
  }
  idLoc[j] = totNobs;
  return idLoc;
}

arma::mat getSimMatById(arma::ivec& idLoc, arma::vec &sim, unsigned int& id,
			unsigned int& K) {
  int nobs = idLoc[id+1]-idLoc[id];
  int totNobs = idLoc[idLoc.size()-1];
  arma::mat ret(nobs, K);
  for (unsigned int i = 0; i < K; ++i) {
    ret.col(i) = sim(arma::span(i*totNobs + idLoc[id], i * totNobs + idLoc[id+1]-1));
  }
  return ret;
}

typedef struct {
  arma::mat matsim;
  arma::mat epredt;
  arma::mat epred;
  arma::mat varsim;
  arma::mat ymat;
  arma::mat ydsim;
  arma::mat yobst;
  arma::mat yobs;
  arma::mat ydobs;
  arma::mat tcomp;
  arma::mat pd;
  arma::mat npde;
  arma::mat eres;
  unsigned int warn = 0;
} calcNpdeInfoId;

arma::mat decorrelateNpdeMat(arma::mat& varsim, unsigned int& warn, unsigned int &id, double &tolChol) {
  arma::mat ch, vYi;
  try {
    ch = chol(varsim);
  } catch (...){
    // Try nearPD from R/RxODE
    if (!(warn & 1)) {
      warn += 1;
    }
    try {
      ch = cholSE__(vYi, tolChol);
    } catch (...){
      Rcpp::stop("The simulated data covariance matrix Cholesky decomposition could not be obtained for %d.", id+1);
    }
  }
  try{
    vYi = trans(inv(trimatu(ch)));
  } catch (...) {
    if (!(warn & 1)) {
      warn += 2;
    }
    try {
      vYi = trans(pinv(trimatu(ch)));
    } catch (...){
      stop("The simulated data covariance matrix Cholesky decomposition inverse or pseudo-inverse could not be obtained for %d.", id+1);
    }
  }
  return vYi;
}



// This is similar to a truncated normal BUT the truncated normal handles the range (low,hi)
// so instead of updating the DV based on cdf method, simply use the truncated normal
// we also don't need to back-calculate the simulated DV value
static inline void handleCensNpdeCdf(calcNpdeInfoId &ret, arma::ivec &cens, unsigned int i, arma::vec &ru2,  arma::vec &ru3, unsigned int& K, bool &ties) {
  // For left censoring NPDE the probability is already calculated with the current pd
  // 1. Replace value with lloq
  // 2. The current pd represents the probability being below the limit of quantitation;
  //    If it is right censored 1-pde[i] represents the probability of being above the limit of quantitation
  // 3. For now `limit` is ignored
  // 4. For blq the pd is replaced with runif(0, p(bloq) or p(paloq))
  // 5. The epred is then replaced with back-calculated uniform value based on sorted tcomp of row
  arma::vec curRow;
  unsigned int j;
  double low, hi;
  switch (cens[i]) {
  case 1:
    ret.pd[i] = ru2[i]*ret.pd[i];
    break;
  case -1:
    ret.pd[i] = 1-ru2[i]*ret.pd[i];
    break;
  default:
    return;
  }
  curRow = sort(ret.matsim.row(i));
  j = trunc(ret.pd[i]*K);
  low = curRow[j];
  if (j+1 == K) hi = 2*low - curRow[j-1];
  else hi = curRow[j+1];
  if (ties) ret.yobst[i] = hi;
  else ret.yobst[i] = low + (hi - low)*ru3[i];
}

calcNpdeInfoId calcNpdeId(arma::ivec& idLoc, arma::vec &sim,
			  arma::vec &dvt, arma::ivec &censIn, unsigned int& id,
			  unsigned int& K, double &tolChol, bool &ties,
			  arma::vec &ruIn, arma::vec &ru2In, arma::vec &ru3In,
			  arma::vec &lambda, arma::vec &yj, arma::vec &hi, arma::vec &low) {
  calcNpdeInfoId ret;
  ret.matsim = getSimMatById(idLoc, sim, id, K);
  // transformed y observations
  ret.yobst = dvt(span(idLoc[id], idLoc[id+1]-1));
  arma::vec ru = ruIn(span(idLoc[id], idLoc[id+1]-1));
  arma::vec ru2 = ru2In(span(idLoc[id], idLoc[id+1]-1));
  arma::vec ru3 = ru3In(span(idLoc[id], idLoc[id+1]-1));
  arma::ivec cens = censIn(span(idLoc[id], idLoc[id+1]-1));
  ret.epredt = mean(ret.matsim, 1);
  ret.varsim = cov(trans(ret.matsim));
  ret.ymat = decorrelateNpdeMat(ret.varsim, ret.warn, id, tolChol);
  arma::mat ymatt = trans(ret.ymat);
  ret.ydsim = ret.matsim;
  ret.ydsim.each_col() -= ret.epredt;
  ret.ydsim = ymatt * ret.ydsim;
  ret.ydobs = ymatt * (ret.yobst - ret.epredt);
  // sim < obs
  ret.tcomp = ret.ydsim;
  for (unsigned int j = K; j--;) {
    for (unsigned int i = ret.ydsim.n_rows; i--;) {
      ret.tcomp(i, j) = ret.tcomp(i, j) < ret.ydobs[i];
    }
  }
  ret.pd = mean(ret.tcomp, 1);
  ret.npde = arma::mat(ret.pd.n_rows, 1);
  if (ties){
    // Ties are allowed
    for (unsigned int j = ret.pd.n_rows; j--;) {
      // handleCensNpdeCdf(ret, cens, j, ru2, ru3, K, ties);
      if (fabs(ret.pd[j]) < DOUBLE_EPS){
	ret.pd[j] = 1 / (2.0 * K);
      } else if (fabs(1-ret.pd[j]) < DOUBLE_EPS){
	ret.pd[j] = 1 - 1 / (2.0 * K);
      }
      ret.npde[j] = Rf_qnorm5(ret.pd[j], 0.0, 1.0, 1, 0);
    }
  } else {
    // Ties are discouraged with jitter
    for (unsigned int j = ret.pd.n_rows; j--;) {
      // handleCensNpdeCdf(ret, cens, j, ru2, ru3, K, ties);
      if (fabs(ret.pd[j]) < DOUBLE_EPS){
	ret.pd[j] = ru[j] / K;
      } else if (fabs(1-ret.pd[j]) < DOUBLE_EPS){
	ret.pd[j] = 1.0 - ru[j] / K;
      } else  {
	ret.pd[j] += ru[j]/K;
      }
      ret.npde[j] = Rf_qnorm5(ret.pd[j], 0.0, 1.0, 1, 0);
    }
  }
  ret.epred = arma::vec(ret.yobst.size());
  ret.yobs = arma::vec(ret.yobst.size());
  ret.eres = arma::vec(ret.yobst.size());
  for (unsigned int j = ret.yobst.size(); j--; ) {
    // Transfer back to original scale ie exp(x) -> log(exp(x))
    ret.yobs[j] = _powerD(ret.yobst[j], lambda[j], (int) yj[j], low[j], hi[j]);
    ret.epred[j] = _powerD(ret.epredt[j], lambda[j], (int) yj[j], low[j], hi[j]);
    ret.eres[j] = ret.yobs[j] - ret.epred[j];
  }
  return ret;
}

extern "C" SEXP _nlmixr_npdeCalc(SEXP npdeSim, SEXP dvIn, SEXP evidIn, SEXP censIn, SEXP limitIn, SEXP npdeOpt) {
  if (TYPEOF(npdeSim) != VECSXP) {
    Rf_errorcall(R_NilValue, "npdeSim needs to be a data.frame");
  }
  int dvLen = Rf_length(dvIn);
  arma::vec dv  = arma::vec(REAL(dvIn), dvLen, false, true);
  //arma::vec npde(REAL(npdeSEXP), dv.size(), false, true);
  int pro = 0;
  SEXP s0 = PROTECT(VECTOR_ELT(npdeSim, 0)); pro++;
  int simLen = Rf_length(s0);
  arma::ivec aSimIdVec(INTEGER(s0), simLen, false, true);
  arma::ivec aIdVec(INTEGER(VECTOR_ELT(npdeSim, 1)), simLen, false, true);
  unsigned int nid, K;
  arma::ivec idLoc = getSimIdLoc(aIdVec, aSimIdVec, nid, K);
  arma::vec sim(REAL(VECTOR_ELT(npdeSim, 3)), simLen, false, true);
  arma::vec lambda(REAL(VECTOR_ELT(npdeSim, 4)), dvLen, false, true);
  arma::vec yj(REAL(VECTOR_ELT(npdeSim, 5)), dvLen, false, true);
  arma::vec low(REAL(VECTOR_ELT(npdeSim, 6)), dvLen, false, true);
  arma::vec hi(REAL(VECTOR_ELT(npdeSim, 7)), dvLen, false, true);
  arma::vec dvt(dvLen);
  // dv -> dv transform
  for (unsigned int i = dvLen; i--; ) {
    // powerDi for log-normal transfers dv = exp(dv)
    dvt[i] = _powerDi(dv[i], lambda[i], (int) yj[i], low[i], hi[i]);
  }
  arma::ivec cens;
  if (Rf_isNull(censIn)) {
    cens = arma::ivec(dvLen, fill::zeros);
  } else {
    cens = as<arma::ivec>(censIn);
  }
  arma::ivec evid;
  if (Rf_isNull(evidIn)) {
    evid = arma::ivec(dvLen, fill::zeros);
  } else {
    evid = as<arma::ivec>(evidIn);
  }
  bool doLimit = false;
  arma::vec ru = randu(dvLen); // Pre-fill uniform random numbers to make sure independent
  arma::vec ru2 = randu(dvLen);
  arma::vec ru3 = randu(dvLen);
  double tolChol = 6.055454e-06;
  bool ties = false;

  SEXP npdeSEXP = PROTECT(Rf_allocVector(REALSXP, dvLen)); pro++;
  SEXP epredSEXP = PROTECT(Rf_allocVector(REALSXP, dvLen)); pro++;
  SEXP dvSEXP = PROTECT(Rf_allocVector(REALSXP, dvLen)); pro++;
  SEXP eresSEXP = PROTECT(Rf_allocVector(REALSXP, dvLen)); pro++;
  arma::vec npde(REAL(npdeSEXP), dvLen, false, true);
  arma::vec epred(REAL(epredSEXP), dvLen, false, true);
  arma::vec dvf(REAL(dvSEXP), dvLen, false, true);
  arma::vec eres(REAL(eresSEXP), dvLen, false, true);
  for (unsigned int curid = 0; curid < idLoc.size()-1; ++curid) {
    calcNpdeInfoId idInfo = calcNpdeId(idLoc, sim, dvt, cens, curid, K, tolChol, ties, ru, ru2, ru3,
				       lambda, yj, hi, low);
    npde(span(idLoc[curid],idLoc[curid+1]-1)) = idInfo.npde;
    epred(span(idLoc[curid], idLoc[curid+1]-1)) = idInfo.epred;
    dvf(span(idLoc[curid], idLoc[curid+1]-1)) = idInfo.yobs;
    eres(span(idLoc[curid], idLoc[curid+1]-1)) = idInfo.eres;
  }
  List ret(4);
  // epred, eres, npde, dv
  ret[0] = epred;
  ret[1] = eres;
  ret[2] = npde;
  ret[3] = dvf;
  Rf_setAttrib(ret, R_ClassSymbol, wrap("data.frame"));
  Rf_setAttrib(ret, R_RowNamesSymbol,
	       IntegerVector::create(NA_INTEGER, -dvLen));
  Rf_setAttrib(ret, R_NamesSymbol, CharacterVector::create("EPRED", "ERES", "NPDE",  "DV"));
  UNPROTECT(pro);
  return ret;
}
