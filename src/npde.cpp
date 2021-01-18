#include "npde.h"

#ifdef _OPENMP
#include <omp.h>
#endif

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
      Rcpp::stop("The simulated data covariance matrix Cholesky decomposition could not be obtained for  id=%d.", id+1);
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
      stop("The simulated data covariance matrix Cholesky decomposition inverse or pseudo-inverse could not be obtained for id=%d.", id+1);
    }
  }
  return vYi;
}

// This is similar to a truncated normal BUT the truncated normal handles the range (low,hi)
// so instead of updating the DV based on cdf method, simply use the truncated normal
// we also don't need to back-calculate the simulated DV value
static inline void handleCensNpdeCdf(calcNpdeInfoId &ret, arma::ivec &cens, arma::vec &limit, int &censMethod, bool &doLimit,
				     unsigned int i, arma::vec &ru2,  arma::vec &ru3, unsigned int& K, bool &ties) {
  if (censMethod != CENS_CDF) return;
  // For left censoring NPDE the probability is already calculated with the current pd
  // 1. Replace value with lloq
  // 2. The current pd represents the probability being below the limit of quantitation;
  //    If it is right censored 1-pde[i] represents the probability of being above the limit of quantitation
  // 4. For blq the pd is replaced with runif(0, p(bloq)) or runif(p(paloq), 1)
  // 5. The epred is then replaced with back-calculated uniform value based on sorted tcomp of row
  arma::vec curRow;
  unsigned int j;
  unsigned int nlimit = 0;
  double low, hi;
  switch (cens[i]) {
  case 1:
    // Sort up here to reduce the calculation burden for the p(lloq)
    curRow = sort(ret.matsim.row(i));
    if (doLimit && R_FINITE(limit[i])) {
      //  in this case (limit, dv) simulate pd between these two probabilities
      while (nlimit < curRow.size() && curRow[nlimit] < limit[i]) nlimit++;
      low = nlimit/curRow.size();
      // Simulate between (limit, dv)
      ret.pd[i] = low +  ru2[i]*(ret.pd[i]-low);
    } else {
      ret.pd[i] = ru2[i]*ret.pd[i];
    }
    break;
  case -1:
    curRow = sort(ret.matsim.row(i));
    if (doLimit && R_FINITE(limit[i])) {
      // In this case (dv, limit) simulate pd between these two probabilities
      while (nlimit < curRow.size() && curRow[curRow.size()-nlimit-1] > limit[i]) nlimit++;
      hi = nlimit/curRow.size();
      low =  1-ret.pd[i];
      ret.pd = low + ru2[i]*(hi-low);
    } else {
      ret.pd[i] = 1-ru2[i]*ret.pd[i];
    }
    break;
  default:
    return;
  }
  // Now back-calculate the EPRED
  j = trunc(ret.pd[i]*K);
  low = curRow[j];
  if (j+1 == K) hi = 2*low - curRow[j-1];
  else hi = curRow[j+1];
  if (ties) ret.yobst[i] = hi;
  else ret.yobst[i] = low + (hi - low)*ru3[i];
}

static inline void handleNpdeNAandCalculateEpred(calcNpdeInfoId& ret, unsigned int& K) {
  ret.namat = umat(ret.matsim.n_rows, ret.matsim.n_cols);
  // Replace NA with zeros
  for (unsigned int j = ret.namat.size(); j--;) {
    ret.namat[j] = ISNA(ret.matsim[j]);
    if (ret.namat[j]) ret.matsim[j] = 0.0;
  }
  // mean = X/K * K/(K-sum(isNa))
  ret.epredt = mean(ret.matsim, 1) % (K/(K-mean(ret.namat, 1)));
  // Now replace NA with epredt for further calculations
  for (unsigned int i =ret.namat.n_rows; i--;) {
    for (unsigned int j = ret.namat.n_cols; j--;) {
      if (ret.namat(i,j)) {
	ret.matsim(i,j) = ret.epredt[i];
      }
    }
  }
}

static inline void calculatePD(calcNpdeInfoId& ret, unsigned int& id, unsigned int &K, double &tolChol) {
  ret.ydsim = ret.matsim.rows(ret.obs);
  ret.varsim = cov(trans(ret.ydsim));
  ret.ymat = decorrelateNpdeMat(ret.varsim, ret.warn, id, tolChol);
  arma::mat ymatt = trans(ret.ymat);
  ret.ydsim.each_col() -= ret.epredt.elem(ret.obs);
  ret.ydsim = ymatt * ret.ydsim;
  ret.ydobs = ymatt * (ret.yobst.elem(ret.obs) - ret.epredt.elem(ret.obs));
  // sim < obs
  ret.tcomp = ret.ydsim;
  for (unsigned int j = K; j--;) {
    for (unsigned int i = ret.ydsim.n_rows; i--;) {
      ret.tcomp(i, j) = ret.tcomp(i, j) < ret.ydobs[i];
    }
  }
  arma::mat pdObs = mean(ret.tcomp, 1);
  ret.pd = arma::mat(ret.matsim.n_rows, 1, fill::zeros);
  ret.pd.rows(ret.obs) = pdObs;
}

static inline void calculateNPDEfromPD(calcNpdeInfoId &ret, arma::ivec &cens, arma::vec &limit, int &censMethod, bool &doLimit,
				       unsigned int &K, bool &ties, arma::vec &ru, arma::vec &ru2, arma::vec& ru3) {
  ret.npde = arma::mat(ret.pd.n_rows, 1);
  if (ties){
    // Ties are allowed
    for (unsigned int j = ret.pd.n_rows; j--;) {
      handleCensNpdeCdf(ret, cens, limit, censMethod, doLimit, j, ru2, ru3, K, ties);
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
      handleCensNpdeCdf(ret, cens, limit, censMethod, doLimit, j, ru2, ru3, K, ties);
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
}

calcNpdeInfoId calcNpdeId(arma::ivec& idLoc, arma::vec &sim,
			  arma::vec &dvt, arma::ivec &evidIn, arma::ivec &censIn, arma::vec &limitIn,
			  int &censMethod, bool &doLimit,
			  unsigned int& id,
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
  arma::vec limit = limitIn(span(idLoc[id], idLoc[id+1]-1));
  arma::ivec evid = evidIn(span(idLoc[id], idLoc[id+1]-1));
  handleNpdeNAandCalculateEpred(ret, K);
  ret.obs = find(evid == 0);
  calculatePD(ret, id, K, tolChol);
  calculateNPDEfromPD(ret, cens, limit, censMethod, doLimit, K, ties, ru, ru2, ru3);
  ret.epred = arma::vec(ret.yobst.size());
  ret.yobs = arma::vec(ret.yobst.size());
  ret.eres = arma::vec(ret.yobst.size());
  for (unsigned int j = ret.yobst.size(); j--; ) {
    // Transfer back to original scale ie log(x) -> exp(log(x))
    ret.epred[j] = _powerDi(ret.epredt[j], lambda[j], (int) yj[j], low[j], hi[j]);
    if (censMethod == CENS_EPRED && cens[j] != 0) {
      ret.yobst[j] = ret.epredt[j];
    } 
    ret.yobs[j] = _powerDi(ret.yobst[j], lambda[j], (int) yj[j], low[j], hi[j]);
    ret.eres[j] = ret.yobs[j] - ret.epred[j];
    if (censMethod == CENS_OMIT && cens[j] != 0) {
      ret.yobs[j] = NA_REAL;
      ret.epred[j] = NA_REAL;
      ret.eres[j] = NA_REAL;
    } else if (evid[j] != 0) {
      ret.yobs[j] = NA_REAL;
      ret.eres[j] = NA_REAL;
      ret.npde[j] = NA_REAL;
    }
  }
  return ret;
}

extern "C" SEXP _nlmixr_npdeCalc(SEXP npdeSim, SEXP dvIn, SEXP evidIn, SEXP censIn, SEXP limitIn, SEXP npdeOpt) {
  BEGIN_RCPP
  if (TYPEOF(npdeSim) != VECSXP) {
    Rf_errorcall(R_NilValue, "npdeSim needs to be a data.frame");
  }
  List opt = as<List>(npdeOpt);
  double tolChol = 6.055454e-06;
  if (opt.containsElementNamed("tolChol")) {
    RObject tmp = opt["tolChol"];
    if (TYPEOF(tmp) == REALSXP) {
      tolChol = as<double>(tmp);
    }
  }
  bool ties = false;
  if (opt.containsElementNamed("ties")) {
    RObject tmp = opt["ties"];
    if (TYPEOF(tmp) == LGLSXP) {
      ties = as<bool>(tmp);
    }
  }
  int censMethod = CENS_TNORM;
  if (opt.containsElementNamed("censMethod")) {
    RObject tmp = opt["censMethod"];
    if (TYPEOF(tmp) == INTSXP) {
       censMethod = as<unsigned int>(tmp);
    }
  }

  List npdeSimL = as<List>(npdeSim);
  int nsim = getPredIndex(npdeSimL);

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
  arma::vec sim(REAL(VECTOR_ELT(npdeSim, nsim)), simLen, false, true);
  arma::vec lambda(REAL(VECTOR_ELT(npdeSim, nsim+1)), simLen, false, true);
  arma::vec yj(REAL(VECTOR_ELT(npdeSim, nsim+2)), simLen, false, true);
  arma::vec low(REAL(VECTOR_ELT(npdeSim, nsim+3)), simLen, false, true);
  arma::vec hi(REAL(VECTOR_ELT(npdeSim, nsim+4)), simLen, false, true);
  arma::vec dvt(dvLen);
  // dv -> dv transform
  for (unsigned int i = dvLen; i--; ) {
    // powerDi for log-normal transfers dv = log(dv)
    dvt[i] = _powerD(dv[i], lambda[i], (int) yj[i], low[i], hi[i]);
  }
  // sim -> sim transform; simulation is on untransformed scale
  for (unsigned int i = simLen; i--;) {
    sim[i] = _powerD(sim[i], lambda[i], (int) yj[i], low[i], hi[i]);
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

  arma::vec limit;
  int hasLimit=0;
  getLimitFromInput(limitIn, dvLen, limit, hasLimit);

  arma::vec ru = randu(dvLen); // Pre-fill uniform random numbers to make sure independent
  arma::vec ru2 = randu(dvLen);
  arma::vec ru3 = randu(dvLen);

  SEXP npdeSEXP = PROTECT(Rf_allocVector(REALSXP, dvLen)); pro++;
  SEXP epredSEXP = PROTECT(Rf_allocVector(REALSXP, dvLen)); pro++;
  SEXP dvSEXP = PROTECT(Rf_allocVector(REALSXP, dvLen)); pro++;
  SEXP eresSEXP = PROTECT(Rf_allocVector(REALSXP, dvLen)); pro++;
  arma::vec npde(REAL(npdeSEXP), dvLen, false, true);
  arma::vec epred(REAL(epredSEXP), dvLen, false, true);
  arma::vec dvf(REAL(dvSEXP), dvLen, false, true);
  arma::vec eres(REAL(eresSEXP), dvLen, false, true);
  
#pragma omp parallel for
  for (unsigned int curid = 0; curid < idLoc.size()-1; ++curid) {
    calcNpdeInfoId idInfo = calcNpdeId(idLoc, sim, dvt, evid, cens, limit, censMethod, doLimit, curid, K, tolChol, ties, ru, ru2, ru3,
				       lambda, yj, hi, low);
    npde(span(idLoc[curid],idLoc[curid+1]-1)) = idInfo.npde;
    epred(span(idLoc[curid], idLoc[curid+1]-1)) = idInfo.epred;
    dvf(span(idLoc[curid], idLoc[curid+1]-1)) = idInfo.yobs;
    eres(span(idLoc[curid], idLoc[curid+1]-1)) = idInfo.eres;
  }
  List ret(3);
  // epred, eres, npde, dv
  ret[0] = epred;
  ret[1] = eres;
  ret[2] = npde;
  Rf_setAttrib(ret, R_ClassSymbol, wrap("data.frame"));
  Rf_setAttrib(ret, R_RowNamesSymbol,
	       IntegerVector::create(NA_INTEGER, -dvLen));
  Rf_setAttrib(ret, R_NamesSymbol, CharacterVector::create("EPRED", "ERES", "NPDE"));
  UNPROTECT(pro);
  return List::create(List::create(_["DV"]=dv),
		      ret);
  END_RCPP
}
