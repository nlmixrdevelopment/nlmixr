#define STRICT_R_HEADER
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

arma::mat decorrelateNpdeEigenMat(arma::mat& varsim, unsigned int& warn) {
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, varsim, "std");
  eigval = sqrt(eigval);
  arma::mat iEigVec;
  try{
    iEigVec = inv(eigvec);
    warn = NPDE_DECORRELATE_EIGEN;
  } catch(...) {
    iEigVec = pinv(eigvec);
    warn = NPDE_DECORRELATE_EIGEN_PINV;
  }
  arma::mat ret =  eigvec * diagmat(eigval) * iEigVec;
  // Try nearPD from R/RxODE
  return ret;
}

arma::mat varNpdMat(arma::mat& varsim) {
  // No decorrelation step
  arma::mat ret(varsim.n_rows, varsim.n_cols, fill::zeros);
  for (unsigned int i = varsim.n_rows; i--;) {
    ret(i,i) = 1/sqrt(varsim(i,i));
  }
  return ret;
}

arma::mat decorrelateNpdeMat(arma::mat& varsim, unsigned int& warn, unsigned int &id, double &tolChol) {
  arma::mat ch, vYi;
  try {
    ch = chol(varsim);
  } catch (...){
    try {
      return decorrelateNpdeEigenMat(varsim, warn);
    } catch (...){
      try {
	ch = cholSE__(varsim, tolChol);
	warn = NPDE_CHOLSE;
      } catch (...) {
	warn = NPDE_NPD;
	return varNpdMat(varsim);
      }
    }
  }
  try{
    vYi = trans(inv(trimatu(ch)));
    if (warn != NPDE_CHOLSE) {
      warn = NPDE_CHOL;
    }
  } catch (...) {
    try {
      vYi = trans(pinv(trimatu(ch)));
      if (warn != NPDE_CHOLSE) {
	warn = NPDE_CHOL_PINV;
      } else {
	warn = NPDE_CHOLSE_PINV;
      }
    } catch (...){
      warn = NPDE_NPD;
      return varNpdMat(varsim);
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
  unsigned int j, j2;
  double low, hi, low2, hi2;
  switch (cens[i]) {
  case 1:
    curRow = sort(trans(ret.matsim.row(i)));
    ret.pd[i] = ru2[i]*ret.pd[i];
    ret.pd2[i] = ru2[i]*ret.pd2[i];
    break;
  case -1:
    curRow = sort(trans(ret.matsim.row(i)));
    ret.pd[i] = 1-ru2[i]*ret.pd[i];
    ret.pd2[i] = 1-ru2[i]*ret.pd2[i];
    break;
  default:
    return;
  }
  // Now back-calculate the EPRED
  j = trunc(ret.pd[i]*K);
  j2 = trunc(ret.pd2[i]*K);
  low = curRow[j];
  low2 = curRow[j2];
  if (j+1 == K) hi = 2*low - curRow[j-1];
  else hi = curRow[j+1];
  if (j2+1 == K) hi2 = 2*low2 - curRow[j2-1];
  else hi2 = curRow[j2+1];
  // check if this makes sense
  // Use npd instead of npde as suggested by Nguyen2017
  // Could be turned on/off by npdeControl()
  if (ties) {
    ret.yobst[i] = hi2;
  }
  else ret.yobst[i] = low2 + (hi2 - low2)*ru3[i];
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
  ret.ymat2 = varNpdMat(ret.varsim);
  arma::mat ymatt = trans(ret.ymat);
  ret.ydsim.each_col() -= ret.epredt.elem(ret.obs);
  ret.ydsim2 = ret.ydsim;
  ret.ydsim = ymatt * ret.ydsim;
  ret.ydsim2 = ret.ymat2 * ret.ydsim2; 
  ret.ydobs = ymatt * (ret.yobst.elem(ret.obs) - ret.epredt.elem(ret.obs));
  ret.ydobs2 = ret.ymat2 * (ret.yobst.elem(ret.obs) - ret.epredt.elem(ret.obs));
  // sim < obs
  ret.tcomp = ret.ydsim;
  ret.tcomp2 = ret.ydsim2;
  for (unsigned int j = K; j--;) {
    for (unsigned int i = ret.ydsim.n_rows; i--;) {
      ret.tcomp(i, j) = ret.tcomp(i, j) < ret.ydobs[i];
      ret.tcomp2(i, j) = ret.tcomp2(i, j) < ret.ydobs2[i];
    }
  }
  arma::mat pdObs = mean(ret.tcomp, 1);
  arma::mat pdObs2 = mean(ret.tcomp2, 1);
  ret.pd = arma::mat(ret.matsim.n_rows, 1, fill::zeros);
  ret.pd2 = arma::mat(ret.matsim.n_rows, 1, fill::zeros);
  ret.pd.rows(ret.obs) = pdObs;
  ret.pd2.rows(ret.obs) = pdObs;
}

static inline void calculateNPDEfromPD(calcNpdeInfoId &ret, arma::ivec &cens, arma::vec &limit, int &censMethod, bool &doLimit,
				       unsigned int &K, bool &ties, arma::vec &ru, arma::vec &ru2, arma::vec& ru3) {
  ret.npde = arma::mat(ret.pd.n_rows, 1);
  ret.npd = arma::mat(ret.pd.n_rows, 1);
  if (ties){
    // Ties are allowed
    for (unsigned int j = ret.pd.n_rows; j--;) {
      handleCensNpdeCdf(ret, cens, limit, censMethod, doLimit, j, ru2, ru3, K, ties);
      if (fabs(ret.pd[j]) < DBL_EPSILON){
	ret.pd[j] = 1 / (2.0 * K);
      } else if (fabs(1-ret.pd[j]) < DBL_EPSILON){
	ret.pd[j] = 1 - 1 / (2.0 * K);
      }
      ret.npde[j] = Rf_qnorm5(ret.pd[j], 0.0, 1.0, 1, 0);
      if (fabs(ret.pd2[j]) < DBL_EPSILON){
	ret.pd2[j] = 1 / (2.0 * K);
      } else if (fabs(1-ret.pd2[j]) < DBL_EPSILON){
	ret.pd2[j] = 1 - 1 / (2.0 * K);
      }
      ret.npd[j] = Rf_qnorm5(ret.pd2[j], 0.0, 1.0, 1, 0);
    }
  } else {
    // Ties are discouraged with jitter
    for (unsigned int j = ret.pd.n_rows; j--;) {
      handleCensNpdeCdf(ret, cens, limit, censMethod, doLimit, j, ru2, ru3, K, ties);
      if (fabs(ret.pd[j]) < DBL_EPSILON){
	ret.pd[j] = ru[j] / K;
      } else if (fabs(1-ret.pd[j]) < DBL_EPSILON){
	ret.pd[j] = 1.0 - ru[j] / K;
      } else  {
	ret.pd[j] += ru[j]/K;
      }
      ret.npde[j] = Rf_qnorm5(ret.pd[j], 0.0, 1.0, 1, 0);
      if (fabs(ret.pd2[j]) < DBL_EPSILON){
	ret.pd2[j] = ru[j] / K;
      } else if (fabs(1-ret.pd2[j]) < DBL_EPSILON){
	ret.pd2[j] = 1.0 - ru[j] / K;
      } else  {
	ret.pd2[j] += ru2[j]/K;
      }
      ret.npd[j] = Rf_qnorm5(ret.pd2[j], 0.0, 1.0, 1, 0);
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
      ret.npd[j] = NA_REAL;
    }
  }
  return ret;
}

rxGetId2_t rxGetId2;

extern "C" SEXP _nlmixr_npdeCalc(SEXP npdeSim, SEXP dvIn, SEXP evidIn, SEXP censIn, SEXP limitIn, SEXP npdeOpt) {
  BEGIN_RCPP
  rxGetId2 = (rxGetId2_t) R_GetCCallable("RxODE", "rxGetId");
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

  if (hasLimit && censMethod == CENS_CDF) {
    for (unsigned int i = 0; i < cens.size(); ++i) {
      switch(cens[i]){
      case 1:
	if (R_FINITE(limit[i])) {
	  warning("limits are ignored for npde back-transformation with 'cdf' method");
	  i= cens.size();
	}
	break;
      case -1:
	if (R_FINITE(limit[i])) {
	  warning("limits are ignored for npde back-transformation with 'cdf' method");
	  i= cens.size();
	}
	break;
      case 0:
	break;
      }
    }
  }

  arma::vec ru = randu(simLen); // Pre-fill uniform random numbers to make sure independent
  arma::vec ru2 = randu(simLen);
  arma::vec ru3 = randu(simLen);

  SEXP npdeSEXP = PROTECT(Rf_allocVector(REALSXP, dvLen)); pro++;
  SEXP npdSEXP = PROTECT(Rf_allocVector(REALSXP, dvLen)); pro++;
  SEXP epredSEXP = PROTECT(Rf_allocVector(REALSXP, dvLen)); pro++;
  SEXP dvSEXP = PROTECT(Rf_allocVector(REALSXP, dvLen)); pro++;
  SEXP eresSEXP = PROTECT(Rf_allocVector(REALSXP, dvLen)); pro++;
  arma::vec npde(REAL(npdeSEXP), dvLen, false, true);
  arma::vec npd(REAL(npdSEXP), dvLen, false, true);
  arma::vec epred(REAL(epredSEXP), dvLen, false, true);
  arma::vec dvf(REAL(dvSEXP), dvLen, false, true);
  arma::vec eres(REAL(eresSEXP), dvLen, false, true);
  arma::ivec warn(idLoc.size()-1);

  int cores = as<int>(opt["cores"]);
  
#pragma omp parallel for num_threads(cores)
  for (unsigned int curid = 0; curid < idLoc.size()-1; ++curid) {
    calcNpdeInfoId idInfo = calcNpdeId(idLoc, sim, dvt, evid, cens, limit, censMethod, doLimit, curid, K, tolChol, ties, ru, ru2, ru3,
				       lambda, yj, hi, low);
    npde(span(idLoc[curid],idLoc[curid+1]-1)) = idInfo.npde;
    npd(span(idLoc[curid], idLoc[curid+1]-1)) = idInfo.npd;
    epred(span(idLoc[curid], idLoc[curid+1]-1)) = idInfo.epred;
    dvf(span(idLoc[curid], idLoc[curid+1]-1)) = idInfo.yobs;
    eres(span(idLoc[curid], idLoc[curid+1]-1)) = idInfo.eres;
    warn[curid] = idInfo.warn;
  }
  std::string sCholPinv = "";
  int nCholPinv = 0;
  std::string sEigen = "";
  int nEigen = 0;
  std::string sEigenPinv = "";
  int nEigenPinv = 0;
  std::string sCholSE = "";
  int nCholSE = 0;
  std::string sCholSEPinv = "";
  int nCholSEPinv = 0;
  std::string sPD = "";
  int nPD = 0;
  for (unsigned int curid = 0; curid < warn.size(); ++curid) {
    switch(warn[curid]) {
    case NPDE_CHOL_PINV:
      if (sCholPinv == "") sCholPinv = rxGetId2(curid);
      else {
	sCholPinv += ", ";
	sCholPinv += rxGetId2(curid);
      }
      nCholPinv++;
      break;
    case NPDE_DECORRELATE_EIGEN:
      if (sEigen == "") sEigen = rxGetId2(curid);
      else {
	sEigen += ", ";
	sEigen += rxGetId2(curid);
      }
      nEigen++;
      break;
    case NPDE_DECORRELATE_EIGEN_PINV:
      if (sEigenPinv == "")  sEigenPinv =  rxGetId2(curid);
      else {
	sEigenPinv += ", ";
	sEigenPinv += rxGetId2(curid);
      }
      nEigenPinv++;
      break;
    case NPDE_CHOLSE:
      if (sCholSE == "") sCholSE = rxGetId2(curid);
      else {
	sCholSE += ", ";
	sCholSE += rxGetId2(curid);
      }
      nCholSE++;
      break;
    case NPDE_CHOLSE_PINV:
      if (sCholSEPinv == "") sCholSEPinv = rxGetId2(curid);
      else {
	sCholSEPinv += ", ";
	sCholSEPinv += rxGetId2(curid);
      }
      nCholSEPinv++;
      break;
    case NPDE_NPD:
      if (sPD == "") sPD = rxGetId2(curid);
      else {
	sPD += ", ";
	sPD += rxGetId2(curid);
      }
      nPD++;
      break;
    }
  }
  double rCholPinv = (double)nCholPinv / (double)warn.size(),
    rEigen = (double)nEigen / (double)warn.size(),
    rEigenPinv = (double)nEigenPinv / (double)warn.size(),
    rCholSE = (double)nCholSE / (double)warn.size(),
    rCholSEPinv = (double)nCholSEPinv / (double)warn.size(),
    rPD = (double)nPD / (double)warn.size();
  if (sCholPinv != "") {
    Rf_warningcall(R_NilValue, _("npde decorrelation used Cholesky pseudo-inverse for %.1f%%, id: %s"), rCholPinv*100, sCholPinv.c_str());
  }
  if (sEigen != "") {
    Rf_warningcall(R_NilValue, _("npde decorrelation used Eigen-values for %.1f%%, id: %s"), rEigen*100, sEigen.c_str());
  }
  if (sEigenPinv != "") {
    Rf_warningcall(R_NilValue, _("npde decorrelation used Eigen-value pseudo-inverse for %.1f%%, id: %s"), rEigenPinv*100, sEigenPinv.c_str());
  }
  if (sCholSE != "") {
    Rf_warningcall(R_NilValue, _("npde decorrelation used generalized Cholesky for %.1f%%, id: %s"), rCholSE*100, sCholSE.c_str());
  }
  if (sCholSEPinv != "") {
    Rf_warningcall(R_NilValue, _("npde decorrelation used generalized Cholesky pseudo inverse for %.1f%%, id: %s"), rCholSEPinv*100, sCholSEPinv.c_str());
  }
  if (sPD != "") {
    Rf_warningcall(R_NilValue, _("npde decorrelation failed (return normalized prediction discrepancies) for %.1f%% id: %s"), rPD*100, sPD.c_str());
  }
  
  List ret(4);
  // epred, eres, npde, dv
  ret[0] = List::create(_["EPRED"]=epred);
  ret[1] = List::create(_["ERES"]=eres);
  ret[2] = List::create(_["NPDE"]=npde);
  ret[3] = List::create(_["NPD"]=npd);
  SEXP ret2 = PROTECT(dfCbindList(wrap(ret))); pro++;
  UNPROTECT(pro);
  return List::create(dvf, ret2);
  END_RCPP
}
