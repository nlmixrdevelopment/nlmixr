#include "armahead.h"

arma::uvec getSimIdLoc(arma::uvec& id, arma::uvec& simId,
		      unsigned int &nid, unsigned int &K) {
  unsigned int i = 0;
  unsigned int lastId=id[0];
  unsigned int lastSimId = simId[0];
  unsigned int totNobs = 0;
  nid = 1;
  while (simId[i] == lastSimId) {
    if (lastId != id[i]) {
      nid++;
      lastId = id[i];
    }
    i++; totNobs++;
  }
  K =   id.size() /  totNobs;

  arma::uvec idLoc(nid+1);
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

arma::mat getSimMatById(arma::uvec& idLoc, arma::vec &sim, unsigned int& id,
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
  arma::mat varsim;
  arma::mat ymat;
  arma::mat ydsim;
  arma::mat yobst;
  arma::mat ydobs;
  arma::mat tcomp;
  arma::mat pd;
  arma::mat npde;
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

calcNpdeInfoId calcNpdeId(arma::uvec& idLoc, arma::vec &sim,
			  arma::vec &dvt,
			  unsigned int& id,
			  unsigned int& K, double &tolChol,
			  bool &ties, arma::vec &ruIn) {
  calcNpdeInfoId ret;
  ret.matsim = getSimMatById(idLoc, sim, id, K);
  // transformed y observations
  ret.yobst = dvt(span(idLoc[id], idLoc[id+1]-1));
  arma::vec ru = ruIn(span(idLoc[id], idLoc[id+1]-1));
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
      if (fabs(ret.pd[j]) < DOUBLE_EPS){
	pd[j] = 1 / (2.0 * K);
      } else if (fabs(1-ret.pd[i]) < DOUBLE_EPS){
	pd[j] = 1 - 1 / (2.0 * K);
      }
      ret.npde[j] = Rf_qnorm5(ret.npd[j], 0.0, 1.0, 1, 0);
    }
  } else {
    // Ties are discouraged with jitter
    for (unsigned int j = ret.pd.n_rows; j--;) {
      if (fabs(ret.pd[j]) < DOUBLE_EPS){
	pd[j] = ru[j] / K;
      } else if (fabs(1-ret.pd[i]) < DOUBLE_EPS){
	pd[j] = 1.0 - ru[j] / K;
      } else  {
	pj[j] += ru[j]/K;
      }
      ret.npde[j] = Rf_qnorm5(ret.npd[j], 0.0, 1.0, 1, 0);
    }
  }
  //
  return ret;
}

extern "C" SEXP _nlmixr_npdeCalc(SEXP id, SEXP simId, SEXP idVec, SEXP simIn, SEXP dvIn) {
  unsigned int curid = INTEGER(id)[0];
  arma::uvec aIdVec = as<arma::uvec>(idVec);
  arma::uvec aSimIdVec = as<arma::uvec>(simId);

  unsigned int nid, K;
  arma::uvec idLoc = getSimIdLoc(aIdVec, aSimIdVec, nid, K);
  arma::vec sim = as<arma::vec>(simIn);
  arma::vec dv  = as<arma::vec>(dvIn);
  arma::vec ru = randu(dv.size()); // Pre-fill uniform random numbers to make sure independent

  double tolChol = 6.055454e-06;

  calcNpdeInfoId idInfo = calcNpdeId(idLoc, sim, dv, curid, K, tolChol);
  
  List ret(7);
  ret[0] = wrap(idInfo.epredt);
  ret[1] = wrap(idInfo.ymat);
  ret[2] = wrap(idInfo.ydsim);
  ret[3] = wrap(idInfo.ydobs);
  ret[4] = wrap(idInfo.tcomp);
  ret[5] = wrap(idInfo.pde);
  ret[6] = wrap(idInfo.npde);
  return ret;
}
