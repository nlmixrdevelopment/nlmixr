#include <stdarg.h>
#define ARMA_DONT_PRINT_ERRORS
#define ARMA_DONT_USE_OPENMP // Known to cause speed problems
#include <RcppArmadillo.h>
#include <R.h>
#include <Rmath.h>
#include "nlmixr_types.h"
#include <RxODE.h>

typedef double (*rxPow)(double x, double lambda, int yj);

double powerDi(double x, double lambda, int yj){
  static rxPow fun=NULL;
  if (fun == NULL) fun = (rxPow) R_GetCCallable("RxODE","powerDi");
  return fun(x, lambda, yj);
}

double powerD(double x, double lambda, int yj){
  static rxPow fun=NULL;
  if (fun == NULL) fun = (rxPow) R_GetCCallable("RxODE","powerD");
  return fun(x, lambda, yj);
}

using namespace Rcpp;
using namespace R;
using namespace arma;


//[[Rcpp::export]]
List nlmixrParameters(NumericVector theta, DataFrame eta){
  DataFrame eta2 = eta;
  eta2.erase(0);
  NumericVector pred(theta.size()+eta2.ncol());
  CharacterVector predn(pred.size());
  NumericMatrix ipred(eta.nrow(),pred.size());
  unsigned int j;
  for (j=theta.size(); j--;){
    pred[j] = theta[j];
    std::fill_n(ipred.begin()+eta.nrow()*j,eta.nrow(),theta[j]);
    predn[j] = "THETA[" + std::to_string(j+1) + "]";
  }
  unsigned int k = theta.size();
  unsigned int n, n1;
  double M1, M2, M3, M4, delta, delta_n, delta_n2, term1;
  // Create data frame for ETA/EPS mean, sd, variance, kurtosis and skewness
  List etadf(eta2.ncol()+1);
  NumericVector cur1(9);
  etadf[eta2.ncol()]=cur1;
  for (j=eta2.ncol();j--;){
    pred[k+j]=0;
    predn[k+j] = "ETA[" + std::to_string(j+1) + "]";
    NumericVector cur = as<NumericVector>(eta2[j]);
    std::copy(cur.begin(),cur.end(),ipred.begin()+(j+k)*eta.nrow());
    // This uses welford's method as described by
    // https://www.johndcook.com/blog/skewness_kurtosis/ adapted to
    // this problem.
    n =0; n1 = 0;
    M1 =  M2 = M3 = M4 =0;
    NumericVector stat(9);
    for (unsigned int i = eta.nrow(); i--;){
      n1=n; n++;
      delta = cur[i] - M1;
      delta_n = delta / n;
      delta_n2 = delta_n * delta_n;
      term1 = delta * delta_n * n1;
      M1 += delta_n;
      M4 += term1 * delta_n2 * (n*n - 3*n + 3) + 6 * delta_n2 * M2 - 4 * delta_n * M3;
      M3 += term1 * delta_n * (n - 2) - 3 * delta_n * M2;
      M2 += term1;
    }
    stat[0] = M1;
    stat[1] = M2/((double)n1);
    stat[2] = sqrt((double)stat[1]);
    stat[3] = sqrt((double)(n)) * M3/ pow(M2, 1.5);
    stat[4] = (double)(n)*M4 / (M2*M2) - 3.0;
    etadf[j] = stat;
  }
  pred.names() = predn;
  ipred.attr("dimnames") = List::create(R_NilValue,predn);
  return List::create(_["pred"]=pred, _["ipred"]=ipred,_["eta.lst"]=etadf);
}

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

// [[Rcpp::export]]
List nlmixrResid(List &innerList, NumericMatrix &omegaMat, NumericVector &cdv,
		 IntegerVector &evid, NumericVector &lambda, NumericVector &yj,
		 IntegerVector &cens, NumericVector &limit,
		 DataFrame etasDf, List etaLst){
  NumericVector dv = clone(cdv);
  bool doCwres = (innerList.size() == 2);
  DataFrame pred = as<DataFrame>(innerList["pred"]);
  IntegerVector ID= as<IntegerVector>(pred[0]);
  IntegerVector TIME= as<IntegerVector>(pred[1]);
  pred.erase(0,2);
  NumericVector prednv = as<NumericVector>(pred[0]);
  NumericVector prednvI(prednv.size());
  pred.erase(0);
  // Individual predictions data frame.
  DataFrame ipred = as<DataFrame>(innerList["ipred"]);
  ipred.erase(0,2);
  NumericVector iprednv = as<NumericVector>(ipred[0]);
  NumericVector iprednvI(iprednv.size());
  NumericVector lowerLim(iprednv.size());
  NumericVector upperLim(iprednv.size());
  bool interestingLim=false;
  ipred.erase(0);
  // Now get fp r, rp
  unsigned int neta = omegaMat.nrow();
  NumericMatrix fpp(iprednv.size(),neta);
  NumericMatrix fpi(iprednv.size(),neta);
  NumericMatrix rpp(iprednv.size(),neta);
  NumericMatrix rpi(iprednv.size(),neta);
  unsigned int i, j;
  double om;
  unsigned int nid=etasDf.nrows();
  double tc = sqrt((double)nid);
  for (j = neta;j--;){
    NumericVector cur = etaLst[j];
    if (doCwres){
      fpp(_,j) = NumericVector(pred[j]);
      fpi(_,j) = NumericVector(ipred[j]);
    }
    // Finalize by adding shrinkage.
    om =omegaMat(j,j);
    cur[5]= (1.0-cur[1]/om)*100.0;
    cur[6]= (1.0-cur[2]/sqrt(om))*100.0;
    cur[7] = cur[0]*tc/(cur[2]);
    cur[8] = 2*Rf_pt(-fabs(cur[7]),(double)(nid-1),1,0);
  }
  CharacterVector dimN= (as<List>(omegaMat.attr("dimnames")))[1];
  CharacterVector dimN2(dimN.size()+1);
  dimN2[dimN.size()] = "IWRES";
  std::copy(dimN.begin(),dimN.end(),dimN2.begin());
  etaLst.attr("names") = dimN2;
  etaLst.attr("row.names") = CharacterVector::create("mean","var","sd","skewness", "kurtosis","var shrinkage (%)", "sd shrinkage (%)","t statistic","p-value");
  etaLst.attr("class") = "data.frame";
  if (doCwres){
    pred.erase(0,neta);
    ipred.erase(0,neta);
  }
  NumericVector rp = pred[0];
  arma::vec rpv = as<arma::vec>(rp);
  NumericVector ri = ipred[0];
  NumericVector dvTBS(dv.size());
  for (i = dv.size(); i--;){
    iprednvI[i]= powerDi(iprednv[i], lambda[i], (int)yj[i]);
    prednvI[i] = powerDi(prednv[i], lambda[i], (int)yj[i]);
    switch (cens[i]){
    case 1:
      // (limit, dv) ; limit could be -inf
      if (R_FINITE(limit[i])){
	double lim0TBS = powerD(limit[i], lambda[i], (int)yj[i]);
	double lim1TBS = powerD(dv[i], lambda[i], (int) yj[i]);
	double sd = sqrt(ri[i]);
	interestingLim=true;
	lowerLim[i] = limit[i];
	upperLim[i] = dv[i];
	dvTBS[i] = truncnorm(iprednv[i], sd, lim0TBS, lim1TBS);
	dv[i] =powerDi(dvTBS[i], lambda[i], (int) yj[i]);
      } else {
	// (-Inf, dv)
	double lim1TBS = powerD(dv[i], lambda[i], (int) yj[i]);
	double sd = sqrt(ri[i]);
	interestingLim=true;
	lowerLim[i] = R_NegInf;
	upperLim[i] = dv[i];
	dvTBS[i] = truncnorm(iprednv[i], sd, R_NegInf, lim1TBS);
	dv[i] =powerDi(dvTBS[i], lambda[i], (int) yj[i]);
      }
      break;
    case -1:
      // (dv, limit); limit could be +inf
      if (R_FINITE(limit[i])){
	//(dv, limit)
	double lim1TBS = powerD(limit[i], lambda[i], (int)yj[i]);
	double lim0TBS = powerD(dv[i], lambda[i], (int) yj[i]);
	double sd = sqrt(ri[i]);
	interestingLim=true;
	lowerLim[i] = dv[i];
	upperLim[i] = limit[i];
	dvTBS[i] = truncnorm(iprednv[i], sd, lim0TBS, lim1TBS);
	dv[i] =powerDi(dvTBS[i], lambda[i], (int) yj[i]);
      } else {
	// (dv, Inf)
	double lim1TBS = powerD(dv[i], lambda[i], (int) yj[i]);
	double sd = sqrt(ri[i]);
	interestingLim=true;
	lowerLim[i] = dv[i];
	upperLim[i] = R_PosInf;
	dvTBS[i] = truncnorm(iprednv[i], sd, lim1TBS, R_PosInf);
	dv[i] =powerDi(dvTBS[i], lambda[i], (int) yj[i]);
      }
      break;
    case 0:
      lowerLim[i] = NA_REAL;
      upperLim[i] = NA_REAL;
      dvTBS[i] = powerD(dv[i], lambda[i], (int)yj[i]);
      break;
    }
  }
  int warn1 = 0;
  int warn2 = 0;
  for (unsigned int j = ri.size(); j--;){
    if (ri[j] < DOUBLE_EPS){
      warn1 = 1;
      ri[j] = 1;
    }
    if (ISNAN(ri[j])){
      warn2 = 1;
      ri[j] = 1;
    }
  }
  arma::vec riv = as<arma::vec>(ri);
  DataFrame etasDf1 = etasDf;
  etasDf1.erase(0);
  NumericMatrix etas1(nid,neta);
  List etasDfFull(etasDf1.size());
  etasDfFull.names()=etasDf1.names();
  etasDfFull.attr("row.names")=IntegerVector::create(NA_INTEGER,-iprednv.size());
  etasDfFull.attr("class") = "data.frame";
  for (j = neta; j--;){
    NumericVector cur(iprednv.size(),NA_REAL);
    etasDfFull[j]= cur;
    etas1(_,j)=NumericVector(etasDf1[j]);
  }
  
  arma::mat etas = as<arma::mat>(etas1);
  arma::mat fppm = as<arma::mat>(fpp);
  arma::mat fpim = as<arma::mat>(fpi);
  arma::mat omegaMatA = as<arma::mat>(omegaMat);
  arma::vec Vfop, Vfoi;
  NumericVector dErr_dEta_i(fppm.n_rows, NA_REAL);
  if (!doCwres){
    int lastId = ID[ID.size()-1], lastCol = nid-1, lastIndex=ID.size()-1;
    int etaFulli = nid-1;
    double curEta=0.0;
    for (j = fppm.n_rows; j--; ){
      if (lastId != ID[j]){
        // Fill in full eta data frame
        for (i = neta; i--;){
          curEta = (as<NumericVector>(etasDf1[i]))[etaFulli];
          NumericVector cur = etasDfFull[i];
          std::fill_n(cur.begin()+j+1,lastIndex-j,curEta);
        }
        etaFulli--;
        lastId=ID[j];
        lastIndex=j;
        lastCol--;
        if (lastCol == 0){
          // Finalize ETA
          for (i = neta; i--;){
            curEta = (as<NumericVector>(etasDf1[i]))[0];
            NumericVector cur = etasDfFull[i];
            std::fill_n(cur.begin(),lastIndex+1,curEta);
          }
          break;
        }
      }
    }
  } else {
    arma::mat V_fo_p = (fppm * omegaMatA * fppm.t()); // From Mentre 2006 p. 352
    arma::mat V_fo_i = (fpim * omegaMatA * fpim.t()); // From Mentre 2006 p. 352
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
    Vfop = V_fo_p.diag();
    Vfoi = V_fo_i.diag();
    // Calculate dErr/dEta
    NumericVector dErr_dEta_p(fppm.n_rows, NA_REAL);
    int lastId = ID[ID.size()-1], lastCol = nid-1, lastIndex=ID.size()-1;
    // List etaFull(neta);
    int etaFulli = nid-1;
    double curEta=0.0;
    for (j = fppm.n_rows; j--; ){
      if (lastId != ID[j]){
        // Fill in full eta data frame
        for (i = neta; i--;){
          curEta = (as<NumericVector>(etasDf1[i]))[etaFulli];
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
          for (i = neta; i--;){
            curEta = (as<NumericVector>(etasDf1[i]))[0];
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
  // For WRES
  // W <- fitted(object, population=TRUE, type="Vfo") + fitted(object, population=TRUE, type="Vi")
  NumericVector resI = dv-prednvI;
  NumericVector res = dvTBS-prednv;
  arma::vec resv = as<arma::vec>(res);
  arma::vec wres, cwres;
  NumericVector cpredI, cresI, cres, cpred;
  if (doCwres){
    wres = resv/sqrt(abs(Vfop+rpv));
    // For CPRED, CRES and CWRES
    cpred = iprednv - dErr_dEta_i;
    cpredI = NumericVector(cpred.size());
    cres = dvTBS-cpred;
    cresI = NumericVector(cres.size());
    for (i = cres.size(); i--;){
      cpredI[i] = powerDi(cpred[i], lambda[i], (int)yj[i]);
      cresI[i]  = dv[i] - cpredI[i];
    }
    //W <- sqrt(fitted(object, population=FALSE, type="Vfo") + fitted(object, population=FALSE, type="Vi"))
    arma::vec cresv=as<arma::vec>(cres);
    cwres = cresv/sqrt(Vfoi+riv);
  }
  // e["Vfo"] = wrap(Vfo);
  // e["dErr_dEta"] = wrap(dErr_dEta);
  // W <- fitted(object, population=TRUE, type="Vfo") + fitted(object, population=TRUE, type="Vi")
  // Efo= f(|eta=0)
  // cov = Vfo|eta=0+Vi|eta=0
  // These WRES are calculated the same as Hooker 2007, but don't seem to match NONMEM.
  //  In theory the decorrelation is done by eigenvector decomposition. However, it doens't seem to affect the WRES values.
  //  W <- fitted(object, population=TRUE, type="Vfo") + fitted(object, population=TRUE, type="Vi")
  // tmp <- aggregate(W,by=list(object$ID),FUN=function(x){e <- eigen(diag(x)); return(diag(e$vectors %*% diag(sqrt(e$values)) %*% rxInv(e$vectors)))});
  // tmp <- data.frame(ID=tmp[,1],stack(data.frame(tmp[,-1])));
  // tmp <- tmp[order(tmp$ID,tmp$ind),];
  // tmp <- sqrt(tmp$values)
  // Even though it does not match NONMEM, it should be adequate metric for WRES...
  // pred.erase(0,neta);
  // ipred.erase(0,neta);
  // rp(_,0) = NumericVector(pred[0]);
  // ri(_,0) = NumericVector(ipred[0]);
  
  NumericVector iwres=(dvTBS-iprednv)/sqrt(ri);
  NumericVector stat=etaLst[neta];
  unsigned int n =0, n1 = 0;
  double M1 =  0, M2 = 0, M3 = 0, M4 =0, term1, delta, delta_n, delta_n2;
  for (unsigned int i = iprednv.size(); i--;){
    if (evid[i] == 0){
      n1=n; n++;
      delta = iwres[i] - M1;
      delta_n = delta / n;
      delta_n2 = delta_n * delta_n;
      term1 = delta * delta_n * n1;
      M1 += delta_n;
      M4 += term1 * delta_n2 * (n*n - 3*n + 3) + 6 * delta_n2 * M2 - 4 * delta_n * M3;
      M3 += term1 * delta_n * (n - 2) - 3 * delta_n * M2;
      M2 += term1;    
    } else {
      iwres[i] = NA_REAL;
    }
  }
  stat[0] = M1;
  stat[1] = M2/((double)n1);
  stat[2] = sqrt(stat[1]);
  stat[3] = sqrt((double)(n)) * M3/ pow(M2, 1.5);
  stat[4] = (double)(n)*M4 / (M2*M2) - 3.0;
  stat[5] = (1-stat[1])*100;
  stat[6] = (1-stat[2])*100;
  stat[7] = M1*sqrt((double)n)/stat[2];
  stat[8] = 2*Rf_pt(stat[7],(double)n1,1,0);
  List ret(4);
  // Now do inverse TBS to backtransform
  if (doCwres){
    NumericVector ires = dv-iprednvI;
    for (i = cwres.size(); i--;){
      if (evid[i] != 0){
	wres[i] = NA_REAL;
	cwres[i] = NA_REAL;
	resI[i] = NA_REAL;
	ires[i] = NA_REAL;
      }
    }
    if (interestingLim){
      ret[0] = DataFrame::create(_["PRED"]=prednvI,
				 _["RES"]=resI,
				 _["WRES"]=wrap(wres),
				 _["IPRED"]=iprednvI,
				 _["IRES"]=ires,
				 _["IWRES"]=iwres,
				 _["CPRED"]=cpredI,
				 _["CRES"]=cresI,
				 _["CWRES"]=cwres,
				 _["lowerLim"] = lowerLim,
				 _["upperLim"] = upperLim);
    } else {
      ret[0] = DataFrame::create(_["PRED"]=prednvI,
				 _["RES"]=resI,
				 _["WRES"]=wrap(wres),
				 _["IPRED"]=iprednvI,
				 _["IRES"]=ires,
				 _["IWRES"]=iwres,
				 _["CPRED"]=cpredI,
				 _["CRES"]=cresI,
				 _["CWRES"]=cwres);
    }
  } else {
    NumericVector ires = dv-iprednvI;
    for (i = ires.size(); i--;){
      if (evid[i] != 0){
	resI[i] = NA_REAL;
	ires[i] = NA_REAL;
      }
    }
    if (interestingLim){
      ret[0] = DataFrame::create(_["PRED"]=prednvI,
				 _["RES"]=resI,
				 _["IPRED"]=iprednvI,
				 _["IRES"]=ires,
				 _["IWRES"]=iwres,
				 _["lowerLim"] = lowerLim,
				 _["upperLim"] = upperLim);
    } else {
      ret[0] = DataFrame::create(_["PRED"]=prednvI,
				 _["RES"]=resI,
				 _["IPRED"]=iprednvI,
				 _["IRES"]=ires,
				 _["IWRES"]=iwres);
    }
  }
  ret[1] = etaLst;
  ret[2] = etasDfFull;
  ret[3] = dv;
  if (warn1){
    warning("Some variances were zero, replaced with 1.");
  }
  if (warn2){
    warning("Some variances were NA/NaN, replaced with 1.");
  }
  return ret;
}

arma::mat cholSE__(arma::mat A, double tol);

// [[Rcpp::export]]
List npde(IntegerVector id, NumericVector dv, IntegerVector evid,
	  NumericVector sim, NumericVector lambda, NumericVector yj,
	  bool ties, double tolChol){
  bool warn1 = false;
  bool warn2 = false;
  bool warn3=false;
  int nobs = dv.size();
  int K = sim.size()/nobs;
  NumericVector epredt(nobs);
  NumericVector epred(nobs);
  NumericVector eres(nobs);
  NumericVector erest(nobs);
  NumericVector simrest(sim.size());
  // NumericVector simY(sim.size());
  // NumericVector Y(nobs);
  NumericVector pd(nobs);
  NumericVector npde(nobs);
  NumericVector nsim(nobs);
  NumericVector dvt(nobs);
  IntegerVector idLoc(nobs);
  // This uses welford's method for means
  int lastId=id[0], nid=1;
  idLoc[0]=0;
  int k2=0;
  for (int i = 0; i < nobs; i++){
    if (lastId != id[i]){
      idLoc[nid] = i;
      nid++;
      lastId = id[i];
    }
    // Calculte EPREDt
    if (ISNA(sim[i])){
      epredt[i] = 0;
      warn3=true;
    } else {
      k2++;
      epredt[i] = sim[i];
    }
    for (int j = 1; j < K; j++){
      if (!ISNA(sim[i+nobs*j])){
	k2++;
	epredt[i] = epredt[i]+(sim[i+nobs*j]-epredt[i])/((double)(k2));
      } else {
	warn3=true;
      }
    }
    // Calculate dvt, EPRED, ERES ERESt
    dvt[i] = powerD(dv[i], lambda[i], (int)yj[i]);
    epred[i] = powerDi(epredt[i], lambda[i], (int)yj[i]);
    eres[i] = dv[i] - epred[i];
    if (evid[i] != 0) eres[i] = NA_REAL;
    erest[i] = dvt[i] - epredt[i];
    // Calculate simRESt
    for (int j = K; j--; ){
      simrest[i+nobs*j] = sim[i+nobs*j]-epredt[i];
    }
  }
  idLoc[nid] = nobs;
  // Calculate Var(Yj)
  // Even though we are using matrix algebra, we are still using Welford's method
  double d;
  for (int i = 0; i < nid; i++){
    arma::vec d1(idLoc[i+1]-idLoc[i]);
    std::copy(simrest.begin()+idLoc[i], simrest.begin()+idLoc[i+1], d1.begin());
    arma::ivec evid1(idLoc[i+1]-idLoc[i]);
    std::copy(evid.begin()+idLoc[i], evid.begin()+idLoc[i+1], evid1.begin());
    arma::uvec evid0 = find(evid1 == 0);
    arma::vec d2 = d1(evid0);
    bool anyNa = false;
    for (int j = d2.size(); j--;){
      if (ISNAN(d2[j])){
        anyNa=true;
	break;
      }
    }
    mat vYi(d2.size(), d2.size(), fill::zeros);
    k2=0;
    if (anyNa) {
      warn3=true;
    } else {
      k2++;
      vYi = d2 * trans(d2);
    }
    mat vYi2;
    for (int j = 1; j < K; j++){
      std::copy(simrest.begin()+idLoc[i]+nobs*j, simrest.begin()+idLoc[i+1]+nobs*j, d1.begin());
      d2 = d1(evid0);
      anyNa=false;
      for (int k = d2.size(); k--;){
        if (ISNAN(d2[k])){
          anyNa=true;
          break;
        }
      }
      if (anyNa){
	warn3=true;
      } else {
        k2++;
        d=k2;
        vYi2=d2 * trans(d2);
        vYi = vYi + (vYi2 - vYi)/d;
      }
    }
    // Now decorrelate covariance
    vec eigval;
    mat eigvec;
    vec d;
    mat u;
    mat v;
    arma::mat ch;
    mat eigvalS(vYi.n_rows, vYi.n_rows, fill::zeros);
    try {
      ch = chol(vYi);
    } catch (...){
      // Try nearPD from R/RxODE
      warn1 = true;
      try {
	ch = cholSE__(vYi, tolChol);
      } catch (...){
	stop("The simulated data covariance matrix Cholesky decomposition could not be obtained for %d.", i+1);
      }
    }
    try{
      vYi = trans(inv(trimatu(ch)));
    } catch (...) {
      warn2=true;
      try {
	vYi = trans(pinv(trimatu(ch)));
      } catch (...){
	stop("The simulated data covariance matrix Cholesky decomposition inverse or pseudo-inverse could not be obtained for %d.", i+1);
      }
    }
    // Now calculate  Y_{i,sim} and Y_i
    std::copy(erest.begin()+idLoc[i], erest.begin()+idLoc[i+1], d1.begin());
    d2 = d1(evid0);
    mat Yi = vYi*d2;
    mat Yis;
    std::copy(simrest.begin()+idLoc[i], simrest.begin()+idLoc[i+1], d1.begin());
    d2 = d1(evid0);
    Yis = vYi*d2;
    int kk = idLoc[i];
    for (int k = idLoc[i]; k < idLoc[i+1]; k++){
      if (evid1[k-idLoc[i]] == 0){
	if (Yi[kk-idLoc[i]] < Yis[kk-idLoc[i]]){
	  pd[k] = 1.0;
	} else {
	  pd[k] = 0.0;
	}
	kk++;
      } else {
	pd[k] = NA_REAL;
      }
    }
    double de;
    k2 = 0;
    for (int j = 1; j < K; j++){
      std::copy(simrest.begin()+idLoc[i]+j*nobs, simrest.begin()+idLoc[i+1]+j*nobs, d1.begin());
      d2 = d1(evid0);
      Yis = vYi*d2;
      int kk = idLoc[i];
      for (int k = idLoc[i]; k < idLoc[i+1]; k++){
	if (evid1[k-idLoc[i]] == 0){
	  if (ISNAN(Yis[kk-idLoc[i]])){
	    warn3=true;
	  } else if (Yi[kk-idLoc[i]] < Yis[k-idLoc[i]]){
	    k2++;
	    de = k2;
	    pd[k] = pd[k]+(1-pd[k])/de;
	  } else {
	    k2++;
	    de = k2;
	    pd[k] = pd[k]-pd[k]/de;
	  }
	  kk++;
	} else {
	  pd[k] = NA_REAL;
	}
      }
    }
  }
  // Now deal with exceptions pd=1 or 0 or when ties=FALSE, put some
  // random noise in.
  if (ties){
    for (int i = pd.size(); i--;){
      if (evid[i] == 0){
	if (fabs(pd[i]) < DOUBLE_EPS){
	  pd[i] = 1 / (2.0 * K);
	} else if (fabs(1-pd[i]) < DOUBLE_EPS){
	  pd[i] = 1 - 1 / (2.0 * K);
	}
	npde[i] = -qnorm(pd[i], 0.0, 1.0, true, false);
      } else {
	pd[i] = NA_REAL;
	npde[i] = NA_REAL;
      }
    }
  } else {
    // This is a little different than the NPDE manual
    // the bounds are 1/(2K) to 1 - 1/(2K) instead of 0 to 1
    for (int i = pd.size(); i--;){
      if (evid[i] == 0){
	if (fabs(pd[i]) < DOUBLE_EPS){
	  pd[i] = Rf_runif(0, 1/K);
	} else if (fabs(1-pd[i]) < DOUBLE_EPS){
	  pd[i] = Rf_runif(1.0 - 1/K, 1.0);
	} else {
	  pd[i] = Rf_runif(0, 1/K) + pd[i];
	}
	if (fabs(pd[i]) < 1 / (2.0 * K)){
	  pd[i] = 1 / (2.0 * K);
	} else if (fabs(1-pd[i]) < 1 / (2.0 * K)){
	  pd[i] = 1 - 1 / (2.0 * K);
	}
	npde[i] = -qnorm(pd[i], 0.0, 1.0, true, false);
      } else {
	npde[i] = NA_REAL;
	pd[i] = NA_REAL;
      }
    }
  }
  List ret(3);
  ret[0] = epred;
  ret[1] = eres;
  ret[2] = npde;
  ret.attr("names") = CharacterVector::create("EPRED", "ERES", "NPDE");
  ret.attr("row.names") = IntegerVector::create(NA_INTEGER, -nobs);
  ret.attr("class") = "data.frame";
  if (warn1){
    warning("Some subjects simulated data covariance matrix used a generalized Cholesky decomposition.");
  }
  if (warn2){
    warning("Some subjects simulated data covariance matrix were inverted by a psuedo-inverse");
  }
  if (warn3){
    warning("Some of the simulated subjects had NA/NaN values, implying failed solves and likely an unstable model");
  }
  return ret;
}

//[[Rcpp::export]]
RObject augPredTrans(NumericVector& pred, NumericVector& ipred, NumericVector& lambda,
		     RObject& yjIn){
  IntegerVector yj = as<IntegerVector>(yjIn);
  for (int i = pred.size(); i--;){
    pred[i] = powerDi(pred[i], lambda[i], yj[i]);
    ipred[i] = powerDi(ipred[i], lambda[i], yj[i]);
  }
  return R_NilValue;
}
