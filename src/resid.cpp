#include <stdarg.h>
#include <RcppArmadillo.h>
#include <R.h>
#include <Rmath.h>
#include "nlmixr_types.h"

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
List nlmixrResid(List &innerList, NumericMatrix &omegaMat, NumericVector &dv, DataFrame etasDf, List etaLst){
  DataFrame pred = as<DataFrame>(innerList["pred"]);
  IntegerVector ID= as<IntegerVector>(pred[0]);
  IntegerVector TIME= as<IntegerVector>(pred[1]);
  pred.erase(0,2);
  NumericVector prednv = as<NumericVector>(pred[0]);
  pred.erase(0);
  DataFrame ipred = as<DataFrame>(innerList["ipred"]);
  ipred.erase(0,2);
  NumericVector iprednv = as<NumericVector>(ipred[0]);
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
    fpp(_,j) = NumericVector(pred[j]);
    fpi(_,j) = NumericVector(ipred[j]);
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
  pred.erase(0,neta);
  ipred.erase(0,neta);
  NumericVector rp = pred[0];
  arma::vec rpv = as<arma::vec>(rp);
  NumericVector ri = ipred[0];
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
  arma::vec Vfop = V_fo_p.diag();
  arma::vec Vfoi = V_fo_i.diag();
  // Calculate dErr/dEta
  NumericVector dErr_dEta_p(fppm.n_rows, NA_REAL);
  NumericVector dErr_dEta_i(fppm.n_rows, NA_REAL);
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
  // For WRES
  // W <- fitted(object, population=TRUE, type="Vfo") + fitted(object, population=TRUE, type="Vi")
  NumericVector res = dv-prednv;
  arma::vec resv = as<arma::vec>(res);
  arma::vec wres = resv/sqrt(Vfop+rpv);
  // For CPRED, CRES and CWRES
  NumericVector cpred = iprednv - dErr_dEta_i;
  NumericVector cres = dv-cpred;
  //W <- sqrt(fitted(object, population=FALSE, type="Vfo") + fitted(object, population=FALSE, type="Vi"))
  arma::vec cresv=as<arma::vec>(cres);
  arma::vec cwres = cresv/sqrt(Vfoi+riv);
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

  NumericVector iwres=(dv-iprednv)/sqrt(ri);
  NumericVector stat=etaLst[neta];
  unsigned int n =0, n1 = 0;
  double M1 =  0, M2 = 0, M3 = 0, M4 =0, term1, delta, delta_n, delta_n2;
  for (unsigned int i = iprednv.size(); i--;){
    n1=n; n++;
    delta = iwres[i] - M1;
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
  stat[2] = sqrt(stat[1]);
  stat[3] = sqrt((double)(n)) * M3/ pow(M2, 1.5);
  stat[4] = (double)(n)*M4 / (M2*M2) - 3.0;
  stat[5] = (1-stat[1])*100;
  stat[6] = (1-stat[2])*100;
  stat[7] = M1*sqrt((double)n)/stat[2];
  stat[8] = 2*Rf_pt(stat[7],(double)n1,1,0);
  List ret(3);
  ret[0] = DataFrame::create(_["PRED"]=prednv,
                             _["RES"]=res,
                             _["WRES"]=wrap(wres),
                             _["IPRED"]=iprednv,
                             _["IRES"]=dv-iprednv,
                             _["IWRES"]=iwres,
                             _["CPRED"]=cpred,
                             _["CRES"]=cres,
                             _["CWRES"]=cwres);
  ret[1] = etaLst;
  ret[2] = etasDfFull;
  return ret;
}
