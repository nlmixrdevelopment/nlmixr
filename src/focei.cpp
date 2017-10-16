// [[Rcpp::depends(RcppArmadillo)]]
#include <stdarg.h>
#include <RcppArmadillo.h>
#include <R.h>
#include "nlmixr_types.h"

extern "C" double RxODE_sumV(int n, ...);
extern "C" double RxODE_prodV(int n, ...);
extern "C" double RxODE_sum (double *input, int n);

#define _prod RxODE_prodV
#define _sum  RxODE_sumV

using namespace Rcpp;
using namespace R;
using namespace arma;

extern "C" SEXP RxODE_ode_solver_focei_eta (SEXP sexp_eta, SEXP sexp_rho);
extern "C" int nEq();
extern "C" unsigned int nLhs();
extern "C" unsigned int nAllTimes();
extern "C" int rxEvid(int i);
extern "C" double rxLhs(int i);
extern "C" void rxCalcLhs(int i);
extern "C" unsigned int nObs();
extern "C" void RxODE_ode_solve_env(SEXP sexp_rho);
extern "C" void RxODE_ode_free();
extern "C" double RxODE_safe_zero(double x);
extern "C" double RxODE_safe_log(double x);
extern "C" double RxODE_sign_exp(double sgn, double x);
extern "C" double RxODE_abs_log(double x);
void rxDetaDomega(SEXP rho);
double stablizeNums(double x){
  if (ISNAN(x)){
    return 0.0;
  } else if (ISNA(x)){
    return 0.0;
  } else if (!R_FINITE(x)){
    if (x < 0){
      return -42e100;
    } else {
      return 42e100;
    }
  } else {
    return x;
  }
}

double Lhs(int i){
  return stablizeNums(rxLhs(i));
}
//' Returns the gradient for FOCEi environment
//' @param rho Environment to calculate the gradient for...
//' @export
//' @keywords internal
// [[Rcpp::export]]
void rxGrad(SEXP rho){
  Environment e = as<Environment>(rho);
  mat err = mat(nObs(),1);
  mat r = mat(nObs(),1);
  int prd = as<int>(e["pred.minus.dv"]);
  // unsigned int neta = as<unsigned int>(e["neta"]);
  unsigned int ntheta = as<unsigned int>(e["ntheta"]);
  NumericVector DV = as<NumericVector>(e["DV"]);
  NumericVector f(nObs());
  mat fpm = mat(nObs(), ntheta);
  mat rp = mat(nObs(),ntheta);
  mat lp = mat(ntheta,1);
  mat B = mat(nObs(),1);
  List c(ntheta);
  // List a(ntheta);
  mat cur;
  // int do_nonmem = as<int>(e["nonmem"]);
  unsigned int i = 0, j = 0, k = 0;
  RxODE_ode_solve_env(rho);
  double *lpi =Calloc(ntheta*nObs()*3,double);
  for (i = 0; i < ntheta; i++){
    lp[i] = 0;
    // a[i] = mat(nObs(),1);
    c[i] = mat(nObs(),1);
  }
  for (i = 0; i < nAllTimes(); i++){
    if (!rxEvid(i)){
      // Rprintf("%d:\n",i);
      // Rprintf("F\n");
      rxCalcLhs(i);
      f[k] = Lhs(0); // Pred
      if (prd == 1){
        err(k, 0) =  f[k] - DV[k];
      } else {
        err(k, 0) =  DV[k] - f[k];
      }
      // Rprintf("dpred\n");
      // d(pred)/d(theta#)
      for (j = 1; j < ntheta+1; j++){
        fpm(k, j-1) = Lhs(j);
        // cur = as<mat>(a[j-1]);
        // if (do_nonmem){
        //   cur(k,0) =  Lhs(j);
        // } else {
        //   cur(k,0) = Lhs(j)-err(k, 0)/RxODE_safe_zero(Lhs(neta+1))*Lhs(j+neta+1);
        // }
        // a[j-1]=cur;
      }
      // R
      // Rprintf("R\n");
      if (Lhs(j) < 0){
        Rprintf("R term is %d\n",j);
        for (j = 0; j < nLhs(); j++){
          Rprintf("Lhs(%d) = %f\n", j, Lhs(j));
        }
        Rprintf("\n");
        // temp = getAttrib(sexp_theta, R_NamesSymbol);
        // for (j = 0; j < length(sexp_theta); j++){
        //   Rprintf("params[%d] = %f\n", j, par_ptr[j]);
        // }
        RxODE_ode_free();
        stop("A covariance term is zero or negative and should remain positive");
      }
      r(k, 0)=Lhs(j); // R always has to be positive.
      B(k, 0)=2/RxODE_safe_zero(Lhs(j));
      // Rprintf("dR\n");
      for (j=ntheta+2; j < nLhs(); j++){
        rp(k,j-ntheta-2) = Lhs(j);
        cur = as<mat>(c[j-ntheta-2]);
        cur(k,0) = Lhs(j)/RxODE_safe_zero(r(k, 0));
        c[j-ntheta-2] = cur;
      }
      for (j = 0; j < ntheta; j++){
        // .5*apply(eps*fp*B + .5*eps^2*B*c - c, 2, sum) - OMGAinv %*% ETA
        // Rprintf("lp[%d] ",j);
        cur = as<mat>(c[j]);
	lpi[k         +3*nObs()*j] = _prod(5, 0.25, err(k, 0), err(k, 0), B(k, 0), cur(k,0));
        lpi[k+nObs()  +3*nObs()*j] = - _prod(2, 0.5, cur(k,0));
        lpi[k+2*nObs()+3*nObs()*j] = - _prod(4, 0.5, err(k, 0), fpm(k, j), B(k, 0));
      }
      // Rprintf("done\n");
      k++;
    }
  }
  // Deferred sum to reduce round-off error.
  for (j = 0; j < ntheta; j++){
    // Rprintf("%f", eta[j]);
    lp[j] = RxODE_sum(lpi+3*nObs()*j,3*nObs());
  }
  // Free
  RxODE_ode_free();
  Free(lpi);
  e["err"] = err;
  e["f"] = f;
  e["dErr"] = fpm;
  e["dR"] = rp;
  e["R"] = r;
  e["lp"] = -lp;
}

// [[Rcpp::export]]
void rxInnerNum(SEXP etanews, SEXP rho){
  Environment e = as<Environment>(rho);
  if (!(e.exists("neta") && e.exists("ntheta") && e.exists("dOmega") &&
        e.exists("DV") && e.exists("nonmem") && e.exists("eta") &&
        e.exists("eta.mat") && e.exists("eta.trans") &&
        e.exists("params")
        )){
    stop("Environment not setup correctly for rxInnerNum.");
  }
  NumericVector par_ptr = as<NumericVector>(e["params"]);
  IntegerVector eta_i = as<IntegerVector>(e["eta.trans"]);
  NumericVector etanew = as<NumericVector>(etanews);
  NumericVector eta = as<NumericVector>(e["eta"]);

  mat etam = as<mat>(e["eta.mat"]);
  unsigned int recalc = 0;
  unsigned int i = 0, j = 0, k = 0;

  if (!e.exists("llik")){
    recalc = 1;
  } else if (eta.size() != etanew.size()){
    stop("Inconsistent eta size for rxInnerNum.");
  } else {
    for (i = 0; i < (unsigned int)(eta.size()); i++){
      if (eta[i] != etanew[i]){
        recalc = 1;
        break;
      }
    }
  }

  mat cur;
  if (recalc){
    // Rprintf("Start!\n");
    for (i = 0; i < (unsigned int)(etanew.size()); i++){
      /* Rprintf("\tpar[%d] from %f to %f\n",eta_i[i],par_ptr[eta_i[i]], eta[i]); */
      par_ptr[eta_i[i]] = etanew[i];
      eta[i] = etanew[i];
      etam(i,0) = eta[i];
    }
    e["params"] = par_ptr;
    e["eta"]    = eta;
    e["eta.mat"] = etam;

    e["params.0"] = par_ptr;
    e["eta.0"]    = eta;
    e["eta.mat.0"] = etam;
    
    RxODE_ode_solve_env(rho);
    // In this case F = lhs(0)
    // and R = lhs(1)
    NumericVector DV = as<NumericVector>(e["DV"]);
  
    NumericVector f(nObs());
    mat err = mat(nObs(),1);
    mat r = mat(nObs(),1);
  
    NumericVector llik(1);
  
    double *llik0 =Calloc(nObs()*2,double);

    // Get the base F, Err and R
    int prd = as<int>(e["pred.minus.dv"]);
    for (i = 0; i < nAllTimes(); i++){
      if (!rxEvid(i)){
	rxCalcLhs(i);
	f[k] = Lhs(0); // Pred
	if (prd == 1){
	  err(k, 0) =  f[k] - DV[k];
	} else {
	  err(k, 0) =  DV[k] - f[k];
	}
	if (Lhs(1) < 0){
	  for (j = 0; j < nLhs(); j++){
	    Rprintf("Lhs(%d) = %f\n", j, Lhs(j));
	  }
	  Rprintf("\n");
	  // temp = getAttrib(sexp_theta, R_NamesSymbol);
	  // for (j = 0; j < length(sexp_theta); j++){
	  //   Rprintf("params[%d] = %f\n", j, par_ptr[j]);
	  // }
	  RxODE_ode_free();
	  stop("A covariance term is zero or negative and should remain positive");
	}
	r(k, 0)=Lhs(1); // R always has to be positive.
	k++;
      }
    }
    // Rprintf("Done getting F and R!\n");
    // Now get the dErr/d(eta); and dR/d(eta) by finitie difference (simple forward difference)
    double eps = 1e-4;
    unsigned int neta = as<unsigned int>(e["neta"]);
  
    mat fpm = mat(nObs(), neta); // d(pred)/d(eta#)
    mat rp = mat(nObs(),neta);

    unsigned int curEta;
    double tmp;
    for (curEta = 0; curEta < neta; curEta++){
      for (i = 0; i < (unsigned int)(etanew.size()); i++){
	/* Rprintf("\tpar[%d] from %f to %f\n",eta_i[i],par_ptr[eta_i[i]], eta[i]); */
	if (i == curEta){
	  // Forward difference
	  par_ptr[eta_i[i]] = etanew[i]+eps;
	} else {
	  par_ptr[eta_i[i]] = etanew[i];
	}
	eta[i] = etanew[i];
	etam(i,0) = eta[i];
      }
      e["params"] = par_ptr;
      e["eta"]    = eta;
      e["eta.mat"] = etam;

      e["params.0"] = par_ptr;
      e["eta.0"]    = eta;
      e["eta.mat.0"] = etam;
    
      RxODE_ode_solve_env(rho);
      // Now populate dpred/deta# and dR/deta#
      k = 0;
      if (!rxEvid(i)){
	rxCalcLhs(i);
	tmp = Lhs(0)-f[k];
	if (tmp == 0){
	  tmp = 0;
	} else if (tmp  < 0){
	  tmp = -exp(log(-tmp)-log(eps));
	} else {
	  tmp = exp(log(tmp)-log(eps));
	}
	fpm(k, curEta) = stablizeNums(tmp);
	tmp = Lhs(1)-r(k, 0);
	if (tmp == 0){
	  tmp  = 0;
        } else if (tmp < 0){
	  tmp  = -exp(log(-tmp)-log(eps));
	} else {
	  tmp  = exp(log(tmp)-log(eps));
        }
	rp(k, curEta) = stablizeNums(tmp);
	  k++;
      }
    }
    // Rprintf("Calculated fpm and rp!\n");
    // Now calculate the others like rxInner with sensitivity equations.
    mat B = mat(nObs(),1);
    List c(neta);
    List a(neta);
  
    mat lp = mat(neta,1);
    
    double *lpi =Calloc(neta*nObs()*3,double);
  
    for (i = 0; i < neta; i++){
      a[i] = mat(nObs(),1);
      c[i] = mat(nObs(),1);
      lp[i] = 0;
    }
    int do_nonmem = as<int>(e["nonmem"]);
    int do_table = as<int>(e["table"]);
    k = 0;
    for (i = 0; i < nAllTimes(); i++){
      if (!rxEvid(i)){
	B(k,0) = 2.0/RxODE_safe_zero(r(k, 0));
	// d(pred)/d(eta#)
	for (j = 0; j < neta; j++){
	  cur = as<mat>(a[j]);
	  if (do_nonmem){
	    cur(k,0) =  fpm(k, j);
	  } else {
	    cur(k,0) = fpm(k, j)-err(k, 0)/RxODE_safe_zero(r(k, 0))*rp(k, j);
	  }
	  a[j]=cur;
	  cur = as<mat>(c[j]);
	  cur(k, 0) = rp(k, j)/r(k, 0);
	  c[j] = cur;
	  lpi[k         +3*nObs()*j] = _prod(5, 0.25, err(k, 0), err(k, 0), B(k, 0), cur(k,0));
	  lpi[k+nObs()  +3*nObs()*j] = - _prod(2, 0.5, cur(k,0));
	  lpi[k+2*nObs()+3*nObs()*j] = - _prod(4, 0.5, err(k, 0), fpm(k, j), B(k, 0));
	}
	llik0[k] = -_prod(4,0.5,err(k, 0),err(k, 0),1.0/RxODE_safe_zero(r(k, 0)));
	llik0[k+nObs()] = -_prod(2, 0.5, RxODE_safe_log(r(k, 0)));
	k++;
      }
    }
    // Rprintf("Calculated B, c, a and llik pieces...\n");
    // Deferred sums (to reduce round-off error)
    for (j = 0; j < neta; j++){
      // Rprintf("%f", eta[j]);
      lp[j] = RxODE_sum(lpi+3*nObs()*j,3*nObs());
    }
    llik[0] = RxODE_sum(llik0,2*nObs());

    
    Free(lpi);
    Free(llik0);

    // Free
    RxODE_ode_free();
  
    mat llikm = mat(1,1);
    
    mat omegaInv = as<mat>(e["omegaInv"]);

    NumericVector llik2(1);
    llikm(0, 0)=  llik[0];
    llikm = -(llikm - 0.5*(etam.t() * omegaInv * etam));
    llik2[0] = llikm(0, 0);

    mat ep2 = -(lp - omegaInv * etam);
    if (do_table){
      mat omega = as<mat>(e["omega"]);
      mat Vfo_full = (fpm * omega * fpm.t()); // From Mentre 2006 p. 352
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
      // The Vfo calculation is separated out so it can be used in CWRES and WRES.
    
      mat Vfo = Vfo_full.diag();
      mat dErr_dEta = fpm * etam;
      e["Vfo"] = wrap(Vfo);
      e["dErr_dEta"] = wrap(dErr_dEta);
    }
    
    // Assign in env
    e["err"] = err;
    e["f"] = f;
    e["dErr"] = fpm;
    e["dR"] = rp;
    e["c"] = c;
    e["R"] = r;
    e["B"] = B;
    e["a"] = a;
    e["llik"] = llik;
    e["lp"] = lp;
    e["llik2"] = wrap(llik2);
    e["ep2"] = wrap(ep2);
  }
}

// [[Rcpp::export]]
void rxInner(SEXP etanews, SEXP rho){
  Environment e = as<Environment>(rho);
  if (!(e.exists("neta") && e.exists("ntheta") && e.exists("dOmega") &&
        e.exists("DV") && e.exists("nonmem") && e.exists("eta") &&
        e.exists("eta.mat") && e.exists("eta.trans") &&
        e.exists("params")
        )){
    stop("Environment not setup correctly for rxInner.");
  }
  if (as<int>(e["numeric"]) == 1){
    rxInnerNum(etanews, rho);
    return;
  }
  NumericVector par_ptr = as<NumericVector>(e["params"]);
  IntegerVector eta_i = as<IntegerVector>(e["eta.trans"]);
  NumericVector etanew = as<NumericVector>(etanews);
  NumericVector eta = as<NumericVector>(e["eta"]);
  mat etam = as<mat>(e["eta.mat"]);
  unsigned int recalc = 0;
  unsigned int i = 0, j = 0, k = 0;
  if (!e.exists("llik")){
    recalc = 1;
  } else if (eta.size() != etanew.size()){
    stop("Inconsistent eta size for rxInner.");
  } else {
    for (i = 0; i < (unsigned int)(eta.size()); i++){
      if (eta[i] != etanew[i]){
        recalc = 1;
        break;
      }
    }
  }
  if (recalc){
    for (i = 0; i < (unsigned int)(etanew.size()); i++){
      /* Rprintf("\tpar[%d] from %f to %f\n",eta_i[i],par_ptr[eta_i[i]], eta[i]); */
      par_ptr[eta_i[i]] = etanew[i];
      eta[i] = etanew[i];
      etam(i,0) = eta[i];
    }
    e["params"] = par_ptr;
    e["eta"]    = eta;
    e["eta.mat"] = etam;
    
    RxODE_ode_solve_env(rho);
    
    unsigned int neta = as<unsigned int>(e["neta"]);
    List dOmega = as<List>(e["dOmega"]);
    NumericVector DV = as<NumericVector>(e["DV"]);
    int do_nonmem = as<int>(e["nonmem"]);
    int do_table = as<int>(e["table"]);
    
    mat fpm = mat(nObs(), neta); // d(pred)/d(eta#)
    
    mat rp = mat(nObs(),neta);
    
    NumericVector f(nObs());
    mat err = mat(nObs(),1);
    mat r = mat(nObs(),1);

    mat B = mat(nObs(),1);
    List c(neta);
    List a(neta);
  
    NumericVector llik(1);
    mat lp = mat(neta,1);
    
    double *lpi =Calloc(neta*nObs()*3,double);
    double *llik0 =Calloc(nObs()*2,double);

    for (i = 0; i < neta; i++){
      a[i] = mat(nObs(),1);
      c[i] = mat(nObs(),1);
      lp[i] = 0;
    }

    
    /* // Now create the pred vector and d(pred)/d(eta) matrix. */
    /* // Assuming Lhs(0) = pred and Lhs(1:n) = d(pred)/d(eta#) */
    // Solve
    mat cur;
    int prd = as<int>(e["pred.minus.dv"]);
    for (i = 0; i < nAllTimes(); i++){
      if (!rxEvid(i)){
        rxCalcLhs(i);
        f[k] = Lhs(0); // Pred
	// if (ISNAN(f[k])){
	//   Rprintf("DV[k]: %f\n",DV[k]);
	//   for (j = 0; j < neta; j++){
	//     Rprintf("\t: ETA[%d] =%f\n",j+1,eta[j]);
	//   }
	//   stop("NaN in pred.");
	// }
        if (prd == 1){
          err(k, 0) =  f[k] - DV[k];
        } else {
          err(k, 0) =  DV[k] - f[k];
        }
        // d(pred)/d(eta#)
        for (j = 1; j < neta+1; j++){
          fpm(k, j-1) = Lhs(j);
          cur = as<mat>(a[j-1]);
          if (do_nonmem){
            cur(k,0) =  Lhs(j);
          } else {
            cur(k,0) = Lhs(j)-err(k, 0)/RxODE_safe_zero(Lhs(neta+1))*Lhs(j+neta+1);
          }
          a[j-1]=cur;
        }
        if (Lhs(j) < 0){
          for (j = 0; j < nLhs(); j++){
            Rprintf("Lhs(%d) = %f\n", j, Lhs(j));
          }
          Rprintf("\n");
          // temp = getAttrib(sexp_theta, R_NamesSymbol);
          // for (j = 0; j < length(sexp_theta); j++){
          //   Rprintf("params[%d] = %f\n", j, par_ptr[j]);
          // }
          RxODE_ode_free();
          stop("A covariance term is zero or negative and should remain positive");
        }
        r(k, 0)=Lhs(j); // R always has to be positive.
        /* logR[k]=log(Lhs(j)); */
        /* Rinv[k]=1/Lhs(j); */
        B(k, 0)=_prod(2, 2.0, 1.0/RxODE_safe_zero(Lhs(j)));
        for (j=neta+2; j < nLhs(); j++){
          /* Rprintf("j: %d; Adj: %d; k: %d\n",j, j-neta-2,k); */
          rp(k,j-neta-2) = Lhs(j);
          cur = as<mat>(c[j-neta-2]);
          cur(k,0) = _prod(2,Lhs(j),1.0/RxODE_safe_zero(r(k, 0)));
          c[j-neta-2] = cur;
        }
        for (j = 0; j < neta; j++){
          // .5*apply(eps*fp*B + .5*eps^2*B*c - c, 2, sum) - OMGAinv %*% ETA
          cur = as<mat>(c[j]);
	  lpi[k         +3*nObs()*j] = _prod(5, 0.25, err(k, 0), err(k, 0), B(k, 0), cur(k,0));
	  lpi[k+nObs()  +3*nObs()*j] = - _prod(2, 0.5, cur(k,0));
	  lpi[k+2*nObs()+3*nObs()*j] = - _prod(4, 0.5, err(k, 0), fpm(k, j), B(k, 0));
        }
        llik0[k] = -_prod(4,0.5,err(k, 0),err(k, 0),1.0/RxODE_safe_zero(r(k, 0)));
        llik0[k+nObs()] = -_prod(2, 0.5, RxODE_safe_log(r(k, 0)));
        k++;
      }
    }
    // Deferred sums (to reduce round-off error)
    for (j = 0; j < neta; j++){
      // Rprintf("%f", eta[j]);
      lp[j] = RxODE_sum(lpi+3*nObs()*j,3*nObs());
    }
    llik[0] = RxODE_sum(llik0,2*nObs());
    Free(lpi);
    Free(llik0);
    // Free
    RxODE_ode_free();
    mat llikm = mat(1,1);
    
    mat omegaInv = as<mat>(e["omegaInv"]);

    NumericVector llik2(1);
    llikm(0, 0)=  llik[0];
    llikm = -(llikm - 0.5*(etam.t() * omegaInv * etam));
    llik2[0] = llikm(0, 0);

    mat ep2 = -(lp - omegaInv * etam);
    if (do_table){
      mat omega = as<mat>(e["omega"]);
      mat Vfo_full = (fpm * omega * fpm.t()); // From Mentre 2006 p. 352
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
      // The Vfo calculation is separated out so it can be used in CWRES and WRES.
    
      mat Vfo = Vfo_full.diag();
      mat dErr_dEta = fpm * etam;
      e["Vfo"] = wrap(Vfo);
      e["dErr_dEta"] = wrap(dErr_dEta);
    }
    
    // Assign in env
    e["err"] = err;
    e["f"] = f;
    e["dErr"] = fpm;
    e["dR"] = rp;
    e["c"] = c;
    e["R"] = r;
    e["B"] = B;
    e["a"] = a;
    e["llik"] = llik;
    e["lp"] = lp;
    e["llik2"] = wrap(llik2);
    e["ep2"] = wrap(ep2);
  }
}

//' Get the Hessian for the environment
//' @param rho environment
//' @export
//' @keywords internal
// [[Rcpp::export]]
void rxHessian(SEXP rho){
  Environment e = as<Environment>(rho);
  int do_nonmem = as<int>(e["nonmem"]);
  int neta = as<int>(e["neta"]);
  mat omegaInv = as<mat>(e["omegaInv"]);
  mat B = as<mat>(e["B"]);
  List c = as<List>(e["c"]);
  List a = as<List>(e["a"]);
  NumericVector f = as<NumericVector>(e["f"]);
  mat H(neta, neta);
  int k, l, j = 0;
  double *sm = NULL;
  double tmp;
  mat al, ak, cl, ck;
  for (k = 0; k < neta; k++){
    for (l = 0; l <= k; l++){
      al = as<mat>(a[l]);
      ak = as<mat>(a[k]);
      cl = as<mat>(c[l]);
      ck = as<mat>(c[k]);
      if (sm == NULL){
        sm = Calloc(2*al.n_rows+1,double);
      }
      for (j = 0; j < (int)al.n_rows; j++){
	// products
	tmp = _prod(4, -0.5, al(j, 0), B(j, 0), ak(j, 0));
	if (ISNAN(tmp)){
	  tmp = 0;
	}
	sm[j] = tmp;
	
	tmp = (do_nonmem ? 1 : -1) * _prod(3, -0.5, cl(j, 0), ck(j, 0));
	if (ISNAN(tmp)){
	  tmp = 0;
	}
	sm[j+al.n_rows] = tmp;
      }
      sm[2*al.n_rows] = - omegaInv(k,l);
      // Deferred sums for accuracy
      H(k,l) = stablizeNums(RxODE_sum(sm, (int)(2*al.n_rows+1)));
      // Fill out the mirror compenent.
      if (l != k){
        H(l,k)=H(k,l);
      }
    }
  }
  if (sm != NULL){
    Free(sm);
  }
  // Rcout << "H:" << std::endl << H << std::endl <<
  //   "omegInv:" << std::endl << omegaInv << std::endl;
  e["H"] = H;
}

void rxInner2(SEXP sexp_eta, SEXP sexp_rho){
  Environment e = as<Environment>(sexp_rho);
  if (as<int>(e["numeric"]) == 1){
    rxInnerNum(sexp_eta, sexp_rho);
  } else {
    rxInner(sexp_eta, sexp_rho);
  }
  if (as<int>(e["switch.solver"]) == 1){
    if (as<int>(e["rc"]) != 0){
      if (as<int>(e["stiff"])){
        e["stiff"] = 0;
      } else {
        e["stiff"] = 1;
      }
      Rprintf("\nWarning: Switched solver.\n");
      e.remove("llik");
      if (as<int>(e["numeric"]) == 1){
        rxInnerNum(sexp_eta, sexp_rho);
      } else {
        rxInner(sexp_eta, sexp_rho);
      }
      // Switch back.
      if (as<int>(e["stiff"])){
        e["stiff"] = 0;
      } else {
        e["stiff"] = 1;
      }
    }
  }
}

//' Get the likelihood for ETA
//' @param sexp_eta ETA for likelihood
//' @param sexp_rho Environment with solving options
//' @export
//' @keywords internal
// [[Rcpp::export]]
NumericVector RxODE_focei_eta_lik(SEXP sexp_eta, SEXP sexp_rho){
  rxInner2(sexp_eta, sexp_rho);
  Environment e = as<Environment>(sexp_rho);
  NumericVector ret = as<NumericVector>(wrap(e["llik2"]));
  return ret;
}
//' Get the likelihood slope for ETA
//' @inheritParams RxODE_focei_eta_lik
//' @export
//' @keywords internal
// [[Rcpp::export]]
NumericVector RxODE_focei_eta_lp(SEXP sexp_eta, SEXP sexp_rho){
  rxInner2(sexp_eta, sexp_rho);
  Environment e = as<Environment>(sexp_rho);
  NumericVector ret = as<NumericVector>(wrap(e["ep2"]));
  return ret;
}
//' Get the Function pointers for the LBJ'S routine
//' @param fstr a string of the function pointer to return.  "lik" for likelihood and "lp" for likelihood gradient.
//' @export
//' @keywords internal
// [[Rcpp::export]]
XPtr<rxFn2> RxODE_focei_eta(std::string fstr){
  if (fstr == "lik")
    return(XPtr<rxFn2>(new rxFn2(&RxODE_focei_eta_lik)));
  else if (fstr == "lp")
    return(XPtr<rxFn2>(new rxFn2(&RxODE_focei_eta_lp)));
  else 
    return XPtr<rxFn2>(R_NilValue); // runtime error as NULL no XPtr
}
//' Finalize likelihood environment
//' @param  rho Environment to finalize.  Returns an individual likelihood. 
//' @export
//' @keywords internal
// [[Rcpp::export]]
NumericVector RxODE_focei_finalize_llik(SEXP rho){
  rxHessian(rho);
  Environment e = as<Environment>(rho);
  // Calculate -1/2 log(det(-H)) by chol.
  mat c;
  NumericVector reset;
  e["corrected"] = 0;
  try{
    c = chol(-as<mat>(e["H"]));
  } catch(...){
    e["reset"] = 0;
    c = -as<mat>(e["H"]);
    Function nearPD = as<Function>(e["nearPD"]);
    c = as<mat>(nearPD(c, rho));
    reset = as<NumericVector>(e["reset"]);
    if (reset[0] != 1){
      Rprintf("Warning: The Hessian is non-positive definite, correcting with nearPD\n");
      e["corrected"] = 1;
    } else {
      stop("Cannot correct Hessian Matrix to be non-positive definte matrix\n");
    }
    c = chol(c);
  }
  vec diag = c.diag();
  vec ldiag = log(diag);
  e["log.det.H.neg.5"]= wrap(sum(ldiag));
  NumericVector ret = -as<NumericVector>(e["llik2"])+ as<NumericVector>(e["log.det.OMGAinv.5"])-as<NumericVector>(e["log.det.H.neg.5"]);
  int do_table = as<int>(e["table"]);
  if (do_table){
    ret.attr("fitted") = as<NumericVector>(e["f"]);
    ret.attr("Vi") = as<NumericVector>(e["R"]); // Sigma
    ret.attr("Vfo") = as<NumericVector>(e["Vfo"]); // FO Variance (See Mentre 2006) Prediction Discrepancies for the Evaluation of Nonlinear Mixed-Effects Models
    ret.attr("dErr_dEta") = as<NumericVector>(e["dErr_dEta"]);
  }
  ret.attr("posthoc") = as<NumericVector>(e["eta"]);
  if (e.exists("c.hess")){
    ret.attr("c.hess") = as<NumericVector>(e["c.hess"]);
  }
  ret.attr("corrected") = as<NumericVector>(e["corrected"]);
  rxDetaDomega(e);
  if (!do_table){
    ret.attr("omega.28") = as<NumericVector>(e["omega.28"]);
  }
  // ret.attr("llik2") = as<NumericVector>(e["llik2"]);
  // ret.attr("log.det.OMGAinv.5") = as<NumericVector>(e["log.det.OMGAinv.5"]);
  // ret.attr("log.det.H.neg.5") = as<NumericVector>(e["log.det.H.neg.5"]);
  e["ret"] = ret;
  return ret;
}

//[[Rcpp::export]]
void rxDetaDomega(SEXP rho){
  // Used in  Eq #28
  Environment e = as<Environment>(rho);
  List dOmega = as<List>(e["omegaInv.dOmega.omegaInv"]);
  NumericVector  omegaInv = as<NumericVector>(e["tr.omegaInv.dOmega.0.5"]);
  List dOmega2 = as<List>(e["omegaInv.dOmega.omegaInv.dEta"]);
  mat eta = as<mat>(e["eta.mat"]);
  mat c,c2;
  vec ret;
  int ntheta = dOmega.length();
  NumericVector dEta(ntheta);
  mat o0 = as<mat>(dOmega[0]);
  int neta = o0.n_rows;
  int i,j;
  mat dEta47 = mat(neta,ntheta);
  for (i = 0; i < ntheta; i++){
    c = 0.5*(eta.t() * as<mat>(dOmega[i]) * eta);
    ret = c.diag();
    dEta[i] = sum(ret)-omegaInv[i];
    List dOmega2T = dOmega2[i];
    for (j = 0; j < neta; j++){
      c2 = eta.t() * as<mat>(dOmega2T[j]);
      dEta47(j, i) = c2(0, 0);
    }
  }
  e["omega.28"] = dEta;
  e["omega.47"] = dEta47;
}

// [[Rcpp::export]]
void rxOuter_ (SEXP rho){
  //Outer problem gradient for lbfgs
  Environment e = as<Environment>(rho);
  unsigned int i, j, k=0, h, n, i0 = 0,e1,e2, e3;
  
  mat omegaInv = as<mat>(e["omegaInv"]);
  
  unsigned int neta = as<unsigned int>(e["neta"]);
  unsigned int ntheta = as<unsigned int>(e["ntheta"]);
  List dOmega = as<List>(e["dOmega"]);
  NumericVector DV = as<NumericVector>(e["DV"]);
  int do_nonmem = as<int>(e["nonmem"]);
  
  unsigned int nomega = (unsigned int)(dOmega.size());

  RxODE_ode_solve_env(rho);
  
  mat fpm = mat(nObs(), neta);
  mat fpt = mat(nObs(),ntheta+nomega);
  List fp2(neta);
  
  mat rp = mat(nObs(),neta);
  mat rpt = mat(nObs(),ntheta+nomega);
  List rp2(neta);

  NumericVector f(nObs());
  mat err = mat(nObs(),1);
  mat r = mat(nObs(),1);

  mat B = mat(nObs(),1);
  List c(neta);
  List a(neta);
  
  List fpte(ntheta+nomega);
  List rpte(ntheta+nomega);

  NumericVector llik(1);
  mat lp = mat(neta,1);
  
  double *lpi =Calloc(neta*nObs()*3,double);
  double *llik0 =Calloc(nObs()*2,double);
  double *lDnDtS = Calloc(neta*ntheta*nObs()*8, double);
  double *lDnS = Calloc(neta*(neta+1)/2*nObs()*8, double);

  mat lDnDt = mat(neta,ntheta+nomega);
  mat lDn = mat(neta,neta);

  for (i = 0; i < neta; i++){
    a[i] = mat(nObs(),1);
    c[i] = mat(nObs(),1);
    fp2[i] = mat(nObs(),neta);
    rp2[i] = mat(nObs(),neta);
    lp(i,0) = 0;
    for (j = 0; j < neta; j++){
      lDn(i, j) = -omegaInv(i, j);
    }
  }
  for (j = 0; j < ntheta+nomega; j++){
    fpte[j] = mat(nObs(),neta);
    rpte[j] = mat(nObs(),neta);
  }

  llik[0]=0;
  
  // Now create the pred vector and d(pred)/d(eta) matrix.
  // Assuming Lhs(0) = pred and Lhs(1:n) = d(pred)/d(eta#)
  mat cur, cur2, cuR;
  for (i = 0; i < nAllTimes(); i++){
    if (!rxEvid(i)){
      rxCalcLhs(i);
      f[k] = Lhs(0); // Pred
      err(k, 0) = f[k] - DV[k];
      // d(pred)/d(eta#)
      // Rprintf("d(pred)/d(eta#)\n");
      for (j = 1; j < neta+1; j++){
        fpm(k,j-1) = Lhs(j);
        cur = as<mat>(a[j-1]);
        cur(k,0) =  Lhs(j);
        a[j-1]=cur;
      }
      /* // d(pred)/d(theta#) */
      // Rprintf("d(pred)/d(theta#)\n");
      i0 = 1+neta;
      for (j = i0; j < i0+ntheta; j++){
        fpt(k,j-i0) = Lhs(j);
      }
      /* // d^2(pred)/d^2(eta#) */
      // Rprintf("d^2(pred)/d^2(eta#)\n");
      i0 += ntheta;
      e1=0; e2=0;
      for (j = i0; j < i0+(neta)*(neta+1)/2; j++){
        /* fp2[(nAllTimes()-ixds)*(j-i0)+k] = Lhs(j); */
        cur = as<mat>(fp2[e1]);
        cur(k,e2) = Lhs(j);
        fp2[e1] = cur;
        if (e1 == e2){
          e1=0;
          e2++;
        } else {
          cur = as<mat>(fp2[e2]);
          cur(k,e1) = Lhs(j);
          fp2[e2]=cur;
          e1++;
        }
      }
      /* // d^2(pred)/(d(eta#)d(theta#)) */
      // Rprintf("d^2(pred)/d2(eta#)d(theta#)\n");
      i0 += neta*(neta+1)/2;
      h = 0;
      for (j = i0; j < i0 + neta*ntheta; ){
        for (n = 0; n < neta; n++){
          cur =as<mat>(fpte[h]);
          cur(k,n) =Lhs(j);
          fpte[h] =cur;
          j++;
        }
        h++;
      }
      i0 += neta*ntheta;
      j=i0;
      // Now
      if (Lhs(j) <= 0){
        Rprintf("R = Lhs(%d) = %f\n\n", j, Lhs(j));
        for (j = 0; j < nLhs(); j++){
          Rprintf("Lhs(%d) = %f\n", j, Lhs(j));
        }
        Rprintf("\n");
        RxODE_ode_free();
        stop("A covariance term is zero or negative and should remain positive.");
      }
      r(k, 0)=Lhs(j); // R always has to be positive.
      /* logR[k]=log(Lhs(j)); */
      /* Rinv[k]=1/Lhs(j); */
      B(k, 0)=2/RxODE_safe_zero(Lhs(j));
      /* // d(R)/d(eta#) */
      // Rprintf("d(R)/d(eta#)\n");
      i0++;
      for (j=i0; j < i0+neta; j++){
        /* Rprintf("j: %d; Adj: %d; k: %d\n",j, j-neta-2,k); */
        rp(k, j-i0) = Lhs(j);
        cur = as<mat>(c[j-i0]);
        // cur(k,0) = RxODE_sign_exp(Lhs(j)*RxODE_safe_zero(r(k, 0)),RxODE_abs_log(Lhs(j))-RxODE_abs_log(RxODE_safe_zero(r(k, 0))));
        cur(k, 0) = _prod(2, Lhs(j), 1.0/RxODE_safe_zero(r(k, 0)));
        c[j-i0] = cur;
        if (!do_nonmem){
          // tmp1[["_sens_rx_pred__ETA_1_"]],ncol=1) - err/R*matrix(tmp1[["_sens_rx_r__ETA_1_"]]
          cur =as<mat>(a[j-i0]);
          cur(k, 0)+= -_prod(3, err(k, 0), 1.0/RxODE_safe_zero(r(k, 0)), Lhs(j));
          a[j-i0] = cur;
        }
      }
      i0 += neta;
      /* // d(R)/d(theta#) */
      // Rprintf("d(R)/d(theta#)\n");
      for (j=i0; j < i0+ntheta; j++){
        /* Rprintf("j: %d; Adj: %d; k: %d\n",j, j-neta-2,k); */
        rpt(k, j-i0) = Lhs(j);
      }
      i0 += ntheta;
      /* // d^2(R)/d^2(eta) */
      // Rprintf("d(R)/d^2(eta#)\n");
      e1=0; e2=0;
      for (j=i0; j < i0+neta*(neta+1)/2; j++){
        /* Rprintf("j: %d; Adj: %d; k: %d\n",j, j-neta-2,k); */
        cur = as<mat>(rp2[e1]);
        cur(k,e2) = Lhs(j);
        rp2[e1] = cur;
        if (e1 == e2){
          e1=0;
          e2++;
        } else {
          cur = as<mat>(rp2[e2]);
          cur(k,e1) = Lhs(j);
          rp2[e2] = cur;
          e1++;
        }
      }
      // d^2(R)/(d(eta#)d(theta#))
      // Rprintf("d^2(R)/d(eta#)d(theta#)\n");
      i0 += neta*(neta+1)/2;
      h = 0;
      for (j = i0; j < i0+ntheta*neta; ){
        for (n = 0; n < neta; n++){
          cur = as<mat>(rpte[h]);
          cur(k,n) = Lhs(j);
          rpte[h] = cur;
          j++;
        }
        h++;
      }
      // Rprintf("lp\n");
      i0 += ntheta*neta;
      for (j = 0; j < neta; j++){
        //.5*apply(eps*fp*B + .5*eps^2*B*c - c, 2, sum) - OMGAinv %*% ETA
        // eq 12
        cur = as<mat>(c[j]);
        lpi[k         +3*nObs()*j] = _prod(5, 0.25, err(k, 0), err(k, 0), B(k, 0), cur(k,0));
        lpi[k+nObs()  +3*nObs()*j] = - _prod(2, 0.5, cur(k,0));
        lpi[k+2*nObs()+3*nObs()*j] = - _prod(4, 0.5, err(k, 0), fpm(k, j), B(k, 0));
      }
      // Rprintf("47\n");
      for (h=0; h < ntheta; h++){
        for (n = 0; n < neta; n++){
          // Eq #47 Almquist 2015
          // fpm = d(err)/d(eta)
          // fpt = d(err)/d(theta)
          // fp2 = d^2(err)/d(eta)^2
          // fpte = d^2(err)/d(eta)d(theta)
          // n = k
          // h = m
          cur = as<mat>(fpte[h]);
          cur2 = as<mat>(rpte[h]);
          
          // lDnDt(n, h) += -_prod(3, fpt(k, h), fpm(k, n), 1.0/RxODE_safe_zero(r(k, 0)))
          //   + _prod(5, err(k, 0), fpm(k, n), rpt(k, h), 1.0/RxODE_safe_zero(r(k, 0)), 1.0/RxODE_safe_zero(r(k, 0)))
          //   - _prod(3, err(k, 0), cur(k, n), 1.0/RxODE_safe_zero(r(k, 0)))
          //   + _prod(6, 0.5,err(k, 0), err(k, 0), cur2(k, n), 1.0/RxODE_safe_zero(r(k, 0)), 1.0/RxODE_safe_zero(r(k, 0)))
          //   - _prod(7, err(k, 0), err(k, 0), 1.0/RxODE_safe_zero(r(k, 0)), 1.0/RxODE_safe_zero(r(k, 0)), 1.0/RxODE_safe_zero(r(k, 0)), rp(k, n), rpt(k, h))
          //   + _prod(5, err(k, 0), rp(k, n), fpt(k, h), 1.0/RxODE_safe_zero(r(k, 0)), 1.0/RxODE_safe_zero(r(k, 0)))
          //   // trace is not needed since R is a scalar, not a vector
          //   - _prod(5, 0.5, rp(k, n), rpt(k, h), 1.0/RxODE_safe_zero(r(k, 0)), 1.0/RxODE_safe_zero(r(k, 0)))
          //   - _prod(3, 0.5, cur2(k, n), 1.0/RxODE_safe_zero(r(k, 0)))
          //   ;
          lDnDtS[k         +8*nObs()*(h+n*ntheta)] = -_prod(3, fpt(k, h), fpm(k, n), 1.0/RxODE_safe_zero(r(k, 0)));
          lDnDtS[k+nObs()  +8*nObs()*(h+n*ntheta)] = + _prod(5, err(k, 0), fpm(k, n), rpt(k, h), 1.0/RxODE_safe_zero(r(k, 0)), 1.0/RxODE_safe_zero(r(k, 0)));
          lDnDtS[k+2*nObs()+8*nObs()*(h+n*ntheta)] = - _prod(3, err(k, 0), cur(k, n), 1.0/RxODE_safe_zero(r(k, 0)));
          lDnDtS[k+3*nObs()+8*nObs()*(h+n*ntheta)] = + _prod(6, 0.5,err(k, 0), err(k, 0), cur2(k, n), 1.0/RxODE_safe_zero(r(k, 0)), 1.0/RxODE_safe_zero(r(k, 0)));
          lDnDtS[k+4*nObs()+8*nObs()*(h+n*ntheta)] = - _prod(7, err(k, 0), err(k, 0), 1.0/RxODE_safe_zero(r(k, 0)), 1.0/RxODE_safe_zero(r(k, 0)), 1.0/RxODE_safe_zero(r(k, 0)), rp(k, n), rpt(k, h));
          lDnDtS[k+5*nObs()+8*nObs()*(h+n*ntheta)] = + _prod(5, err(k, 0), rp(k, n), fpt(k, h), 1.0/RxODE_safe_zero(r(k, 0)), 1.0/RxODE_safe_zero(r(k, 0)));
          lDnDtS[k+6*nObs()+8*nObs()*(h+n*ntheta)] = - _prod(5, 0.5, rp(k, n), rpt(k, h), 1.0/RxODE_safe_zero(r(k, 0)), 1.0/RxODE_safe_zero(r(k, 0)));
          lDnDtS[k+7*nObs()+8*nObs()*(h+n*ntheta)] = - _prod(3, 0.5, cur2(k, n), 1.0/RxODE_safe_zero(r(k, 0)));
        }
      }
      // Rprintf("13\n");
      // Eq #13 Almquist 2015
      e3 = 0;
      for (e1 = 0; e1 < neta; e1++){
        for (e2 = 0; e2 <= e1; e2++){
          // fpm = d(err)/d(eta)
          // fpt = d(err)/d(theta)
          // fp2 = d^2(err)/d(eta)^2
          // fpte = d^2(err)/d(eta)d(theta)
          // e1 = k
          // e2 = l
          cur = as<mat>(fp2[e1]);
          cuR = as<mat>(rp2[e1]);
          lDnS[k         +8*nObs()*e3] = - _prod(3, fpm(k, e1), fpm(k, e2), 1.0/RxODE_safe_zero(r(k, 0)));
          lDnS[k+nObs()  +8*nObs()*e3] = + _prod(5, err(k, 0),  rp(k, e2), fpm(k,e1), 1.0/RxODE_safe_zero(r(k, 0)),  1.0/RxODE_safe_zero(r(k, 0)));
          lDnS[k+2*nObs()+8*nObs()*e3] = - _prod(3, err(k, 0),  cur(k,e2),  1.0/RxODE_safe_zero(r(k, 0)));
          lDnS[k+3*nObs()+8*nObs()*e3] = + _prod(6, 0.5, err(k, 0), err(k, 0), cuR(k, e2),  1.0/RxODE_safe_zero(r(k, 0)),  1.0/RxODE_safe_zero(r(k, 0)));
          lDnS[k+4*nObs()+8*nObs()*e3] = - _prod(7, err(k, 0), err(k, 0), rp(k, e1), rp(k, e2), 1.0/RxODE_safe_zero(r(k, 0)), 1.0/RxODE_safe_zero(r(k, 0)), 1.0/RxODE_safe_zero(r(k, 0)));
          lDnS[k+5*nObs()+8*nObs()*e3] = + _prod(5, err(k, 0), rp(k, e1), fpm(k, e2), 1.0/RxODE_safe_zero(r(k, 0)), 1.0/RxODE_safe_zero(r(k, 0)));
          lDnS[k+6*nObs()+8*nObs()*e3] = + _prod(5, 0.5, rp(k, e1), rp(k, e2), 1.0/RxODE_safe_zero(r(k, 0)), 1.0/RxODE_safe_zero(r(k, 0)));
          lDnS[k+7*nObs()+8*nObs()*e3] = - _prod(3, 0.5, cuR(k, e2), 1.0/RxODE_safe_zero(r(k, 0)));
          e3++;
          // lDn(e1, e2) += -_prod(3, fpm(k, e1), fpm(k, e2), 1.0/RxODE_safe_zero(r(k, 0)));
          // lDn(e1, e2) += + _prod(5, err(k, 0),  rp(k, e2), fpm(k,e1), 1.0/RxODE_safe_zero(r(k, 0)),  1.0/RxODE_safe_zero(r(k, 0)));
          // lDn(e1, e2) += - _prod(3, err(k, 0),  cur(k,e2),  1.0/RxODE_safe_zero(r(k, 0)));
          // lDn(e1, e2) += + _prod(6, 0.5, err(k, 0), err(k, 0), cuR(k, e2),  1.0/RxODE_safe_zero(r(k, 0)),  1.0/RxODE_safe_zero(r(k, 0)));
          // lDn(e1, e2) += - _prod(7, err(k, 0), err(k, 0), rp(k, e1), rp(k, e2), 1.0/RxODE_safe_zero(r(k, 0)), 1.0/RxODE_safe_zero(r(k, 0)), 1.0/RxODE_safe_zero(r(k, 0)));
          // lDn(e1, e2) += + _prod(5, err(k, 0), rp(k, e1), fpm(k, e2), 1.0/RxODE_safe_zero(r(k, 0)), 1.0/RxODE_safe_zero(r(k, 0)));
          // lDn(e1, e2) += + _prod(5, 0.5, rp(k, e1), rp(k, e2), 1.0/RxODE_safe_zero(r(k, 0)), 1.0/RxODE_safe_zero(r(k, 0)));
          // lDn(e1, e2) += - _prod(3, 0.5, cuR(k, e2), 1.0/RxODE_safe_zero(r(k, 0)));
          // lDn(e1, e2) += -(// d(err)/d*(e1+1)/2(eta,l)*R^-1*d(err)/d(eta,k)
          //               RxODE_sign_exp(fpm(k, e1) * fpm(k, e2) * RxODE_safe_zero(r(k, 0)),
          //                              RxODE_abs_log(fpm(k, e1)) +
          //                              RxODE_abs_log(fpm(k, e2)) -
          //                              RxODE_abs_log(RxODE_safe_zero(r(k, 0))))-
          //               // err/R^2*dR/d(eta,l)*derr/d(eta,k)
          //                  RxODE_sign_exp(err(k, 0) *  rp(k, e2) * fpm(k,e1),
          //                              RxODE_abs_log(err(k, 0)) +
          //                              RxODE_abs_log(rp(k, e2)) + // dR/d(eta,l)
          //                              RxODE_abs_log(fpm(k,e1)) - // derr/d(eta,k)
          //                              2*RxODE_abs_log(RxODE_safe_zero(r(k, 0))))+
          //               // err/R*derr^2/(d(eta,k)*d(eta,j))
          //                  RxODE_sign_exp(err(k, 0) * cur(k,e2) * RxODE_safe_zero(r(k, 0)),
          //                              RxODE_abs_log(err(k, 0)) +
          //                              RxODE_abs_log(cur(k,e2)) - //derr^2/deta(k,l)
          //                              RxODE_abs_log(RxODE_safe_zero(r(k, 0))))-
          //               // 0.5*err^2/R^2*dR/(deta(k,j))
          //               RxODE_sign_exp(cuR(k, e2),
          //                              log(0.5)+2*RxODE_abs_log(err(k, 0)) +
          //                              RxODE_abs_log(cuR(k, e2)) - 2*RxODE_abs_log(RxODE_safe_zero(r(k, 0)))) +
          //               // err^2/R^3*dR/(deta(k))*dR/(deta(l))
          //               RxODE_sign_exp(rp(k, e1) * rp(k, e2) * RxODE_safe_zero(r(k, 0)),
          //                              2*RxODE_abs_log(err(k, 0)) + RxODE_abs_log(rp(k, e1)) +
          //                              RxODE_abs_log(rp(k, e2)) -3* RxODE_abs_log(RxODE_safe_zero(r(k, 0)))) -
          //               // err/R^2*dR/(deta(k))*derr/deta(l)
          //                  RxODE_sign_exp(err(k, 0) * rp(k, e1) * fpm(k, e2),
          //                              RxODE_abs_log(err(k, 0)) + RxODE_abs_log(rp(k, e1)) + RxODE_abs_log(fpm(k, e2)) -
          //                              2*RxODE_abs_log(r(k, 0))) -
          //               // 0.5*trace(dR/deta(k)*dR/deta(l)/R^2)
          //                  RxODE_sign_exp(rp(k, e1) * rp(k, e2),
          //                              log(0.5)+RxODE_abs_log(rp(k, e1)) +
          //                              RxODE_abs_log(rp(k, e2)) - 2*RxODE_abs_log(RxODE_safe_zero(r(k, 0)))) +
          //               //0.5*trace(dR^2/d(eta(k,l))*R^-1)
          //               RxODE_sign_exp(cuR(k, e2)*RxODE_safe_zero(r(k, 0)),
          //                              log(0.5)+RxODE_abs_log(cuR(k, e2))- RxODE_abs_log(RxODE_safe_zero(r(k, 0)))));
          
        }
      }
      llik0[k] = -_prod(4,0.5,err(k, 0),err(k, 0),1.0/RxODE_safe_zero(r(k, 0)));
      llik0[k+nObs()] = -_prod(2, 0.5, RxODE_safe_log(r(k, 0)));
      k++;
    }
  }
  // Deferred sums (to reduce round-off error)
  for (j = 0; j < neta; j++){
    lp[j] = RxODE_sum(lpi+3*nObs()*j,3*nObs());
  }
  for (h=0; h < ntheta; h++){
    for (n = 0; n < neta; n++){
      lDnDt(n, h) = RxODE_sum(lDnDtS+8*nObs()*(h+n*ntheta), 8*nObs());
    }
  }
  e3 = 0;
  for (e1 = 0; e1 < neta; e1++){
    for (e2 = 0; e2 <= e1; e2++){
      // FIXME add omegaInf to RxODE_sum
      lDn(e1, e2) += RxODE_sum(lDnS+8*nObs()*e3, 8*nObs());
      lDn(e2, e1) = lDn(e1, e2);
      e3++;
    }
  }
  llik[0] = RxODE_sum(llik0,2*nObs());
  Free(lpi);
  Free(llik0);
  Free(lDnDtS);
  Free(lDnS);
  /* Finalize Eq #47 in Almquist 2015*/
  mat omega47 = as<mat>(e["omega.47"]);
  for (h=ntheta; h < ntheta+nomega; h++){
    for (n = 0; n < neta; n++){
      lDnDt(n, h) += -omega47(n,h-ntheta);
      // Finalize Eta2 and R2.
      for (i = 0; i < nObs(); i++){
        cur = as<mat>(fpte[h]);
        cur(i,n) = 0;
        fpte[h] = cur;
        cur = as<mat>(rpte[h]);
        cur(i,n) = 0;
        rpte[h] = cur;
      }
    }
    // Finalize dErr.dTheta to contain 0 for omega terms.
    for (i = 0; i < nObs(); i++){
      fpt(i, h) = 0;
      rpt(i, h) = 0;
    }
  }
  for (e1 = 0; e1 < neta; e1++){
    for (e2 = 0; e2 <= e1; e2++){
      lDn(e2, e1) = lDn(e1, e2);
    }
  }
  /* llik = -.5*sum(eps^2/(f^2*sig2) + log(f^2*sig2)) - .5*t(ETA) %*% OMGAinv %*% ETA */

  int do_table = as<int>(e["table"]);
  if (do_table){
    mat omega = as<mat>(e["omega"]);
    mat Vfo_full = (fpm * omega * fpm.t()); // From Mentre 2006 p. 352
    mat Vfo = Vfo_full.diag();
    mat etam = as<mat>(e["eta.mat"]);
    mat dErr_dEta = fpm * etam;
    mat dR_dEta = rp * etam;
    e["Vfo"] = wrap(Vfo);
    e["dErr_dEta"] = wrap(dErr_dEta);
  }
  
  e["f"] = f;
  e["err"] = err;
  
  e["dErr"] = fpm;
  e["dErr2"] = fp2;
  e["dErr.dTheta"] = fpt;
  e["dErr.dEta.dTheta"] = fpte;

  e["R"] = r;
  
  e["dR"] = rp;
  e["dR.dTheta"] = rpt;
  e["dR2"] = rp2;
  e["dR.dEta.dTheta"] = rpte;

  e["a"] = a;
  e["B"] = B;
  e["c"] = c;

  e["llik"] = llik;
  e["lp"] = lp;

  e["l.dEta.dTheta"] = lDnDt;
  e["H2"] = lDn;
  RxODE_ode_free();
}

// [[Rcpp::export]]
void rxDetaDtheta(SEXP rho){
  int i,h,n;
  Environment e = as<Environment>(rho);
  rxHessian(e); // Calculate Hessian
  mat H2 = as<mat>(e["H2"]);
  // mat H2 = as<mat>(e["H"]); // Although the paper prescribes using H2, using H gives a more accuate gradient...
  // Now caluclate dl(eta)/dTheta (Eq 28) and add to the overall dl/dTheta
  mat lDnDt = as<mat>(e["l.dEta.dTheta"]);
  mat iH2;
  e["reset"] = 0;
  NumericVector reset = as<NumericVector>(e["reset"]);
  try{
    iH2 = inv(H2);
  } catch(...){
    Rprintf("Warning: Hessian (H) seems singular; Using pseudo-inverse\n");
    iH2 = pinv(H2);
  }
  // # 46
  mat DnDt = -iH2 * lDnDt;
  e["dEta.dTheta"] = DnDt;
  // Now  (dErr/dTheta)*; #33
  mat dErrdTheta= as<mat>(e["dErr.dTheta"]);
  mat dErr = as<mat>(e["dErr"]); // derr/deta
  int ntheta = dErrdTheta.n_cols;
  // matrix(tmp2$dErr.dTheta[,1]) + tmp2$dErr %*% matrix(tmp2$dEta.dTheta[,1])
  mat dErrdTheta_ = mat(dErrdTheta.n_rows,0);
  for (i = 0; i < ntheta; i++){
    mat cur = dErrdTheta.col(i)+dErr * DnDt.col(i);
    dErrdTheta_ = join_rows(dErrdTheta_,cur);
  }
  e["dErr.dTheta."] = dErrdTheta_;
  // Now (dR/dTheta)*
  mat dRdTheta= as<mat>(e["dR.dTheta"]);
  mat dR = as<mat>(e["dR"]);
  // matrix(tmp2$dR.dTheta[,1]) + tmp2$dR %*% matrix(tmp2$dEta.dTheta[,1])
  mat dRdTheta_ = mat(dRdTheta.n_rows,0);
  mat cur;
  for (i = 0; i < ntheta; i++){
    cur = dRdTheta.col(i)+dR * DnDt.col(i);
    dRdTheta_ = join_rows(dRdTheta_,cur);
  }
  e["dR.dTheta."] = dRdTheta_;
  // Now #37
  // tmp2$dErr.dEta.dTheta[[theta]][,eta]-sum(over eta1,matrix(rowSums(tmp2$dErr2[[eta]][,eta1]*tmp2$dEta.dTheta[eta1,theta])))
  // tmp2$dErr.dEta.dTheta[[2]][,2]-matrix(rowSums(tmp2$dErr2[[1]]*tmp2$dEta.dTheta[2,2]))
  int neta = DnDt.n_rows;
  List dErrdEtadTheta_(neta);
  List dErrdEtadTheta = as<List>(e["dErr.dEta.dTheta"]); // dErr.dEta.dTheta
  List dErr2 = as<List>(e["dErr2"]);
  mat mat0, mat1, mat2;
  for (i = 0 ; i < neta; i++){
    cur = mat(dErrdTheta.n_rows,0);
    for(h = 0; h < ntheta; h++){
      // dErr2 =d^2(err)/d(eta)^2
      cur = join_rows(cur, (as<mat>(dErrdEtadTheta[h])).col(i) +
                      (as<mat>(dErr2[i])) * DnDt.col(h));
    }
    dErrdEtadTheta_[i]=cur;
  }
  e["dErr.dEta.dTheta."] = dErrdEtadTheta_;
  // And #37 equavialent for dR.
  List dRdEtadTheta_(neta);
  List dRdEtadTheta = as<List>(e["dR.dEta.dTheta"]);
  List dR2 = as<List>(e["dR2"]);
  for (i = 0 ; i < neta; i++){
    cur = mat(dRdTheta.n_rows,0);
    for(h = 0; h < ntheta; h++){
      cur = join_rows(cur, (as<mat>(dRdEtadTheta[h])).col(i) +
                      (as<mat>(dR2[i])) * DnDt.col(h));
    }
    dRdEtadTheta_[i]=cur;
  }
  e["dR.dEta.dTheta."] = dRdEtadTheta_;
  // Now dc*/dTheta #32
  // as.matrix(tmp2$dR.dTheta.[,theta]) * tmp2$dR[,eta]/(tmp2$R*tmp2$R)+tmp2$dR.dEta.dTheta.[[eta]][,theta]/tmp2$R
  List DcDh(neta);
  mat R = as<mat>(e["R"]);
  for (n = 0; n < neta; n++){
    cur = mat(dRdTheta.n_rows,0);
    for (h = 0; h < ntheta; h++){
      mat1 = as<mat>(dRdEtadTheta_[n]);
      mat2 = -dRdTheta_.col(h)  % dR.col(n)/(R % R)+mat1.col(h)/R;
      cur = join_rows(cur, mat2);
    }
    DcDh[n]=cur;
  }
  e["dc.dTheta"] = DcDh;
  // Now dB*/dTheta #31
  mat DbDh = mat(dRdTheta.n_rows,0);
  for (h = 0; h < ntheta; h++){
    mat1 = -2*dRdTheta_.col(h)/(R % R);
    DbDh = join_rows(DbDh, mat1);
  }
  e["dB.dTheta"] = DbDh;
  // Now da*/dTheta #30
  // matrix(tmp2$dErr.dEta.dTheta.[[eta]][,theta]) + matrix(tmp2$dErr.dTheta.[,theta])*matrix(tmp2$dR[,eta])/tmp2$R + matrix(tmp2$err)*matrix(tmp2$dR.dTheta.[,theta])*matrix(tmp2$dR[,eta])/(tmp2$R*tmp2$R) - tmp2$err*matrix(tmp2$tmp2$dR.dEta.dTheta.[[eta]][,theta])/tmp2$R
  // matrix(tmp2$dErr.dEta.dTheta.[[1]][,1]) + matrix(tmp2$dErr.dTheta.[,1]) * matrix(tmp2$dR[,1])/tmp2$R + matrix(tmp2$err)*matrix(tmp2$dR.dTheta.[,1])*matrix(tmp2$dR[,1])/(tmp2$R*tmp2$R) - tmp2$err*matrix(tmp2$dR.dEta.dTheta.[[1]][,1])/tmp2$R
  List DaDh(neta);
  mat mat3;
  mat err =as<mat>(e["err"]);
  NumericVector do_nonmem_v = as<NumericVector>(e["nonmem"]);
  int do_nonmem = (int)(do_nonmem_v[0]);
  for (n = 0; n < neta; n++){
    cur = mat(dRdTheta.n_rows,0);
    for (h = 0; h < ntheta; h++){
      mat1 = as<mat>(dRdEtadTheta_[n]);
      mat2 = as<mat>(dErrdEtadTheta_[n]);
      if (do_nonmem){
        mat3 = mat2.col(h);
      } else {
        mat3 = mat2.col(h) - dErrdTheta_.col(h) % dR.col(n) /R + err % dRdTheta_.col(h) % dR.col(n)/(R % R) - err % mat1.col(h)/R;
      }
      cur = join_rows(cur, mat3);
    }
    DaDh[n]=cur;
  }
  e["da.dTheta"] = DaDh;
  // Now calculate dH/dTheta (Eq 29)
  List DhDh(ntheta);
  int k, l;
  mat al, ak, dal, dak, cl, ck, dcl, dck;
  List a = as<List>(e["a"]);
  List c = as<List>(e["c"]);
  mat B = as<mat>(e["B"]);
  int ptheta = as<int>(e["ntheta"]);
  List dOmegainv = as<List>(e["dOmegaInv"]);
  for (h = 0; h < ntheta; h++){
    mat1 = mat(neta, neta);
    for (k = 0; k < neta; k++){
      for (l = 0; l <= k; l++){
        al = as<mat>(a[l]);
        ak = as<mat>(a[k]);
        dal = as<mat>(DaDh[l]);
        dak = as<mat>(DaDh[k]);
        cl = as<mat>(c[l]);
        ck = as<mat>(c[k]);
        dcl = as<mat>(DcDh[l]);
        dck = as<mat>(DcDh[k]);
        if (do_nonmem){
          mat1(k,l) = -0.5*sum(dal.col(h) % B % ak + al % DbDh.col(h) % ak + al % B % dak.col(h) + dcl.col(h) % ck + cl % dck.col(h));
        } else {
          mat1(k,l) = -0.5*sum(dal.col(h) % B % ak + al % DbDh.col(h) % ak + al % B % dak.col(h) - dcl.col(h) % ck - cl % dck.col(h));
        }
        if (h >= ptheta){
          // Put in dOmega^-1/dTheta term.
          mat2 = as<mat>(dOmegainv[h-ptheta]);
          mat1(k,l) = mat1(k,l)-mat2(k,l);
        }
        mat1(l,k) = mat1(k,l);
      }
    }
    DhDh[h]=mat1;
  }
  e["dH.dTheta"] = DhDh;
  mat H = as<mat>(e["H"]);
  mat Hinv;
  try{
    Hinv = inv(H);
  } catch(...){
    Rprintf("Warning: Hessian (H) seems singular; Using pseudo-inverse\n");
    Hinv = pinv(H);
  }
  e["Hinv"] = Hinv;
  NumericVector dEta = as<NumericVector>(e["omega.28"]);
  NumericVector dLdTheta(ntheta);
  for (h = 0; h < ntheta; h++){
    mat1 = -0.5*sum(2 * err % dErrdTheta.col(h) / R - err % err % dRdTheta.col(h) / (R % R) + dRdTheta.col(h) / R);
    if (h >= ptheta){
      mat1 = mat1+ dEta[h-ptheta];
    }
    // Now add -1/2*tr(Hinv*DhDh[h])
    mat2 = as<mat>(DhDh[h]);
    mat3 = Hinv * mat2;
    dLdTheta[h] = mat1(0,0)-0.5*sum(mat3.diag());
  }
  e["l.dTheta"] = dLdTheta;
  mat omegaInv = as<mat>(e["omegaInv"]);
  mat eta = as<mat>(e["eta.mat"]);
  mat aret = -(as<vec>(e["llik"])-0.5*(eta.t() * omegaInv * eta));
  e["llik2"] = aret(0,0);
  e["corrected"] = 0;
  mat cH;
  vec diag;
  vec ldiag;
  NumericVector dH;
  try{
    cH = chol(-as<mat>(e["H"]));
    diag = cH.diag();
    ldiag = log(diag);
    dH = wrap(sum(ldiag));
  } catch(...){
    cH = -as<mat>(e["H"]);
    dH = wrap(det(-cH));
    if (dH[0] > 0){
      dH=0.5*log(dH);
      Rprintf("Warning: The Hessian chol() decomposition failed, but we can use the (slower) determinant instead\n");
    } else {
      Function nearPD = as<Function>(e["nearPD"]);
      cH = as<mat>(nearPD(cH, rho));
      reset = as<NumericVector>(e["reset"]);
      if (reset[0] != 1){
        cH = chol(cH);
        diag = cH.diag();
        ldiag = log(diag);
        dH = wrap(sum(ldiag));
        Rprintf("Warning: The Hessian is non-positive definite, correcting with nearPD\n");
        e["corrected"] = 1;
      }
    }
  }
  if (reset[0] != 1){
    NumericVector ret(1);
    ret = -aret(0,0);
    // log(det(omegaInv^1/2)) = 1/2*log(det(omegaInv))
    ret += as<NumericVector>(e["log.det.OMGAinv.5"]);
    ret += -dH;
    ret.attr("fitted") = as<NumericVector>(e["f"]);
    mat etam = as<mat>(e["eta.mat"]);
    ret.attr("posthoc") = as<NumericVector>(wrap(e["eta.mat"]));
    if (e.exists("c.hess")){
      ret.attr("c.hess") = as<NumericVector>(e["c.hess"]);
    }
    if (e.exists("inits.vec")){
      // This calculation is done on the non-scaled parameters, but
      // needs to be changed to the scaled parameters.
      // This assumes scaling is to 1.0
      //
      NumericVector ini = as<NumericVector>(e["inits.vec"]);
      double scaleTo = as<double>(e["scale.to"]);
      if (ini.size() != dLdTheta.size()){
        stop("Inconsistent gradient and inits.vec size.");
      }
      NumericVector dLdThetaS(ntheta);
      // Even though the parameters are scaled for the overall
      // optimization, by the time the eta parameters are updated, the
      // parameters are unscaled. so dEta.dTheta are not scaled.
      for (h = 0; h < ntheta; h++){
	// f(theta(scaled))
	// The next is dscale/dTheta
	// L(theta(scaled));
	// We have dL/dtheta
	// theta(scaled) = scaled*ini/scaleTo or
	// scaled(theta) = theta/(ini/scaleTo)
	// dscaled / dtheta = 1/ini/scaleTo
        dLdThetaS[h] = dLdTheta[h]/(ini[h]/scaleTo);
      }
      e["l.dTheta.s"]=dLdThetaS;
      ret.attr("grad") = as<NumericVector>(wrap(dLdThetaS));
    } else {
      ret.attr("grad") = as<NumericVector>(wrap(dLdTheta));
    }
    ret.attr("dEta.dTheta") = as<NumericVector>(wrap(DnDt));
    mat omega = as<mat>(e["omega"]);
    mat fpm = as<mat>(e["dErr"]);
    int do_table = as<int>(e["table"]);
    if (do_table){
      mat Vfo_full = (fpm * omega * fpm.t()); // From Mentre 2006 p. 352
      mat Vfo = Vfo_full.diag();
      mat dErr_dEta = fpm * etam;
      ret.attr("Vfo") = wrap(Vfo);
      ret.attr("dErr_dEta") = wrap(dErr_dEta);
      ret.attr("Vi") = as<NumericVector>(e["R"]); // Sigma
    }
    ret.attr("corrected") = as<NumericVector>(e["corrected"]);
    e["ret"] = ret;
  } else {
    e["ret"] = NA_REAL;
  }
}

// [[Rcpp::export]]
NumericVector rxOuter(SEXP rho){
  rxDetaDomega(rho); // setup omega.28 and omega.47
  rxOuter_(rho);
  rxDetaDtheta(rho);
  Environment e = as<Environment>(rho);
  NumericVector ret = as<NumericVector>(e["ret"]);
  return ret;
}
//' Update ETAs based on d(eta)/d(theta)
//'
//' This updates the ETA initial estimates based on the knowledge of d(eta)/d(theta)
//'
//' @param DnDhS This is the d(eta)/d(theta) list where there is a
//'   d(eta)/d(theta) matrix for each subject
//'
//' @param DhS This is the change in theta observed between steps.
//'
//' @param initS This is the ETA initial condition matrix
//'
//' @param acceptNS Acceptance criteria for the new eta.  |eta| < acceptNS for the new eta
//'   to be accepted.
//'
//' @keywords internal
//'
//' @export
// [[Rcpp::export]]
NumericVector rxUpdateEtas(SEXP DnDhS, SEXP DhS, SEXP initS, SEXP acceptNS){
  int i = 0;
  uword j = 0;
  List DnDh = as<List>(DnDhS); // e["dEta.dTheta"]
  mat Dh = as<mat>(DhS); // e["dTheta"]
  mat inits = as<mat>(initS); // e["inits.mat"]
  mat cur;
  mat prod;
  double acceptN = as<double>(acceptNS);
  int accept;
  for (i = 0; i < DnDh.size(); i++){
    cur = as<mat>(DnDh[i]);
    prod = cur * Dh;
    accept = 1;
    for (j = 0; j < prod.n_rows; j++){
      if (abs(prod(j,0)+inits(i, j)) > acceptN){
        accept = 0;
        break;
      }
    }
    if (accept){
      for (j = 0; j < prod.n_rows; j++){
        inits(i, j) += prod(j, 0);
      }
    }
  }
  NumericVector ret = as<NumericVector>(wrap(inits));
  return ret;
}

// [[Rcpp::export]]
List foceiDataSetup(const DataFrame &df){
  // Purpose: get positions of each id and the length of each id's observations
  IntegerVector id    = df["ID"];
  IntegerVector evid  = df["EVID"];
  NumericVector dv    = df["DV"];
  NumericVector time0 = df["TIME"];
  NumericVector amt   = df["AMT"];
  int ids = id.size();
  int lastId = id[0]-1;
  // Get the number of subjects
  // Get the number of observations
  // Get the number of doses
  int nSub = 0, nObs = 0, nDoses = 0, i = 0, j = 0, k=0;
  for (i = 0; i < ids; i++){
    if (lastId != id[i]){
      nSub++;
      lastId=id[i];
    }
    if (evid[i]){
      nDoses++;
    } else {
      nObs++;
    }
  }
  // Now create data frames of observations and events
  NumericVector newDv(nObs);
  NumericVector newTimeO(nObs);
  
  IntegerVector newEvid(nDoses);
  NumericVector newAmt(nDoses);
  NumericVector newTimeA(nDoses);

  lastId = id[0]-1;
  IntegerVector newId(nSub);
  IntegerVector posDose(nSub);
  IntegerVector posObs(nSub);
  IntegerVector nDose(nSub);
  IntegerVector nObsN(nSub);
    
  int m = 0;
  for (i = 0; i < ids; i++){
    if (lastId != id[i]){
      lastId     = id[i];
      newId[m]   = id[i];
      posDose[m] = j;
      posObs[m]  = k;
      if (m != 0){
	nDose[m-1] = nDoses;
	nObsN[m-1]  = nObs;
      }
      nDoses = 0;
      nObs = 0;
      m++;
    }
    if (evid[i]){
      // Dose
      newEvid[j]  = evid[i];
      newTimeA[j] = time0[i];
      newAmt[j]   = amt[i];
      nDoses++;
      j++;
    } else {
      // Observation
      newDv[k]    = dv[i];
      newTimeO[k] = time0[i];
      nObs++;
      k++;
    }
  }
  nDose[m-1]=nDoses;
  nObsN[m-1]=nObs;
  return List::create(_["dose"]=DataFrame::create(_["evid"]   = newEvid,
						  _["time"]   = newTimeA,
						  _["amt"]    = newAmt),
		      _["obs"]=DataFrame::create(_["dv"]      = newDv,
						 _["time"]    = newTimeO),
		      _["ids"]=DataFrame::create(_["id"]      = newId,
						 _["posDose"] = posDose,
						 _["posObs"]  = posObs,
						 _["nDose"]   = nDose,
						 _["nObs"]    = nObsN));
}
