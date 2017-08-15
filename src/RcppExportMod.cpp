#include "evaluate.h"
#include "lin_cmt.h"

Rcpp::EvalBase *ev = NULL;                  // pointer to abstract base class
int NPAR=0;

typedef void (*S_fp) (double *, double *);
void nelder_(S_fp func, int n, double *start, double *step,
	     int itmax, double ftol_rel, double rcoef, double ecoef, double ccoef,
  int *iconv, int *it, int *nfcall, double *ynewlo, double *xmin,
  int *iprint);

void nmfn_wrap(double *x, double *fx)
{
	int i;
	Rcpp::NumericVector par(NPAR);
	for (i = 0; i < NPAR; i++) par[i] = x[i];
	*fx = as<double>(ev->eval(par));
}

RcppExport SEXP
neldermead_wrap(SEXP fcallSEXP, SEXP rhoSEXP, SEXP nparSEXP, SEXP startSEXP, SEXP stepSEXP,
                SEXP itmaxSEXP, SEXP ftol_relSEXP, SEXP rcoefSEXP, SEXP ecoefSEXP,
                SEXP ccoefSEXP,
                SEXP iprintSEXP)
{
BEGIN_RCPP

    ev = new Rcpp::EvalStandard(fcallSEXP, rhoSEXP);    // assign R function and environment
    int i, Iconv, Itnum, Nfcall, Itmax, iprint;
    double ftol_rel, rcoef, ecoef, ccoef;
    double Start[99], Xmin[99], Ynewlo, Step[99];

    NPAR = INTEGER(nparSEXP)[0];
    for (i=0; i<NPAR; i++) Start[i] = REAL(startSEXP)[i];
    for (i=0; i<NPAR; i++) Step[i]  = REAL(stepSEXP )[i];
    Itmax = INTEGER(itmaxSEXP)[0];
    ftol_rel = REAL(ftol_relSEXP)[0];
	rcoef = REAL(rcoefSEXP)[0];
	ecoef = REAL(ecoefSEXP)[0];
    ccoef = REAL(ccoefSEXP)[0];
    iprint = INTEGER(iprintSEXP)[0];

    nelder_(nmfn_wrap, NPAR, Start, Step,
            Itmax, ftol_rel, rcoef, ecoef, ccoef,
            &Iconv, &Itnum, &Nfcall, &Ynewlo, Xmin,
            &iprint);
  //Rcpp::Rcout <<ev->getNbEvals() <<std::endl;

 	Rcpp::NumericVector par(NPAR);
    for (i=0; i<NPAR; i++) par[i]  = Xmin[i];
	return Rcpp::List::create(Rcpp::Named("convergence") = Iconv,
	                          Rcpp::Named("Itnum") = Itnum,
	                          Rcpp::Named("iter") = Nfcall,
	                          Rcpp::Named("value") = Ynewlo,
	                          Rcpp::Named("par") = par);

END_RCPP
}


//------------
double uni_slice(double x0, double (*g)(double), double w, int m, double lower, double upper);

double slcfn_wrap(double x)
{
	Rcpp::NumericVector par(1);
	par[0] = x;
	return as<double>(ev->eval(par));
}

RcppExport SEXP
slice_wrap(SEXP fcallSEXP, SEXP rhoSEXP, SEXP x0SEXP, SEXP wSEXP, SEXP mSEXP, SEXP lowerSEXP, SEXP upperSEXP)
{
BEGIN_RCPP

    ev = new Rcpp::EvalStandard(fcallSEXP, rhoSEXP);    // assign R function and environment
    double x0=REAL(x0SEXP)[0], x1;
    double w=REAL(wSEXP)[0];
    int m=INTEGER(mSEXP)[0];
    double lower=REAL(lowerSEXP)[0];
    double upper=REAL(upperSEXP)[0];
    x1 = uni_slice(x0, slcfn_wrap, w, m, lower, upper);

	return Rcpp::List::create(Rcpp::Named("x1") = x1);

END_RCPP
}


RcppExport SEXP nlmixr_lin_cmt( SEXP obs_timeSEXP, SEXP dose_timeSEXP, SEXP doseSEXP, SEXP TinfSEXP,
	SEXP paramsSEXP, SEXP oralSEXP, SEXP infusionSEXP, SEXP ncmtSEXP, SEXP parameterizationSEXP ) {
BEGIN_RCPP

    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;

    using Eigen::VectorXd;
    Rcpp::traits::input_parameter< const VectorXd& >::type obs_time(obs_timeSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type dose_time(dose_timeSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type dose(doseSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type Tinf(TinfSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< const int >::type ncmt(ncmtSEXP);
    Rcpp::traits::input_parameter< const int >::type infusion(infusionSEXP);
    Rcpp::traits::input_parameter< const int >::type oral(oralSEXP);
    Rcpp::traits::input_parameter< const int >::type parameterization(parameterizationSEXP);

    VectorXd g;
    g = generic_cmt_interface<double>(
      obs_time,
      dose_time,
      dose,
      Tinf,
      params,
      ncmt,
      oral,
      infusion,
      parameterization);

    __result = Rcpp::wrap(g);
    return __result;

END_RCPP
}

Rcpp::EvalBase *fev = NULL;                  // pointer to abstract base class
Rcpp::EvalBase *gev = NULL;                  // pointer to abstract base class

typedef void (*S2_fp) (int *, int *, double *, double *, double *, int *, float *, double *);
extern "C" void n1qn1_ (S2_fp simul, int n[], double x[], double f[], double g[], double var[], double eps[],
                        int mode[], int niter[], int nsim[], int imp[], int lp[], double zm[], int izs[], float rzs[], double dzs[]);

unsigned int n1qn1_calls = 0;
int n1qn1_fprint = 0;
static void fwrap(int *ind, int *n, double *x, double *f, double *g, int *ti, float *tr, double *td)
{
  int i;
  Rcpp::NumericVector par(*n), ret(*n);
  for (i = 0; i < *n; i++) par[i] = x[i];
        
  if (*ind==2 || *ind==4) {
    n1qn1_calls++;
    ret = fev->eval(par);
    if (n1qn1_fprint){
      Rprintf("%3d:%#14.8g:", n1qn1_calls, ret[0]);
      for (i = 0; i < *n; i++) Rprintf(" %#8g", x[i]);
      Rprintf("\n");
    }
    *f = ret[0];
  }
  if (*ind==3 || *ind==4) {
    ret = gev->eval(par);
    for (i = 0; i < *n; i++) g[i] = ret[i];
  }
}


RcppExport SEXP
n1qn1_wrap(
           SEXP fSEXP, SEXP gSEXP, SEXP rhoSEXP, SEXP xSEXP, SEXP epsSEXP, 
           SEXP nSEXP, SEXP modeSEXP, SEXP niterSEXP, SEXP nsimSEXP, SEXP impSEXP,
           SEXP nzmSEXP, SEXP zmSEXP, SEXP fprint_sexp) {
  BEGIN_RCPP
    n1qn1_calls=0;
  n1qn1_fprint = as<int>(fprint_sexp);
  if (TYPEOF(fSEXP) == EXTPTRSXP){
    fev = new Rcpp::EvalCompiled(fSEXP, rhoSEXP); // xptr
  } else {
    fev = new Rcpp::EvalStandard(fSEXP, rhoSEXP); // Standard evaulation
  }
  if (TYPEOF(gSEXP) == EXTPTRSXP){
    gev = new Rcpp::EvalCompiled(gSEXP, rhoSEXP); // xptr
  } else {
    gev = new Rcpp::EvalStandard(gSEXP, rhoSEXP); // Standard evaulation
  }
  int i, j, k =0;
    
  int n, mode, niter, nsim, imp, lp=6, nzm;
  n = INTEGER(nSEXP)[0];
  mode = INTEGER(modeSEXP)[0];
  niter = INTEGER(niterSEXP)[0];
  nsim = INTEGER(nsimSEXP)[0];
  imp = INTEGER(impSEXP)[0];
  nzm = INTEGER(nzmSEXP)[0];
  double x[n], f, g[n], var[n], eps, zm[nzm];
  int izs[1]; float rzs[1]; double dzs[1];
  for (i=0; i<n; i++) x[i] = REAL(xSEXP)[i];
  for (i=0; i<nzm; i++) zm[i] = REAL(zmSEXP)[i];
  eps = REAL(epsSEXP)[0];
  for (i=0; i<n; i++) var[i] = .1;

  n1qn1_(fwrap,&n,x,&f,g,var,&eps,
         &mode,&niter,&nsim,&imp,&lp,zm,izs,rzs,dzs);
        
  Rcpp::NumericVector par(n);
  for (i=0; i<n; i++) par[i] = x[i];
  Rcpp::NumericVector hess(nzm);
  // On input this is hessian
  // On output this is H = LDL'
  // Triangular matrix is paramterized by column instead of row.
  using namespace arma;
  mat L = mat(n,n);
  mat D = mat(n,n);
  mat H = mat(n,n);
  // NumericVector zms(nzm);
  // for (i = 0; i < nzm; i++) zms[i]=zm[i];
  L.zeros();
  D.zeros();
  k =0;
  for (i=0; i<n; i++){
    for (j=i; j<n; j++){
      if (i == j){
	D(i,i)=zm[k];
	L(i,i)=1;
      } else {
	L(j,i)=zm[k];
      }
      k++;
    }
  }
  H = L*D*L.t();
  k = 0;
  for (i=0; i<n; i++){
    for (j=i; j<n; j++){
      hess[k]=H(j,i);
      k++;
    }
  }
  return Rcpp::List::create(Rcpp::Named("value") = f,
                            Rcpp::Named("par") = par,
			    // Rcpp::Named("L") = L,
			    // Rcpp::Named("D") = D,
			    Rcpp::Named("H") = H,
			    // Rcpp::Named("zm")=zms,
			    Rcpp::Named("c.hess") = hess);
        
  END_RCPP
}

extern "C" void qnbd_ (int indqn,S2_fp simul,int n[], double x[],double f[], double g[],int imp[],int io[],double zero[],
		       int napmax[], int itmax[],double epsf[],double epsg[],double epsx[],double df0[],double binf[],double bsup[],
		       int nfac[], double trav[],int ntrav[],double itrav[],double nitrav[], int izs[], int rzs[], int dzs[]);

/*
inqn = indicator of qnbd

Notes from google translate from french.

Input: 1 = standard
2 =  dh and initialized at the beginning of work and Ifac, f, g initialize

output:
  If <0 incapable of calculating a point better than the initial point
  If = 0 stop by the user
  If> 0 one provides a better point than the starting point
  <-10 unsuitable input parameters
  = -6 stop when calculating the direction of descent and iter = 1
  = -5 stop when calculating the approximation of the hessian iter = 1
  = -3 Simul anomaly: negative sign at a point or F and g have been previously computed
  = -2 failure of linear search at the first iteration
  = -1 f not defined at the initial point
  = 1 stop on epsg
  = 2 stop on epsf
  = 3 stop on epsx
  = 4 napmax
  = 5 itmax
  = 6 slope in the opposite direction to the gradient too small
  = 7 stop when calculating the direction of descent
  = 8 stop when calculating the Hessian approximation
  = 10 stop by linear search failure, cause not specified
  = 11 idem with indsim <0
  = 12 a step too small close to a step too big 
    This may result from an error in the gradient
  = 13 too many calls in a linear search

n = Dimension of x
Binf, bsup terminals inf, sup, dim n e?
X = variables to optimize (control) es
F value of criterion s
G gradient of f s
Zero close zero machine e
Napmax maximum number of simulate calls
Itmax maximum number of descent itineraries e
Itrav vect work dim nitrav = 2n, decomposes into indic and izig
Nfac number of factorized variables (e if indqn = 2) s
Imp printing factor
    - Varies from 0 (no impressions) to 3 (many impressions)
Io number of the results file e
Epsx vect dim n precision on x e
Epsf critere stop on f e
Epsg stop if sup a norm2 (g +) / n e
Work vect work dim ntrav
  It is necessary to ntrav> n (n + 1) / 2 + 6n
  Df0> 0 decrement f prevue (take 1. by default) e
    Izs, rzs, dzs: cf modulopt standards

    Indications on internal variables a qnbd and zqnbd
    Izig is used for storing constraints (active if izag> 1)
    If i does not change d ens, remove 1 a izig (positive)
    Otherwise we add izag
    Factorization only if izig is zero
    Dh Hessian estimate dim n (n + 1) / 2 ranked in three parts
    Indic (i) new index of index i
    Indic vect dim n order of storage of indices
    No need to initialize it if indqn = 1
 */

/*

fr <- function(x) {   ## Rosenbrock Banana function
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
rho = environment(fr)
start=c(2,2); step = -.2*start; itmax=999
ftol_rel=1e-6; rcoef=1.; ecoef=2.; ccoef=.5

require(Rcpp)
dyn.load("nelder_wrap4.dll")
is.loaded("neldermead_wrap")
.Call("neldermead_wrap", fr, rho, length(start), start, step,
      as.integer(itmax), ftol_rel, rcoef, ecoef, ccoef)
dyn.unload("nelder_wrap4.dll")

*/

