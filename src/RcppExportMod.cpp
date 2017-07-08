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
    *fx = ev->eval(par);
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
    return ev->eval(par);
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

