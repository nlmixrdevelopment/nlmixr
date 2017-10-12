#include <stdio.h>
#include <stdarg.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#define JAC_Rprintf Rprintf
#define JAC0_Rprintf if (_jac_counter_val() == 0) Rprintf
#define ODE_Rprintf Rprintf
#define ODE0_Rprintf if (_dadt_counter_val() == 0) Rprintf
#define LHS_Rprintf Rprintf
#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))
#define R_pow Rx_pow
#define R_pow_di Rx_pow_di

// Types for par pointers.r
typedef void (*RxODE_update_par_ptr)(double t);
typedef double (*RxODE_transit3)(double t, double n, double mtt);
typedef double (*RxODE_fn) (double x);
typedef double (*RxODE_fn2) (double x, double y);
typedef double (*RxODE_fn2i) (double x, int i);
typedef double (*RxODE_transit4)(double t, double n, double mtt, double bio);
typedef double (*RxODE_vec) (int val);
typedef long (*RxODE_cnt) ();
typedef void (*RxODE_inc) ();
typedef double (*RxODE_val) ();
typedef SEXP (*RxODE_ode_solver) (SEXP sexp_theta, SEXP sexp_inits, SEXP sexp_lhs, SEXP sexp_time, SEXP sexp_evid,SEXP sexp_dose, SEXP sexp_pcov, SEXP sexp_cov, SEXP sexp_locf, SEXP sexp_atol, SEXP sexp_rtol, SEXP sexp_hmin, SEXP sexp_hmax, SEXP sexp_h0, SEXP sexp_mxordn, SEXP sexp_mxords, SEXP sexp_mx,SEXP sexp_stiff, SEXP sexp_transit_abs, SEXP sexp_object, SEXP sexp_extra_args, SEXP sexp_matrix, SEXP sexp_add_cov);
typedef void (*RxODE_assign_fn_pointers)(void (*fun_dydt)(unsigned int, double, double *, double *),void (*fun_calc_lhs)(double, double *, double *),void (*fun_calc_jac)(unsigned int, double, double *, double *, unsigned int),void (*fun_update_inis)(SEXP _ini_sexp),int fun_jt,int fun_mf, int fun_debug);

typedef void (*RxODE_ode_solver_old_c)(int *neq,double *theta,double *time,int *evid,int *ntime,double *inits,double *dose,double *ret,double *atol,double *rtol,int *stiff,int *transit_abs,int *nlhs,double *lhs,int *rc);
typedef void (*RxODE_ode_solver_0_6_c)(int *neq,double *theta,double *time,int *evid,int *ntime,double *inits,double *dose,double *ret,double *atol,double *rtol,int *stiff,int *transit_abs,int *nlhs,double *lhs,int *rc,double hmin, double hmax,double h0,int mxordn,int mxords,int mxstep);
typedef double (*RxODE_solveLinB)(double t, int linCmt, int diff1, int diff2, double A, double alpha, double B, double beta, double C, double gamma, double ka, double tlag);
typedef double (*RxODE_sum_prod)(double *input, int n);
// Give par pointers
RxODE_vec _par_ptr, _InfusionRate;
RxODE_update_par_ptr _update_par_ptr;
RxODE_cnt _dadt_counter_val, _jac_counter_val;
RxODE_inc _dadt_counter_inc, _jac_counter_inc;
RxODE_val podo, tlast;
RxODE_transit4 _transit4;
RxODE_transit3 _transit3;
RxODE_fn _safe_log, safe_zero, factorial, _as_zero, abs_log, abs_log1p;
RxODE_fn2 sign_exp, Rx_pow;
RxODE_fn2i Rx_pow_di;
RxODE_assign_fn_pointers _assign_fn_pointers;
RxODE_ode_solver_old_c _old_c;
RxODE_ode_solver_0_6_c _c_0_6;
RxODE_solveLinB solveLinB;
RxODE_sum_prod _sum1, _prod1;

double _sum(int n, ...){
  va_list valist;
  va_start(valist, n);
  double *p = Calloc(n, double);
  for (unsigned int i = 0; i < n; i++){
    p[i] = va_arg(valist, double);
  }
  va_end(valist);
  double s = _sum1(p, n);
  Free(p);
  return s;
}

double _prod(int n, ...){
  va_list valist;
  va_start(valist, n);
  double *p = Calloc(n, double);
  for (unsigned int i = 0; i < n; i++){
    p[i] = va_arg(valist, double);
  }
  va_end(valist);
  double s = _prod1(p, n);
  Free(p);
  return s;
}

extern double _sign(unsigned int n, ...){
  va_list valist;
  va_start(valist, n);
  double s = 1;
  for (unsigned int i = 0; i < n; i++){
    s = sign(va_arg(valist, double))*s;
    if (s == 0){
      break;
    }
  }
  va_end(valist);
  return s;
}

extern void rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver_ptr();



// prj-specific differential eqns
void rx_9bb7d7958abef0092fc2389e956d74b4_x64_dydt(unsigned int _neq, double t, double *__zzStateVar__, double *__DDtStateVar__)
{
double 
		depot,
	_ETA_1_,
	_THETA_1_,
	center,
	_ETA_2_,
	_THETA_2_,
	_ETA_3_,
	_THETA_3_,
	rx_pred_,
	rx_r_,
	_THETA_4_;

	(void)t;	(void)depot;
	(void)_ETA_1_;
	(void)_THETA_1_;
	(void)center;
	(void)_ETA_2_;
	(void)_THETA_2_;
	(void)_ETA_3_;
	(void)_THETA_3_;
	(void)rx_pred_;
	(void)rx_r_;
	(void)_THETA_4_;

	_update_par_ptr(t);
	_ETA_1_ = _par_ptr(0);
	_THETA_1_ = _par_ptr(1);
	_ETA_2_ = _par_ptr(2);
	_THETA_2_ = _par_ptr(3);
	_ETA_3_ = _par_ptr(4);
	_THETA_3_ = _par_ptr(5);
	_THETA_4_ = _par_ptr(6);

	depot = __zzStateVar__[0];
	center = __zzStateVar__[1];

	__DDtStateVar__[0]=_InfusionRate(0)+0.0+-depot*exp(_ETA_1_+_THETA_1_);
	__DDtStateVar__[1]=_InfusionRate(1)+0.0+-center*exp(_ETA_2_+_THETA_2_)*exp(-_ETA_3_-_THETA_3_)+depot*exp(_ETA_1_+_THETA_1_);
	rx_pred_=center*exp(-_ETA_3_-_THETA_3_);
	rx_r_=Rx_pow_di(_THETA_4_,2);
    _dadt_counter_inc();
}

// Jacobian derived vars
void rx_9bb7d7958abef0092fc2389e956d74b4_x64_calc_jac(unsigned int _neq, double t, double *__zzStateVar__, double *__PDStateVar__, unsigned int __NROWPD__) {
  _jac_counter_inc();
}
// Functional based initial conditions.
void rx_9bb7d7958abef0092fc2389e956d74b4_x64_inis(SEXP _ini_sexp){
	double *__zzStateVar__ = REAL(_ini_sexp);
	double t=0;
double 
		depot,
	_ETA_1_,
	_THETA_1_,
	center,
	_ETA_2_,
	_THETA_2_,
	_ETA_3_,
	_THETA_3_,
	rx_pred_,
	rx_r_,
	_THETA_4_;

	(void)t;	(void)depot;
	(void)_ETA_1_;
	(void)_THETA_1_;
	(void)center;
	(void)_ETA_2_;
	(void)_THETA_2_;
	(void)_ETA_3_;
	(void)_THETA_3_;
	(void)rx_pred_;
	(void)rx_r_;
	(void)_THETA_4_;

	_update_par_ptr(0.0);
	_ETA_1_ = _par_ptr(0);
	_THETA_1_ = _par_ptr(1);
	_ETA_2_ = _par_ptr(2);
	_THETA_2_ = _par_ptr(3);
	_ETA_3_ = _par_ptr(4);
	_THETA_3_ = _par_ptr(5);
	_THETA_4_ = _par_ptr(6);

	depot = __zzStateVar__[0];
	center = __zzStateVar__[1];

	rx_pred_=center*exp(-_ETA_3_-_THETA_3_);
	rx_r_=Rx_pow_di(_THETA_4_,2);
  __zzStateVar__[0]=depot;
  __zzStateVar__[1]=center;
}
// prj-specific derived vars
void rx_9bb7d7958abef0092fc2389e956d74b4_x64_calc_lhs(double t, double *__zzStateVar__, double *_lhs) {
double 
		__DDtStateVar_0__,
	__DDtStateVar_1__,
	depot,
	_ETA_1_,
	_THETA_1_,
	center,
	_ETA_2_,
	_THETA_2_,
	_ETA_3_,
	_THETA_3_,
	rx_pred_,
	rx_r_,
	_THETA_4_;

	(void)t;	(void)__DDtStateVar_0__;
	(void)__DDtStateVar_1__;
	(void)depot;
	(void)_ETA_1_;
	(void)_THETA_1_;
	(void)center;
	(void)_ETA_2_;
	(void)_THETA_2_;
	(void)_ETA_3_;
	(void)_THETA_3_;
	(void)rx_pred_;
	(void)rx_r_;
	(void)_THETA_4_;

	_update_par_ptr(t);
	_ETA_1_ = _par_ptr(0);
	_THETA_1_ = _par_ptr(1);
	_ETA_2_ = _par_ptr(2);
	_THETA_2_ = _par_ptr(3);
	_ETA_3_ = _par_ptr(4);
	_THETA_3_ = _par_ptr(5);
	_THETA_4_ = _par_ptr(6);

	depot = __zzStateVar__[0];
	center = __zzStateVar__[1];

	__DDtStateVar_0__=_InfusionRate(0)+0.0+-depot*exp(_ETA_1_+_THETA_1_);
	__DDtStateVar_1__=_InfusionRate(1)+0.0+-center*exp(_ETA_2_+_THETA_2_)*exp(-_ETA_3_-_THETA_3_)+depot*exp(_ETA_1_+_THETA_1_);
	rx_pred_=center*exp(-_ETA_3_-_THETA_3_);
	rx_r_=Rx_pow_di(_THETA_4_,2);

	_lhs[0]=rx_pred_;
	_lhs[1]=rx_r_;
}
extern SEXP rx_9bb7d7958abef0092fc2389e956d74b4_x64_model_vars(){
	SEXP lst      = PROTECT(allocVector(VECSXP, 12));
	SEXP names    = PROTECT(allocVector(STRSXP, 12));
	SEXP params   = PROTECT(allocVector(STRSXP, 7));
	SEXP lhs      = PROTECT(allocVector(STRSXP, 2));
	SEXP state    = PROTECT(allocVector(STRSXP, 2));
	SEXP stateRmS = PROTECT(allocVector(INTSXP, 2));
	int *stateRm  = INTEGER(stateRmS);
	SEXP sens     = PROTECT(allocVector(STRSXP, 0));
	SEXP fn_ini   = PROTECT(allocVector(STRSXP, 0));
	SEXP dfdy     = PROTECT(allocVector(STRSXP, 0));
	SEXP tran     = PROTECT(allocVector(STRSXP, 12));
	SEXP trann    = PROTECT(allocVector(STRSXP, 12));
	SEXP mmd5     = PROTECT(allocVector(STRSXP, 2));
	SEXP mmd5n    = PROTECT(allocVector(STRSXP, 2));
	SEXP model    = PROTECT(allocVector(STRSXP, 4));
	SEXP modeln   = PROTECT(allocVector(STRSXP, 4));
	SET_STRING_ELT(params,0,mkChar("ETA[1]"));
	SET_STRING_ELT(params,1,mkChar("THETA[1]"));
	SET_STRING_ELT(params,2,mkChar("ETA[2]"));
	SET_STRING_ELT(params,3,mkChar("THETA[2]"));
	SET_STRING_ELT(params,4,mkChar("ETA[3]"));
	SET_STRING_ELT(params,5,mkChar("THETA[3]"));
	SET_STRING_ELT(lhs,0,mkChar("rx_pred_"));
	SET_STRING_ELT(lhs,1,mkChar("rx_r_"));
	SET_STRING_ELT(params,6,mkChar("THETA[4]"));
	SET_STRING_ELT(state,0,mkChar("depot"));
	stateRm[0] = 0;
	SET_STRING_ELT(state,1,mkChar("center"));
	stateRm[1] = 0;
	SET_STRING_ELT(modeln,0,mkChar("model"));
	SET_STRING_ELT(model,0,mkChar("d/dt(depot)=rxRate(depot)+-depot*exp(ETA[1]+THETA[1]);\nd/dt(center)=rxRate(center)+-center*exp(ETA[2]+THETA[2])*exp(-ETA[3]-THETA[3])+depot*exp(ETA[1]+THETA[1]);\nrx_pred_=center*exp(-ETA[3]-THETA[3]);\nrx_r_=Rx_pow_di(THETA[4],2);\n"));
	SET_STRING_ELT(modeln,1,mkChar("normModel"));
	SET_STRING_ELT(model,1,mkChar("d/dt(depot)=rxRate(depot)+-depot*exp(ETA[1]+THETA[1]);\nd/dt(center)=rxRate(center)+-center*exp(ETA[2]+THETA[2])*exp(-ETA[3]-THETA[3])+depot*exp(ETA[1]+THETA[1]);\nrx_pred_=center*exp(-ETA[3]-THETA[3]);\nrx_r_=Rx_pow_di(THETA[4],2);\n"));
	SET_STRING_ELT(modeln,2,mkChar("parseModel"));
	SET_STRING_ELT(model,2,mkChar("__DDtStateVar__[0] = _InfusionRate(0) + 0.0+-depot*exp(_ETA_1_+_THETA_1_);\n__DDtStateVar__[1] = _InfusionRate(1) + 0.0+-center*exp(_ETA_2_+_THETA_2_)*exp(-_ETA_3_-_THETA_3_)+depot*exp(_ETA_1_+_THETA_1_);\nrx_pred_=center*exp(-_ETA_3_-_THETA_3_);\nrx_r_=Rx_pow_di(_THETA_4_,2);\n"));
	SET_STRING_ELT(modeln,3,mkChar("expandModel"));
	SET_STRING_ELT(model,3,mkChar("d/dt(depot)=rxRate(depot)+-depot*exp(ETA[1]+THETA[1]);\nd/dt(center)=rxRate(center)+-center*exp(ETA[2]+THETA[2])*exp(-ETA[3]-THETA[3])+depot*exp(ETA[1]+THETA[1]);\nrx_pred_=center*exp(-ETA[3]-THETA[3]);\nrx_r_=Rx_pow_di(THETA[4],2);\n"));
	SEXP ini    = PROTECT(allocVector(REALSXP,0));
	SEXP inin   = PROTECT(allocVector(STRSXP, 0));
	SET_STRING_ELT(names,0,mkChar("params"));
	SET_VECTOR_ELT(lst,  0,params);
	SET_STRING_ELT(names,1,mkChar("lhs"));
	SET_VECTOR_ELT(lst,  1,lhs);
	SET_STRING_ELT(names,2,mkChar("state"));
	SET_VECTOR_ELT(lst,  2,state);
	SET_STRING_ELT(names,3,mkChar("trans"));
	SET_VECTOR_ELT(lst,  3,tran);
	SET_STRING_ELT(names,5,mkChar("model"));
	SET_VECTOR_ELT(lst,  5,model);
	SET_STRING_ELT(names,4,mkChar("ini"));
	SET_VECTOR_ELT(lst,  4,ini);
	SET_STRING_ELT(names,6,mkChar("md5"));
	SET_VECTOR_ELT(lst,  6,mmd5);
	SET_STRING_ELT(names,7,mkChar("podo"));
	SET_VECTOR_ELT(lst,  7,ScalarLogical(0));
	SET_STRING_ELT(names,8,mkChar("dfdy"));
	SET_VECTOR_ELT(lst,  8,dfdy);
	SET_STRING_ELT(names,9,mkChar("sens"));
	SET_VECTOR_ELT(lst,  9,sens);
	SET_STRING_ELT(names,10,mkChar("fn.ini"));
	SET_VECTOR_ELT(lst,  10,fn_ini);
	SET_STRING_ELT(names,11,mkChar("state.ignore"));
	SET_VECTOR_ELT(lst,  11,stateRmS);
	SET_STRING_ELT(mmd5n,0,mkChar("file_md5"));
	SET_STRING_ELT(mmd5,0,mkChar("74a7df9b878113db4b69ffada8b0208b"));
	SET_STRING_ELT(mmd5n,1,mkChar("parsed_md5"));
	SET_STRING_ELT(mmd5,1,mkChar("9bb7d7958abef0092fc2389e956d74b4"));
	SET_STRING_ELT(trann,0,mkChar("jac"));
	SET_STRING_ELT(tran,0,mkChar("fullint"));
	SET_STRING_ELT(trann,1,mkChar("prefix"));
	SET_STRING_ELT(tran, 1,mkChar("rx_9bb7d7958abef0092fc2389e956d74b4_x64_"));
	SET_STRING_ELT(trann,2,mkChar("dydt"));
	SET_STRING_ELT(tran, 2,mkChar("rx_9bb7d7958abef0092fc2389e956d74b4_x64_dydt"));
	SET_STRING_ELT(trann,3,mkChar("calc_jac"));
	SET_STRING_ELT(tran, 3,mkChar("rx_9bb7d7958abef0092fc2389e956d74b4_x64_calc_jac"));
	SET_STRING_ELT(trann,4,mkChar("calc_lhs"));
	SET_STRING_ELT(tran, 4,mkChar("rx_9bb7d7958abef0092fc2389e956d74b4_x64_calc_lhs"));
	SET_STRING_ELT(trann,5,mkChar("model_vars"));
	SET_STRING_ELT(tran, 5,mkChar("rx_9bb7d7958abef0092fc2389e956d74b4_x64_model_vars"));
	SET_STRING_ELT(trann,6,mkChar("ode_solver"));
	SET_STRING_ELT(tran, 6,mkChar("rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver"));
	SET_STRING_ELT(trann,7,mkChar("ode_solver_sexp"));
	SET_STRING_ELT(tran, 7,mkChar("rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver_sexp"));
	SET_STRING_ELT(trann,8,mkChar("ode_solver_0_6"));
	SET_STRING_ELT(tran, 8,mkChar("rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver_0_6"));
	SET_STRING_ELT(trann,9,mkChar("ode_solver_focei_eta"));
	SET_STRING_ELT(tran, 9,mkChar("rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver_focei_eta"));
	SET_STRING_ELT(trann,10,mkChar("ode_solver_ptr"));
	SET_STRING_ELT(tran, 10,mkChar("rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver_ptr"));
	SET_STRING_ELT(trann,11,mkChar("inis"));
	SET_STRING_ELT(tran, 11,mkChar("rx_9bb7d7958abef0092fc2389e956d74b4_x64_inis"));
	setAttrib(tran, R_NamesSymbol, trann);
	setAttrib(mmd5, R_NamesSymbol, mmd5n);
	setAttrib(model, R_NamesSymbol, modeln);
	setAttrib(ini, R_NamesSymbol, inin);
	setAttrib(lst, R_NamesSymbol, names);
	UNPROTECT(17);
	return lst;
}
extern void rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver(
                    int *neq,
                    double *theta,      //order:
                    double *time,
                    int *evid,
                    int *ntime,
                    double *inits,
                    double *dose,
                    double *ret,
                    double *atol,
                    double *rtol,
                    int *stiff,
                    int *transit_abs,
                    int *nlhs,
                    double *lhs,
                    int *rc
                    ){
  // Backward compatible ode solver for 0.5* C interface
  rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver_ptr();
  _old_c(neq, theta, time, evid, ntime, inits, dose, ret, atol, rtol, stiff, transit_abs, nlhs, lhs, rc);
}

void rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver_0_6(int *neq,
                        double *theta,  //order:
                        double *time,
                        int *evid,
                        int *ntime,
                        double *inits,
                        double *dose,
                        double *ret,
                        double *atol,
                        double *rtol,
                        int *stiff,
                        int *transit_abs,
                        int *nlhs,
                        double *lhs,
                        int *rc,
                        double hmin,
                        double hmax,
                        double h0,
                        int mxordn,
                        int mxords,
                        int mxstep) {
  // Backward compatible ode solver for 0.5* C interface
  rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver_ptr();
  _c_0_6(neq, theta, time, evid, ntime, inits, dose, ret, atol, rtol, stiff, transit_abs, nlhs, lhs, rc,
	hmin, hmax, h0, mxordn, mxords, mxstep);
}

extern void rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver_ptr  (){
  _assign_fn_pointers(rx_9bb7d7958abef0092fc2389e956d74b4_x64_dydt , rx_9bb7d7958abef0092fc2389e956d74b4_x64_calc_lhs , rx_9bb7d7958abef0092fc2389e956d74b4_x64_calc_jac, rx_9bb7d7958abef0092fc2389e956d74b4_x64_inis, 2 , 22,
#ifdef __DEBUG__
                      1
#else
                      0
#endif
                      );
}

extern SEXP rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver_sexp (// Parameters
                                 SEXP sexp_theta,
                                 SEXP sexp_inits,
                                 SEXP sexp_lhs,
				 // Events
				 SEXP sexp_time,
				 SEXP sexp_evid,
				 SEXP sexp_dose,
				 // Covariates
				 SEXP sexp_pcov,
				 SEXP sexp_cov,
				 SEXP sexp_locf,
				 // Solver Options
				 SEXP sexp_atol,
				 SEXP sexp_rtol,
				 SEXP sexp_hmin,
				 SEXP sexp_hmax,
				 SEXP sexp_h0,
				 SEXP sexp_mxordn,
				 SEXP sexp_mxords,
				 SEXP sexp_mx,
				 SEXP sexp_stiff,
				 SEXP sexp_transit_abs,
				 // Object Creation
				 SEXP sexp_object,
				 SEXP sexp_extra_args,
				 SEXP sexp_matrix,
				 SEXP sexp_add_cov){
  rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver_ptr();
  RxODE_ode_solver ode_solver = (RxODE_ode_solver) R_GetCCallable("RxODE","RxODE_ode_solver");
  return ode_solver(sexp_theta,sexp_inits,sexp_lhs,sexp_time,sexp_evid,sexp_dose,sexp_pcov,sexp_cov,sexp_locf,sexp_atol,
		    sexp_rtol,sexp_hmin,sexp_hmax,sexp_h0,sexp_mxordn,sexp_mxords,sexp_mx,sexp_stiff,sexp_transit_abs,
		    sexp_object,sexp_extra_args,sexp_matrix,sexp_add_cov);
}

static R_NativePrimitiveArgType rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solverrx_t[] = {
  //*neq, *theta, *time,  *evid, *ntime, *inits,   *dose,   *ret,     *atol,  *rtol,   *stiff, *transit_abs, *nlhs, *lhs, *rc
  INTSXP,REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP
};

static R_NativePrimitiveArgType  rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver_0_6rx_t[] = {
  //*neq, *theta, *time,  *evid, *ntime, *inits,   *dose,   *ret,     *atol,  *rtol,   *stiff, *transit_abs, *nlhs, *lhs, *rc
  INTSXP,REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP
};

//Initilize the dll to match RxODE's calls
void R_init_rx_9bb7d7958abef0092fc2389e956d74b4_x64 (DllInfo *info){
  // Get the RxODE calling interfaces
  _InfusionRate   = (RxODE_vec) R_GetCCallable("RxODE","RxODE_InfusionRate");
  _update_par_ptr = (RxODE_update_par_ptr) R_GetCCallable("RxODE","RxODE_update_par_ptr");
  _par_ptr = (RxODE_vec) R_GetCCallable("RxODE","RxODE_par_ptr");
  _dadt_counter_val = (RxODE_cnt) R_GetCCallable("RxODE","RxODE_dadt_counter_val");
  _jac_counter_val  = (RxODE_cnt) R_GetCCallable("RxODE","RxODE_jac_counter_val");
  _dadt_counter_inc = (RxODE_inc) R_GetCCallable("RxODE","RxODE_dadt_counter_inc");
  _jac_counter_inc  = (RxODE_inc) R_GetCCallable("RxODE","RxODE_jac_counter_inc");
  podo  = (RxODE_val) R_GetCCallable("RxODE","RxODE_podo");
  tlast = (RxODE_val) R_GetCCallable("RxODE","RxODE_tlast");
  factorial=(RxODE_fn) R_GetCCallable("RxODE","RxODE_factorial");
  _transit3 = (RxODE_transit3) R_GetCCallable("RxODE","RxODE_transit3");
  _transit4 = (RxODE_transit4) R_GetCCallable("RxODE","RxODE_transit4");
  _safe_log =(RxODE_fn) R_GetCCallable("RxODE","RxODE_safe_log");
  safe_zero =(RxODE_fn) R_GetCCallable("RxODE","RxODE_safe_zero");
  _as_zero = (RxODE_fn) R_GetCCallable("RxODE","RxODE_as_zero");
  _assign_fn_pointers=(RxODE_assign_fn_pointers) R_GetCCallable("RxODE","RxODE_assign_fn_pointers");
  _old_c = (RxODE_ode_solver_old_c) R_GetCCallable("RxODE","RxODE_ode_solver_old_c");
  _c_0_6 = (RxODE_ode_solver_0_6_c)R_GetCCallable("RxODE","RxODE_ode_solver_0_6_c");
  _sum1   = (RxODE_sum_prod)R_GetCCallable("RxODE","RxODE_sum");
  _prod1 = (RxODE_sum_prod) R_GetCCallable("RxODE","RxODE_prod");
  sign_exp = (RxODE_fn2) R_GetCCallable("RxODE","RxODE_sign_exp");
  Rx_pow = (RxODE_fn2) R_GetCCallable("RxODE","RxODE_pow");
  Rx_pow_di = (RxODE_fn2i) R_GetCCallable("RxODE","RxODE_pow_di");
  abs_log = (RxODE_fn) R_GetCCallable("RxODE","RxODE_abs_log");
  abs_log1p = (RxODE_fn) R_GetCCallable("RxODE","RxODE_abs_log1p");
  solveLinB = (RxODE_solveLinB) R_GetCCallable("RxODE","RxODE_solveLinB");
  // Register the outside functions
  R_RegisterCCallable("rx_9bb7d7958abef0092fc2389e956d74b4_x64","rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver",       (DL_FUNC) rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver);
  R_RegisterCCallable("rx_9bb7d7958abef0092fc2389e956d74b4_x64","rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver_sexp",  (DL_FUNC) rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver_sexp);
  R_RegisterCCallable("rx_9bb7d7958abef0092fc2389e956d74b4_x64","rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver_0_6",   (DL_FUNC) rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver_0_6);
  R_RegisterCCallable("rx_9bb7d7958abef0092fc2389e956d74b4_x64","rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver_ptr",   (DL_FUNC) rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver_ptr);

  /* R_registerRoutines(info, NULL, NULL, NULL, NULL); */
  /* R_useDynamicSymbols(info,TRUE); */

  static const R_CMethodDef cMethods[] = {
    {"rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver", (DL_FUNC) &rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver, 15, rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solverrx_t},
    {"rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver_0_6", (DL_FUNC) &rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver_0_6, 21, rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver_0_6rx_t},
    {NULL, NULL, 0, NULL}
  };
  
  R_CallMethodDef callMethods[]  = {
    {"rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver_ptr", (DL_FUNC) &rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver_ptr, 0},
    {"rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver_sexp", (DL_FUNC) &rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver_sexp, 23},
    {"rx_9bb7d7958abef0092fc2389e956d74b4_x64_model_vars", (DL_FUNC) &rx_9bb7d7958abef0092fc2389e956d74b4_x64_model_vars, 0},
    {"rx_9bb7d7958abef0092fc2389e956d74b4_x64_inis", (DL_FUNC) &rx_9bb7d7958abef0092fc2389e956d74b4_x64_inis, 1},
    {NULL, NULL, 0}
  };
  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info,FALSE);
  // Register the function pointers so if someone directly calls the
  // ode solvers directly, they use the last loaded RxODE model.
  rx_9bb7d7958abef0092fc2389e956d74b4_x64_ode_solver_ptr();
}
