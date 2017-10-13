#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rconfig.h>
#include <R_ext/Rdynload.h>

typedef int (*rxode_int0)();
typedef unsigned int (*rxode_uint0)();
typedef int (*rxode_int1)(int i);
typedef double (*rxode_double1)(int i);
typedef void (*rxode_void0)();
typedef void (*rxode_void1)(int i);
typedef void (*rxode_void1s)(SEXP s);

typedef double (*rxode_double1d)(double x);
typedef double (*rxode_double2d)(double x, double y);

int nEq(){
  static rxode_int0 nEq = NULL;
  if (nEq == NULL) nEq = (rxode_int0) R_GetCCallable("RxODE","nEq");
  return nEq();
}

unsigned int nLhs(){
  static rxode_uint0 nLhs = NULL;
  if (nLhs == NULL) nLhs = (rxode_uint0) R_GetCCallable("RxODE","nLhs");
  return nLhs();
}

unsigned int nAllTimes(){
  static rxode_uint0 nAllTimes = NULL;
  if (nAllTimes == NULL) nAllTimes = (rxode_uint0) R_GetCCallable("RxODE","nAllTimes");
  return nAllTimes();
}

int rxEvid(int i){
  static rxode_int1 rxEvid = NULL;
  if (rxEvid == NULL) rxEvid = (rxode_int1) R_GetCCallable("RxODE","rxEvid");
  return rxEvid(i);
}

double rxLhs(int i){
  static rxode_double1 rxLhs = NULL;
  if (rxLhs == NULL) rxLhs = (rxode_double1) R_GetCCallable("RxODE","rxLhs");
  return rxLhs(i);
}

void rxCalcLhs(int i){
  static rxode_void1 rxCalcLhs = NULL;
  if (rxCalcLhs == NULL) rxCalcLhs = (rxode_void1) R_GetCCallable("RxODE","rxCalcLhs");
  rxCalcLhs(i);
}

unsigned int nObs(){
  static rxode_uint0 nObs = NULL;
  if (nObs == NULL) nObs = (rxode_uint0) R_GetCCallable("RxODE","nObs");
  return nObs();
}

void RxODE_ode_solve_env(SEXP sexp_rho){
  static rxode_void1s RxODE_ode_solve_env = NULL;
  if (RxODE_ode_solve_env == NULL) RxODE_ode_solve_env = (rxode_void1s) R_GetCCallable("RxODE","RxODE_ode_solve_env");
  RxODE_ode_solve_env(sexp_rho);
}

void RxODE_ode_free(){
  static rxode_void0 RxODE_ode_free = NULL;
  if (RxODE_ode_free == NULL) RxODE_ode_free = (rxode_void0) R_GetCCallable("RxODE","RxODE_ode_free");
  RxODE_ode_free();
}

double RxODE_safe_zero(double x){
  static rxode_double1d RxODE_safe_zero = NULL;
  if (RxODE_safe_zero == NULL) RxODE_safe_zero = (rxode_double1d) R_GetCCallable("RxODE","RxODE_safe_zero");
  return RxODE_safe_zero(x);
}

double RxODE_safe_log(double x){
  static rxode_double1d RxODE_safe_log = NULL;
  if (RxODE_safe_log == NULL) RxODE_safe_log = (rxode_double1d) R_GetCCallable("RxODE","RxODE_safe_log");
  return RxODE_safe_log(x);
}

double RxODE_abs_log(double x){
  static rxode_double1d RxODE_abs_log = NULL;
  if (RxODE_abs_log == NULL) RxODE_abs_log = (rxode_double1d) R_GetCCallable("RxODE","RxODE_abs_log");
  return RxODE_abs_log(x);
}

double RxODE_sign_exp(double sgn, double x){
  static rxode_double2d RxODE_sign_exp = NULL;
  if (RxODE_sign_exp == NULL) RxODE_sign_exp = (rxode_double2d) R_GetCCallable("RxODE","RxODE_sign_exp");
  return RxODE_sign_exp(sgn, x);
}

typedef double (*rxode_sum_fn)(double *input, int n);

extern double RxODE_sum(double *input, int n){
  static rxode_sum_fn RxODE_sum = NULL;
  if (RxODE_sum == NULL) RxODE_sum = (rxode_sum_fn) R_GetCCallable("RxODE","RxODE_sum");
  return RxODE_sum(input, n);
}

extern double RxODE_prod(double *input, int n){
  static rxode_sum_fn RxODE_prod = NULL;
  if (RxODE_prod == NULL) RxODE_prod = (rxode_sum_fn) R_GetCCallable("RxODE","RxODE_prod");
  return RxODE_prod(input, n);
}

extern double RxODE_sumV(int n, ...){
  va_list valist;
  va_start(valist, n);
  double *p = Calloc(n, double);
  for (unsigned int i = 0; i < n; i++){
    p[i] = va_arg(valist, double);
  }
  va_end(valist);
  double s = RxODE_sum(p, n);
  Free(p);
  return s;
}


extern double RxODE_prodV(int n, ...){
  va_list valist;
  va_start(valist, n);
  double *p = Calloc(n, double);
  for (int i = 0; i < n; i++){
    p[i] = va_arg(valist, double);
  }
  va_end(valist);
  double s = RxODE_prod(p, n);
  Free(p);
  return s;
}
