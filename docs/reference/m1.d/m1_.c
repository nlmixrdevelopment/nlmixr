#include <RxODE_model_shared.h>
#define __MAX_PROD__ 0
#define _CMT CMT
#define _SYNC_simeps_ for (int _svari=_solveData->neps; _svari--;){  if (_solveData->op->svar[_svari] == 0) {KA = _PP[0];};   if (_solveData->op->svar[_svari] == 1) {KE = _PP[1];}; }
#define _SYNC_simeta_ for (int _ovari=_solveData->neta; _ovari--;){  if (_solveData->op->ovar[_ovari] == 0) {KA = _PP[0];};   if (_solveData->op->ovar[_ovari] == 1) {KE = _PP[1];}; }
#include "extraC.h"
#include <RxODE_model_shared.c>
extern void  m1__ode_solver_solvedata (rx_solve *solve){
  _solveData = solve;
}
extern rx_solve *m1__ode_solver_get_solvedata(){
  return _solveData;
}
SEXP m1__model_vars();


// prj-specific differential eqns
void m1__dydt(int *_neq, double __t, double *__zzStateVar__, double *__DDtStateVar__)
{
  int _itwhile = 0;
  (void)_itwhile;
  int _cSub = _neq[1];
  double t = __t + _solveData->subjects[_neq[1]].curShift;
  (void)t;
    double depot;
  double KA;
  double centr;
  double KE;

  (void)t;
  (void)depot;
  (void)KA;
  (void)centr;
  (void)KE;


  _update_par_ptr(__t, _cSub, _solveData, _idx);
  KA = _PP[0];
  KE = _PP[1];

  depot = __zzStateVar__[0]*((double)(_ON[0]));
  centr = __zzStateVar__[1]*((double)(_ON[1]));

  __DDtStateVar__[0] = ((double)(_ON[0]))*(_IR[0] -KA*depot);
  __DDtStateVar__[1] = ((double)(_ON[1]))*(_IR[1] + KA*depot-KE*centr);
  (&_solveData->subjects[_cSub])->dadt_counter[0]++;
}

// Jacobian derived vars
void m1__calc_jac(int *_neq, double __t, double *__zzStateVar__, double *__PDStateVar__, unsigned int __NROWPD__) {
  int _itwhile = 0;
  (void)_itwhile;
    int _cSub=_neq[1];
  double t = __t + _solveData->subjects[_neq[1]].curShift;
  (void)t;
    (&_solveData->subjects[_cSub])->jac_counter[0]++;
}
// Functional based initial conditions.
void m1__inis(int _cSub, double *__zzStateVar__){
  int _itwhile = 0;
  (void)_itwhile;
  
}
// prj-specific derived vars
void m1__calc_lhs(int _cSub, double __t, double *__zzStateVar__, double *_lhs) {
    int _itwhile = 0;
  (void)_itwhile;
  double t = __t + _solveData->subjects[_cSub].curShift;
  (void)t;
  }
// Functional based bioavailability
double m1__F(int _cSub,  int _cmt, double _amt, double __t, double *__zzStateVar__){
 return _amt;
}
// Functional based absorption lag
double m1__Lag(int _cSub,  int _cmt, double __t, double *__zzStateVar__){
 return __t;
}
// Modeled zero-order rate
double m1__Rate(int _cSub,  int _cmt, double _amt, double __t, double *__zzStateVar__){
 return 0.0;
}
// Modeled zero-order duration
double m1__Dur(int _cSub,  int _cmt, double _amt, double __t){
 return 0.0;
}
// Model Times
void m1__mtime(int _cSub, double *_mtime){
}
// Matrix Exponential (0)
void m1__ME(int _cSub, double _t, double __t, double *_mat, const double *__zzStateVar__){
  int _itwhile = 0;
  (void)_itwhile;
  double t = __t + _solveData->subjects[_cSub].curShift;
  (void)t;
  }
// Inductive linearization Matf
void m1__IndF(int _cSub, double _t, double __t, double *_matf){
 int _itwhile = 0;
  (void)_itwhile;
  double t = __t + _solveData->subjects[_cSub].curShift;
  (void)t;
  }
extern SEXP m1__model_vars(){
  int pro=0;
  SEXP _mv = PROTECT(_rxGetModelLib("m1__model_vars"));pro++;
  if (!_rxIsCurrentC(_mv)){
    SEXP hash    = PROTECT(allocVector(STRSXP, 1));pro++;
#define __doBuf__  sprintf(buf, "un]\"BAAA@QBtHACAAAAAAAN/rAAAv7#aT)>S4.KA&eZ(ftyxTANR*5cOQQ+Wkc$w)F}QId{0XhI6IM8w5~`Mr<fpJ<RfmMz%%PJ:gwahfE6.zsmB&ws@Y^abBPBSKr\?k_ysClK_4R=d6$L.iZ@#arw%%1oAF,}TZ.oRV<4//k})HTnH.{bc;5sqo[r8]UyPd8~9JV3@3c{d!n6ZeB4/N;g8!SP]8>g:Uyh`f>|^[tfRIx@1}):/HD03[{,\"EEme\"JX3(Iw3{dg[$,VN2^cSrd!2G:;/$}(qNc4)^%%_Iu#Zbx7Y=)5LH4M50c~B,I9[lXRh:wx2%%Qc4kG3nA#gj7yKO+i@Y|B08hU6=a\"Q2d\?AoU,BUW2kU$M=AP\?A)VbkQT!o2^/H4beAWQDQZ|8$[zjG=>._lyS<kDH}{DH4}N]RK)e6]pagP[O(q&nh\"iJI$%%)[8V$Tu%%v=Y_[N;hg!U#OY<f`Wo\"CO5GG^Quc2:[t,}LlPiJ4,g^75>!vwPRa=IK,9N\"1@!#KMx]eUs+3d**,>;Ei\?\?lB2[^blR5W3k!9):\?XXKw]WJ;;%%:S{zb5gb$03V+w*K9QSP*~~AvMr|3%%Q4wC50h{7@},Nswsbd1H5*9KnZ(dgF:N8im5<be)9^gpp0#M(usLR4S3I\"Qmb~z.d.E&64z/@#4kGb[51}(AkI!p4cE4YC\"CYnQDt)7+BXtxwPDCDS!Xk;bqP`0]CfL_|bUhi)7m4brJZX7hlmSY)Z8iY5S|C|Be]Oo#n61;I04ut<D>wdR2~aIn~S}%%X&(_x307(<x]$_#A/*+v/geu;cLV8!<WWC_Fs[qAsiBp,XxtonZ#kj@9/XOqvILEujA4XNmy*a/ABxNw^D{\?/OK;(Oo2j{s%%zp|fC~Z~EY)eA");
    char buf[901];
    __doBuf__
#undef __doBuf__
    SET_STRING_ELT(hash, 0, mkChar(buf));
    SEXP lst      = PROTECT(_rxQr(hash));pro++;
    _assign_ptr(lst);
    UNPROTECT(pro);
    return lst;
  } else {
    UNPROTECT(pro);
    return _mv;
  }
}
extern void m1__dydt_lsoda(int *neq, double *t, double *A, double *DADT)
{
  m1__dydt(neq, *t, A, DADT);
}
extern int m1__dydt_liblsoda(double __t, double *y, double *ydot, void *data)
{
  int *neq = (int*)(data);
  m1__dydt(neq, __t, y, ydot);
  return(0);
}
extern void m1__calc_jac_lsoda(int *neq, double *t, double *A,int *ml, int *mu, double *JAC, int *nrowpd){
  // Update all covariate parameters
  m1__calc_jac(neq, *t, A, JAC, *nrowpd);
}

//Create function to call from R's main thread that assigns the required functions. Sometimes they don't get assigned.
extern void m1__assignFuns(){
  _assignFuns();
}

//Initialize the dll to match RxODE's calls
void R_init0_m1_(){
  // Get C callables on load; Otherwise it isn't thread safe
  _assignFuns();
  R_RegisterCCallable("m1_","m1__assignFuns", (DL_FUNC) m1__assignFuns);
  R_RegisterCCallable("m1_","m1__inis",(DL_FUNC) m1__inis);
  R_RegisterCCallable("m1_","m1__dydt",(DL_FUNC) m1__dydt);
  R_RegisterCCallable("m1_","m1__calc_lhs",(DL_FUNC) m1__calc_lhs);
  R_RegisterCCallable("m1_","m1__calc_jac",(DL_FUNC) m1__calc_jac);
  R_RegisterCCallable("m1_","m1__dydt_lsoda", (DL_FUNC) m1__dydt_lsoda);
  R_RegisterCCallable("m1_","m1__calc_jac_lsoda", (DL_FUNC) m1__calc_jac_lsoda);
  R_RegisterCCallable("m1_","m1__ode_solver_solvedata", (DL_FUNC) m1__ode_solver_solvedata);
  R_RegisterCCallable("m1_","m1__ode_solver_get_solvedata", (DL_FUNC) m1__ode_solver_get_solvedata);
  R_RegisterCCallable("m1_","m1__F", (DL_FUNC) m1__F);
  R_RegisterCCallable("m1_","m1__Lag", (DL_FUNC) m1__Lag);
  R_RegisterCCallable("m1_","m1__Rate", (DL_FUNC) m1__Rate);
  R_RegisterCCallable("m1_","m1__Dur", (DL_FUNC) m1__Dur);
  R_RegisterCCallable("m1_","m1__mtime", (DL_FUNC) m1__mtime);
  R_RegisterCCallable("m1_","m1__ME", (DL_FUNC) m1__ME);
  R_RegisterCCallable("m1_","m1__IndF", (DL_FUNC) m1__IndF);
  R_RegisterCCallable("m1_","m1__dydt_liblsoda", (DL_FUNC) m1__dydt_liblsoda);
}
//Initialize the dll to match RxODE's calls
void R_init_m1_(DllInfo *info){
  // Get C callables on load; Otherwise it isn't thread safe
  R_init0_m1_();
  static const R_CallMethodDef callMethods[]  = {
    {"m1__model_vars", (DL_FUNC) &m1__model_vars, 0},
    {NULL, NULL, 0}
  };

  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info,FALSE);
}

void R_unload_m1_ (DllInfo *info){
  // Free resources required for single subject solve.
  SEXP _mv = PROTECT(_rxGetModelLib("m1__model_vars"));
  if (!isNull(_mv)){
    _rxRmModelLib("m1__model_vars");
  }
  UNPROTECT(1);
}
