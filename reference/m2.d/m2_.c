#define _getRxSolve_ _getRxSolve_1614204903
#define simeps simeps1614204903
#define simeta simeta1614204903
#define _solveData _solveData1614204903
#define _assign_ptr _assign_ptr1614204903
#define _rxRmModelLib _rxRmModelLib1614204903
#define _rxGetModelLib _rxGetModelLib1614204903
#define _old_c _old_c1614204903
#define _ptrid _ptrid1614204903
#define _rxIsCurrentC _rxIsCurrentC1614204903
#define _sumPS _sumPS1614204903
#define _prodPS _prodPS1614204903
#define _prodType _prodType1614204903
#define _sumType _sumType1614204903
#define _update_par_ptr _update_par_ptr1614204903
#define _getParCov _getParCov1614204903
#define linCmtA linCmtA1614204903
#define linCmtC linCmtC1614204903
#define linCmtB linCmtB1614204903
#define _RxODE_rxAssignPtr _RxODE_rxAssignPtr1614204903
#define _rxQr _rxQr1614204903
#define phi phi1614204903
#define logit logit1614204903
#define expit expit1614204903
#define gammap gammap1614204903
#define gammaq gammaq1614204903
#define lowergamma lowergamma1614204903
#define uppergamma uppergamma1614204903
#define gammapInv gammapInv1614204903
#define gammapDer gammapDer1614204903
#define gammapInva gammapInva1614204903
#define gammaqInv gammaqInv1614204903
#define gammaqInva gammaqInva1614204903
#define rxnorm rxnorm1614204903
#define rxnormV rxnormV1614204903
#define rxbinom rxbinom1614204903
#define rxcauchy rxcauchy1614204903
#define rxchisq rxchisq1614204903
#define rxexp rxexp1614204903
#define rxf rxf1614204903
#define rxgeom rxgeom1614204903
#define rxgamma rxgamma1614204903
#define rxbeta rxbeta1614204903
#define rxpois rxpois1614204903
#define rxt_ rxt_1614204903
#define rxunif rxunif1614204903
#define rxweibull rxweibull1614204903
#define rinorm rinorm1614204903
#define rinormV rinormV1614204903
#define ribinom ribinom1614204903
#define ricauchy ricauchy1614204903
#define richisq richisq1614204903
#define riexp riexp1614204903
#define rif rif1614204903
#define rigeom rigeom1614204903
#define rigamma rigamma1614204903
#define ribeta ribeta1614204903
#define ripois ripois1614204903
#define rit_ rit_1614204903
#define riunif riunif1614204903
#define riweibull riweibull1614204903
#define _compareFactorVal _compareFactorVal1614204903
#define _sum _sum1614204903
#define _sign _sign1614204903
#define _prod _prod1614204903
#define _max _max1614204903
#define _min _min1614204903
#define _transit4P _transit4P1614204903
#define _transit3P _transit3P1614204903
#define _assignFuns0 _assignFuns01614204903
#define _assignFuns _assignFuns1614204903
#define _getRxSolve_ _getRxSolve_1614204903
#include <RxODE_model_shared.h>
#define __MAX_PROD__ 0
#define _CMT CMT
#define _SYNC_simeps_ for (int _svari=_solveData->neps; _svari--;){  if (_solveData->svar[_svari] == 0) {V = _PP[0];};   if (_solveData->svar[_svari] == 1) {KA = _PP[1];};   if (_solveData->svar[_svari] == 2) {KE = _PP[2];}; }
#define _SYNC_simeta_ for (int _ovari=_solveData->neta; _ovari--;){  if (_solveData->ovar[_ovari] == 0) {V = _PP[0];};   if (_solveData->ovar[_ovari] == 1) {KA = _PP[1];};   if (_solveData->ovar[_ovari] == 2) {KE = _PP[2];}; }
#include "extraC.h"
_getRxSolve_t _getRxSolve_;
_simfun simeps;
_simfun simeta;
rx_solve *_solveData=NULL;
RxODE_assign_ptr _assign_ptr=NULL;
_rxRmModelLibType _rxRmModelLib=NULL;
_rxGetModelLibType _rxGetModelLib=NULL;
RxODE_ode_solver_old_c _old_c=NULL;
RxODE_fn0i _ptrid=NULL;
_rxIsCurrentC_type _rxIsCurrentC=NULL;
_rxSumType _sumPS=NULL;
_rxProdType _prodPS=NULL;
RxODE_fn0i _prodType=NULL;
RxODE_fn0i _sumType=NULL;
_update_par_ptr_p _update_par_ptr=NULL;
_getParCov_p _getParCov=NULL;
linCmtA_p linCmtA;
linCmtA_p linCmtC;
linCmtB_p linCmtB;
_rx_asgn _RxODE_rxAssignPtr=NULL;
_rx_asgn _rxQr=NULL;
RxODE_fn phi;
RxODE_fn3 logit;
RxODE_fn3 expit;
RxODE_fn2 gammap;
RxODE_fn2 gammaq;
RxODE_fn2 lowergamma;
RxODE_fn2 uppergamma;
RxODE_fn2 gammapInv;
RxODE_fn2 gammapDer;
RxODE_fn2 gammapInva;
RxODE_fn2 gammaqInv;
RxODE_fn2 gammaqInva;
RxODEi_fn2 rxnorm;
RxODEi_fn2 rxnormV;
RxODEi_rxbinom rxbinom;
RxODEi_fn2 rxcauchy;
RxODEi_fn rxchisq;
RxODEi_fn rxexp;
RxODEi_fn2 rxf;
RxODEi_ifn rxgeom;
RxODEi_fn2 rxgamma;
RxODEi_fn2 rxbeta;
RxODEi_ifn rxpois;
RxODEi_fn rxt_;
RxODEi_fn2 rxunif;
RxODEi_fn2 rxweibull;
RxODEi2_fn2 rinorm;
RxODEi2_fn2 rinormV;
RxODEi2_ribinom ribinom;
RxODEi2_fn2 ricauchy;
RxODEi2_fn richisq;
RxODEi2_fn riexp;
RxODEi2_fn2 rif;
RxODEi2_ifn rigeom;
RxODEi2_fn2 rigamma;
RxODEi2_fn2 ribeta;
RxODEi2_ifn ripois;
RxODEi2_fn rit_;
RxODEi2_fn2 riunif;
RxODEi2_fn2 riweibull;
RxODE_compareFactorVal_fn _compareFactorVal;
double _prod(double *input, double *p, int type, int n, ...){
  va_list valist;
  va_start(valist, n);
  for (unsigned int i = 0; i < n; i++){
    input[i] = va_arg(valist, double);
  }
  va_end(valist);
  return _prodPS(input, p, n, type);
}
double _sum(double *input, double *pld, int m, int type, int n, ...){
  va_list valist;
  va_start(valist, n);
  for (unsigned int i = 0; i < n; i++){
    input[i] = va_arg(valist, double);
  }
  va_end(valist);
  double ret = _sumPS(input, n, pld, m, type);
  if (type == 2 && m < 0){
    for (int i = -m; i--;){
      pld[i] = 0.0;
    }
  }
  return ret;
}
double _sign(unsigned int n, ...) {
  va_list valist;
  va_start(valist, n);
  double s = 1;
  for (unsigned int i = 0; i < n; i++) {
    s = sign(va_arg(valist, double))*s;
    if (s == 0){
      break;
    }
  }
  va_end(valist);
  return s;
}
double _max(unsigned int n, ...) {
  va_list valist;
  va_start(valist, n);
  double mx = NA_REAL;
  double tmp = 0;
  if (n >= 1){
    mx = va_arg(valist, double);
    for (unsigned int i = 1; i < n; i++) {
      tmp = va_arg(valist, double);
      if (tmp>mx) mx=tmp;
    }
    va_end(valist);
  }
  return mx;
}
double _min(unsigned int n, ...){
  va_list valist;
  va_start(valist, n);
  double mn = NA_REAL;
  double tmp = 0;
  if (n >= 1){
    mn = va_arg(valist, double);
    for (unsigned int i = 1; i < n; i++){
      tmp = va_arg(valist, double);
      if (tmp<mn) mn=tmp;
    }
    va_end(valist);
  }
  return mn;
}
double _transit4P(double t, unsigned int id, double n, double mtt, double bio){
  double ktr = (n+1)/mtt;
  double lktr = log(n+1)-log(mtt);
  double tc = (t-(_solveData->subjects[id].tlast));
  return exp(log(bio*(_solveData->subjects[id].podo))+lktr+n*(lktr+log(tc))-ktr*(tc)-lgamma1p(n));
}
double _transit3P(double t, unsigned int id, double n, double mtt){
  double ktr = (n+1)/mtt;
  double lktr = log(n+1)-log(mtt);
  double tc = (t-(_solveData->subjects[id].tlast));
  return exp(log(_solveData->subjects[id].podo)+lktr+n*(lktr+log(tc))-ktr*(tc)-lgamma1p(n));
}
void _assignFuns0() {
  _getRxSolve_ = (_getRxSolve_t) R_GetCCallable("RxODE","getRxSolve_");
  _assign_ptr=(RxODE_assign_ptr) R_GetCCallable("RxODE","RxODE_assign_fn_pointers");
  _rxRmModelLib=(_rxRmModelLibType) R_GetCCallable("RxODE","rxRmModelLib");
  _rxGetModelLib=(_rxGetModelLibType) R_GetCCallable("RxODE","rxGetModelLib");
  _RxODE_rxAssignPtr=(_rx_asgn)R_GetCCallable("RxODE","_RxODE_rxAssignPtr");
  _rxQr=(_rx_asgn)R_GetCCallable("RxODE","_RxODE_rxQr");
  _rxIsCurrentC = (_rxIsCurrentC_type)R_GetCCallable("RxODE","rxIsCurrentC");
  _sumPS  = (_rxSumType) R_GetCCallable("PreciseSums","PreciseSums_sum_r");
  _prodPS = (_rxProdType) R_GetCCallable("PreciseSums","PreciseSums_prod_r");
  _prodType=(RxODE_fn0i)R_GetCCallable("PreciseSums", "PreciseSums_prod_get");
  _sumType=(RxODE_fn0i)R_GetCCallable("PreciseSums", "PreciseSums_sum_get");
  _ptrid=(RxODE_fn0i)R_GetCCallable("RxODE", "RxODE_current_fn_pointer_id");
  linCmtA=(linCmtA_p)R_GetCCallable("RxODE", "linCmtA");
  linCmtB=(linCmtB_p)R_GetCCallable("RxODE", "linCmtB");
  linCmtC=(linCmtA_p)R_GetCCallable("RxODE", "linCmtC");
    
  rxnorm = (RxODEi_fn2)R_GetCCallable("RxODE", "rxnorm");
  rxnormV = (RxODEi_fn2)R_GetCCallable("RxODE", "rxnormV");
  rxbinom = (RxODEi_rxbinom)R_GetCCallable("RxODE","rxbinom") ;
  rxcauchy = (RxODEi_fn2)R_GetCCallable("RxODE","rxcauchy") ;
  rxchisq = (RxODEi_fn)R_GetCCallable("RxODE","rxchisq") ;
  rxexp = (RxODEi_fn)R_GetCCallable("RxODE","rxexp");
  rxf = (RxODEi_fn2)R_GetCCallable("RxODE","rxf") ;
  rxgeom = (RxODEi_ifn)R_GetCCallable("RxODE","rxgeom") ;
  rxgamma = (RxODEi_fn2)R_GetCCallable("RxODE","rxgamma") ;
  rxbeta = (RxODEi_fn2)R_GetCCallable("RxODE","rxbeta") ;
  rxpois = (RxODEi_ifn)R_GetCCallable("RxODE","rxpois") ;
  rxt_ = (RxODEi_fn)R_GetCCallable("RxODE","rxt_") ;
  rxunif = (RxODEi_fn2)R_GetCCallable("RxODE","rxunif") ;
  rxweibull = (RxODEi_fn2)R_GetCCallable("RxODE","rxweibull");
  rinorm = (RxODEi2_fn2)R_GetCCallable("RxODE", "rinorm");
  rinormV = (RxODEi2_fn2)R_GetCCallable("RxODE", "rinormV");
  ribinom = (RxODEi2_ribinom)R_GetCCallable("RxODE","ribinom") ;
  ricauchy = (RxODEi2_fn2)R_GetCCallable("RxODE","ricauchy") ;
  richisq = (RxODEi2_fn)R_GetCCallable("RxODE","richisq") ;
  riexp = (RxODEi2_fn)R_GetCCallable("RxODE","riexp");
  rif = (RxODEi2_fn2)R_GetCCallable("RxODE","rif") ;
  rigeom = (RxODEi2_ifn)R_GetCCallable("RxODE","rigeom") ;
  rigamma = (RxODEi2_fn2)R_GetCCallable("RxODE","rigamma") ;
  ribeta = (RxODEi2_fn2)R_GetCCallable("RxODE","ribeta") ;
  ripois = (RxODEi2_ifn)R_GetCCallable("RxODE","ripois") ;
  rit_ = (RxODEi2_fn)R_GetCCallable("RxODE","rit_") ;
  riunif = (RxODEi2_fn2)R_GetCCallable("RxODE","riunif") ;
  riweibull = (RxODEi2_fn2)R_GetCCallable("RxODE","riweibull");
    
  phi = (RxODE_fn)R_GetCCallable("RxODE","phi");
  gammap = (RxODE_fn2) R_GetCCallable("RxODE","gammap");
  gammaq = (RxODE_fn2) R_GetCCallable("RxODE","gammaq");
  gammapInv = (RxODE_fn2) R_GetCCallable("RxODE","gammapInv");
  gammapInva = (RxODE_fn2) R_GetCCallable("RxODE","gammapInva");
  gammaqInv = (RxODE_fn2) R_GetCCallable("RxODE","gammaqInv");
  gammaqInva = (RxODE_fn2) R_GetCCallable("RxODE","gammaqInva");
  uppergamma = (RxODE_fn2) R_GetCCallable("RxODE","uppergamma");
  lowergamma = (RxODE_fn2) R_GetCCallable("RxODE","lowergamma");
  gammapDer  = (RxODE_fn2) R_GetCCallable("RxODE","gammapDer");
  logit = (RxODE_fn3) R_GetCCallable("RxODE", "logit");
  expit = (RxODE_fn3) R_GetCCallable("RxODE", "expit");
  simeta =(_simfun) R_GetCCallable("RxODE", "simeta");
  simeps =(_simfun) R_GetCCallable("RxODE", "simeps");
  _compareFactorVal=(RxODE_compareFactorVal_fn) R_GetCCallable("RxODE", "compareFactorVal");
  _update_par_ptr = (_update_par_ptr_p) R_GetCCallable("RxODE","_update_par_ptr");
  _getParCov = (_getParCov_p) R_GetCCallable("RxODE","_getParCov");
  _solveData = _getRxSolve_();
}
void _assignFuns() {
  if (_assign_ptr == NULL){
    _assignFuns0();
  }
}
extern void  m2__ode_solver_solvedata (rx_solve *solve){
  _solveData = solve;
}
extern rx_solve *m2__ode_solver_get_solvedata(){
  return _solveData;
}
SEXP m2__model_vars();


// prj-specific differential eqns
void m2__dydt(int *_neq, double __t, double *__zzStateVar__, double *__DDtStateVar__)
{
  int _itwhile = 0;
  (void)_itwhile;
  int _cSub = _neq[1];
  double t = __t + _solveData->subjects[_neq[1]].curShift;
  (void)t;
    double C2;
  double centr;
  double V;
  double depot;
  double KA;
  double KE;

  (void)t;
  (void)C2;
  (void)centr;
  (void)V;
  (void)depot;
  (void)KA;
  (void)KE;

  C2 = _PL[0];

  _update_par_ptr(__t, _cSub, _solveData, _idx);
  V = _PP[0];
  KA = _PP[1];
  KE = _PP[2];

  depot = __zzStateVar__[0]*((double)(_ON[0]));
  centr = __zzStateVar__[1]*((double)(_ON[1]));

  C2=centr/safe_zero(V);
  __DDtStateVar__[0] = ((double)(_ON[0]))*(_IR[0] -KA*depot);
  __DDtStateVar__[1] = ((double)(_ON[1]))*(_IR[1] + KA*depot-KE*centr);
  (&_solveData->subjects[_cSub])->dadt_counter[0]++;
}

// Jacobian derived vars
void m2__calc_jac(int *_neq, double __t, double *__zzStateVar__, double *__PDStateVar__, unsigned int __NROWPD__) {
  int _itwhile = 0;
  (void)_itwhile;
    int _cSub=_neq[1];
  double t = __t + _solveData->subjects[_neq[1]].curShift;
  (void)t;
    (&_solveData->subjects[_cSub])->jac_counter[0]++;
}
// Functional based initial conditions.
void m2__inis(int _cSub, double *__zzStateVar__){
  int _itwhile = 0;
  (void)_itwhile;
  
}
// prj-specific derived vars
void m2__calc_lhs(int _cSub, double __t, double *__zzStateVar__, double *_lhs) {
    int _itwhile = 0;
  (void)_itwhile;
  double t = __t + _solveData->subjects[_cSub].curShift;
  (void)t;
    double  __DDtStateVar_0__;
  double  __DDtStateVar_1__;
  double C2;
  double centr;
  double V;
  double depot;
  double KA;
  double KE;

  (void)t;
  (void)__DDtStateVar_0__;
  (void)__DDtStateVar_1__;
  (void)C2;
  (void)centr;
  (void)V;
  (void)depot;
  (void)KA;
  (void)KE;

  C2 = _PL[0];

  _update_par_ptr(__t, _cSub, _solveData, _idx);
  V = _PP[0];
  KA = _PP[1];
  KE = _PP[2];

  depot = __zzStateVar__[0]*((double)(_ON[0]));
  centr = __zzStateVar__[1]*((double)(_ON[1]));

  C2=centr/safe_zero(V);
  __DDtStateVar_0__ = ((double)(_ON[0]))*(_IR[0] -KA*depot);
  __DDtStateVar_1__ = ((double)(_ON[1]))*(_IR[1] + KA*depot-KE*centr);

  _lhs[0]=C2;
}
// Functional based bioavailability
double m2__F(int _cSub,  int _cmt, double _amt, double __t, double *__zzStateVar__){
 return _amt;
}
// Functional based absorption lag
double m2__Lag(int _cSub,  int _cmt, double __t, double *__zzStateVar__){
 return __t;
}
// Modeled zero-order rate
double m2__Rate(int _cSub,  int _cmt, double _amt, double __t, double *__zzStateVar__){
 return 0.0;
}
// Modeled zero-order duration
double m2__Dur(int _cSub,  int _cmt, double _amt, double __t){
 return 0.0;
}
// Model Times
void m2__mtime(int _cSub, double *_mtime){
}
// Matrix Exponential (0)
void m2__ME(int _cSub, double _t, double __t, double *_mat, const double *__zzStateVar__){
  int _itwhile = 0;
  (void)_itwhile;
  double t = __t + _solveData->subjects[_cSub].curShift;
  (void)t;
  }
// Inductive linearization Matf
void m2__IndF(int _cSub, double _t, double __t, double *_matf){
 int _itwhile = 0;
  (void)_itwhile;
  double t = __t + _solveData->subjects[_cSub].curShift;
  (void)t;
  }
extern SEXP m2__model_vars(){
  int pro=0;
  SEXP _mv = PROTECT(_rxGetModelLib("m2__model_vars"));pro++;
  if (!_rxIsCurrentC(_mv)){
    SEXP hash    = PROTECT(allocVector(STRSXP, 1));pro++;
#define __doBuf__  sprintf(buf, "un]\"BAAA@QBtHACAAAAAAAk_tAAAv7#aT)^TmeLANt%%X<W>FN\"ZS|:e\"*B8Ms]kX9s#W@wBb3~)5[~1u%%6]{_$;%%M#:EBm[F`*\?jLmU9U0GmZ2ff(g84UgeUT|YLR\"\?OjC`;H9Ao+vWVM!.vOP4i`*rp9;1y\"L$=:F_F0nGsufTxP/4^wFDN7&{z<k>N47Ql;Ud+rD8mZVMq1H1K}(_u%%yJXTjKv2\?r:x7@UD]Ri336U|i&2ELw0=0$9:4OB6*YKwI7t3\?NJ^Ms_~szohVY<5tTXp4F=\"P{/8[8Tvf_Pp0\?v1}\?,PMJ0v(towqN)N8j,Dl``*M>Moz4W2\"Q^qu&b!38a%%b{_xtT9xUEfo@1XWD01nzgH4tjNjfy6Pu;IvctqF%%;=BKycvJ$GxC4($)6clQbzdXfSg~GlA^<aH(XzH&!]npa$t=D/VdKNoxevTCXnb$wg`^~Y!Y_$htd#O>IVGg{k8UC9$Oh$F!BwH3Ro>*\"cQ%%Qu45hc5:reC^G6CIg@X/2U`B>Moz$Tff``+uQ.^]6!%%h>0V*wFCYmH@gY3Khy*ka/HRqw;]!Mg<E7($Ve/FA{Qfd00T)0KeRwPPlb,`WnlePwl5T@`kq}8!ae)in_(1J9s)f\?G(AJWTpaFgH:N|c$p`7umU[C4IJr)}[kg7$+fc+|g/\?T;t:lcR^ZuEPx)51ykJ\?\"E]^op:kd^fkj6hWxDUzlRA*gEAkVDKDmE+Axdp%%MG\"]DG{9UEvBMHJARsXd#RhFCguu/G]~#&jrF@vzsD^W\"0&IU7|\?EeX>_7u#|X~cAKKG1aZX2K6C>i.zB<\?_O,VgP]pBsvFT}kSw0]Ye2`:~j[`d/u.#>=szpHD1O@^br{4hv[f4%%PzP`Wh$xN%%Y3dI*uu])TR3[+|=g\"s+,lACN:Eoe2S4trI,^|1F/J:$3M\"AlK\"E<D");
    char buf[943];
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
extern void m2__dydt_lsoda(int *neq, double *t, double *A, double *DADT)
{
  m2__dydt(neq, *t, A, DADT);
}
extern int m2__dydt_liblsoda(double __t, double *y, double *ydot, void *data)
{
  int *neq = (int*)(data);
  m2__dydt(neq, __t, y, ydot);
  return(0);
}
extern void m2__calc_jac_lsoda(int *neq, double *t, double *A,int *ml, int *mu, double *JAC, int *nrowpd){
  // Update all covariate parameters
  m2__calc_jac(neq, *t, A, JAC, *nrowpd);
}

//Create function to call from R's main thread that assigns the required functions. Sometimes they don't get assigned.
extern void m2__assignFuns(){
  _assignFuns();
}

//Initialize the dll to match RxODE's calls
void R_init0_m2_(){
  // Get C callables on load; Otherwise it isn't thread safe
  R_RegisterCCallable("m2_","m2__assignFuns", (DL_FUNC) m2__assignFuns);
  R_RegisterCCallable("m2_","m2__inis",(DL_FUNC) m2__inis);
  R_RegisterCCallable("m2_","m2__dydt",(DL_FUNC) m2__dydt);
  R_RegisterCCallable("m2_","m2__calc_lhs",(DL_FUNC) m2__calc_lhs);
  R_RegisterCCallable("m2_","m2__calc_jac",(DL_FUNC) m2__calc_jac);
  R_RegisterCCallable("m2_","m2__dydt_lsoda", (DL_FUNC) m2__dydt_lsoda);
  R_RegisterCCallable("m2_","m2__calc_jac_lsoda", (DL_FUNC) m2__calc_jac_lsoda);
  R_RegisterCCallable("m2_","m2__ode_solver_solvedata", (DL_FUNC) m2__ode_solver_solvedata);
  R_RegisterCCallable("m2_","m2__ode_solver_get_solvedata", (DL_FUNC) m2__ode_solver_get_solvedata);
  R_RegisterCCallable("m2_","m2__F", (DL_FUNC) m2__F);
  R_RegisterCCallable("m2_","m2__Lag", (DL_FUNC) m2__Lag);
  R_RegisterCCallable("m2_","m2__Rate", (DL_FUNC) m2__Rate);
  R_RegisterCCallable("m2_","m2__Dur", (DL_FUNC) m2__Dur);
  R_RegisterCCallable("m2_","m2__mtime", (DL_FUNC) m2__mtime);
  R_RegisterCCallable("m2_","m2__ME", (DL_FUNC) m2__ME);
  R_RegisterCCallable("m2_","m2__IndF", (DL_FUNC) m2__IndF);
  R_RegisterCCallable("m2_","m2__dydt_liblsoda", (DL_FUNC) m2__dydt_liblsoda);
}
//Initialize the dll to match RxODE's calls
void R_init_m2_(DllInfo *info){
  // Get C callables on load; Otherwise it isn't thread safe
  R_init0_m2_();
  static const R_CallMethodDef callMethods[]  = {
    {"m2__model_vars", (DL_FUNC) &m2__model_vars, 0},
    {NULL, NULL, 0}
  };

  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info,FALSE);
  _assignFuns0();

}

void R_unload_m2_ (DllInfo *info){
  // Free resources required for single subject solve.
  SEXP _mv = PROTECT(_rxGetModelLib("m2__model_vars"));
  if (!isNull(_mv)){
    _rxRmModelLib("m2__model_vars");
  }
  UNPROTECT(1);
}
