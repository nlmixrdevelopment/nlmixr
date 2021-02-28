#define _getRxSolve_ _rx4333335fb0170d3ffb11e0dbcff737c40
#define simeps _rx4333335fb0170d3ffb11e0dbcff737c41
#define simeta _rx4333335fb0170d3ffb11e0dbcff737c42
#define _solveData _rx4333335fb0170d3ffb11e0dbcff737c43
#define _assign_ptr _rx4333335fb0170d3ffb11e0dbcff737c44
#define _rxRmModelLib _rx4333335fb0170d3ffb11e0dbcff737c45
#define _rxGetModelLib _rx4333335fb0170d3ffb11e0dbcff737c46
#define _old_c _rx4333335fb0170d3ffb11e0dbcff737c47
#define _ptrid _rx4333335fb0170d3ffb11e0dbcff737c48
#define _rxIsCurrentC _rx4333335fb0170d3ffb11e0dbcff737c49
#define _sumPS _rx4333335fb0170d3ffb11e0dbcff737c410
#define _prodPS _rx4333335fb0170d3ffb11e0dbcff737c411
#define _prodType _rx4333335fb0170d3ffb11e0dbcff737c412
#define _sumType _rx4333335fb0170d3ffb11e0dbcff737c413
#define _update_par_ptr _rx4333335fb0170d3ffb11e0dbcff737c414
#define _getParCov _rx4333335fb0170d3ffb11e0dbcff737c415
#define linCmtA _rx4333335fb0170d3ffb11e0dbcff737c416
#define linCmtC _rx4333335fb0170d3ffb11e0dbcff737c417
#define linCmtB _rx4333335fb0170d3ffb11e0dbcff737c418
#define _RxODE_rxAssignPtr _rx4333335fb0170d3ffb11e0dbcff737c419
#define _rxQr _rx4333335fb0170d3ffb11e0dbcff737c420
#define phi _rx4333335fb0170d3ffb11e0dbcff737c421
#define logit _rx4333335fb0170d3ffb11e0dbcff737c422
#define expit _rx4333335fb0170d3ffb11e0dbcff737c423
#define gammap _rx4333335fb0170d3ffb11e0dbcff737c424
#define gammaq _rx4333335fb0170d3ffb11e0dbcff737c425
#define lowergamma _rx4333335fb0170d3ffb11e0dbcff737c426
#define uppergamma _rx4333335fb0170d3ffb11e0dbcff737c427
#define gammapInv _rx4333335fb0170d3ffb11e0dbcff737c428
#define gammapDer _rx4333335fb0170d3ffb11e0dbcff737c429
#define gammapInva _rx4333335fb0170d3ffb11e0dbcff737c430
#define gammaqInv _rx4333335fb0170d3ffb11e0dbcff737c431
#define gammaqInva _rx4333335fb0170d3ffb11e0dbcff737c432
#define rxnorm _rx4333335fb0170d3ffb11e0dbcff737c433
#define rxnormV _rx4333335fb0170d3ffb11e0dbcff737c434
#define rxbinom _rx4333335fb0170d3ffb11e0dbcff737c435
#define rxcauchy _rx4333335fb0170d3ffb11e0dbcff737c436
#define rxchisq _rx4333335fb0170d3ffb11e0dbcff737c437
#define rxexp _rx4333335fb0170d3ffb11e0dbcff737c438
#define rxf _rx4333335fb0170d3ffb11e0dbcff737c439
#define rxgeom _rx4333335fb0170d3ffb11e0dbcff737c440
#define rxgamma _rx4333335fb0170d3ffb11e0dbcff737c441
#define rxbeta _rx4333335fb0170d3ffb11e0dbcff737c442
#define rxpois _rx4333335fb0170d3ffb11e0dbcff737c443
#define rxt_ _rx4333335fb0170d3ffb11e0dbcff737c444
#define rxunif _rx4333335fb0170d3ffb11e0dbcff737c445
#define rxweibull _rx4333335fb0170d3ffb11e0dbcff737c446
#define rinorm _rx4333335fb0170d3ffb11e0dbcff737c447
#define rinormV _rx4333335fb0170d3ffb11e0dbcff737c448
#define ribinom _rx4333335fb0170d3ffb11e0dbcff737c449
#define ricauchy _rx4333335fb0170d3ffb11e0dbcff737c450
#define richisq _rx4333335fb0170d3ffb11e0dbcff737c451
#define riexp _rx4333335fb0170d3ffb11e0dbcff737c452
#define rif _rx4333335fb0170d3ffb11e0dbcff737c453
#define rigeom _rx4333335fb0170d3ffb11e0dbcff737c454
#define rigamma _rx4333335fb0170d3ffb11e0dbcff737c455
#define ribeta _rx4333335fb0170d3ffb11e0dbcff737c456
#define ripois _rx4333335fb0170d3ffb11e0dbcff737c457
#define rit_ _rx4333335fb0170d3ffb11e0dbcff737c458
#define riunif _rx4333335fb0170d3ffb11e0dbcff737c459
#define riweibull _rx4333335fb0170d3ffb11e0dbcff737c460
#define _compareFactorVal _rx4333335fb0170d3ffb11e0dbcff737c461
#define _sum _rx4333335fb0170d3ffb11e0dbcff737c462
#define _sign _rx4333335fb0170d3ffb11e0dbcff737c463
#define _prod _rx4333335fb0170d3ffb11e0dbcff737c464
#define _max _rx4333335fb0170d3ffb11e0dbcff737c465
#define _min _rx4333335fb0170d3ffb11e0dbcff737c466
#define _transit4P _rx4333335fb0170d3ffb11e0dbcff737c467
#define _transit3P _rx4333335fb0170d3ffb11e0dbcff737c468
#define _assignFuns0 _rx4333335fb0170d3ffb11e0dbcff737c469
#define _assignFuns _rx4333335fb0170d3ffb11e0dbcff737c470
#include <RxODE_model_shared.h>
#define __MAX_PROD__ 0
#define _CMT CMT
#define _SYNC_simeps_ for (int _svari=_solveData->neps; _svari--;){  if (_solveData->svar[_svari] == 0) {KA = _PP[0];};   if (_solveData->svar[_svari] == 1) {KE = _PP[1];}; }
#define _SYNC_simeta_ for (int _ovari=_solveData->neta; _ovari--;){  if (_solveData->ovar[_ovari] == 0) {KA = _PP[0];};   if (_solveData->ovar[_ovari] == 1) {KE = _PP[1];}; }
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
#define __doBuf__  sprintf(buf, "un]\"BAAA@QBtHACAAAAAAAXLtAAAv7#aT)xTMMLAUoJX<W{FnA(PMe;RB2qNh9+dyciJoZ]#jD,Ry/7=M]2r|YlRhO>+#XatUXNY69me<Ed*)M]gS4IT*8Pn+>F\"Q\"R>o0RB}lMd\"$FY)GL!CZ]D4iL(rp9;EFAM$=W\"dwiSe&)PSX2[:{f_$w~=zY]OVIKlQl;U^)rD9^yUMqf{{EY(2YxykWTjav}/4,<6@UvO8RiQv0UGi.^~U,A+od|ru<r\?)EbLCEV\?#0\?HS>SJ&v;ILH7vG0k}y5A}n&WrRPX6)wvQC3&q#5+!JVb[_0,0H<qx:YO>\"30c/[/V8Fl31n2fK5>w_G&K<h006N1EFKE)+)VD21nzgHc[GO6OiDzPT1=C,GB1=5Ct/h{g3N=HrZ>1K2>sdeX}>c[9&xo\?+{a5ZMTWJyY##3+:Y<=I=bX5kr4AUj{BI\"RL<~21st:!5\?Xpq)pMoMn7*zU})UV8U}h=MP;<k9ED#b%%7z\"7;8p^t[vaY:FWGdM#v>VCxZaA+y3l*$ry3iS,J\",~guQ#O{!#&3I{&9(1mz3d\"z(fb;cOw/Td*IR66W08ptdceKsu>_.$RBgqh0H\"s/+%%<b|4/@2@Fv^:b8nB4p%%Ao1#xGDm6TL}).bHrZQ\?J(z!gW9^`0Kh)Vh01nMKI[}jG<b\?E]aQ5SlSivPmSnNc9xn,@LXd~1`M#XgAv(IA$w#;*Aehmc\?Dw4&yHzUc>n9LMu@WwDq`=I5yeIjhxMSm`O&FVaixZ0A%%Qhd^T`,zcxQP`\?C@@whV=k1|(,7[7d\"t#\"gdO:bpMFsP3Z#+(rVOiT=U=W1n\?8z]|gLb[O]Y6\?)bH&9C@e/h9e8cE4}o%%{DyWKa`acoi#T8y6{QQ,$U|fBi@[%%\?PkkTB^h2|/!qW^5I31AAdj\?X/D");
    char buf[922];
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
  _assignFuns0();

}

void R_unload_m1_ (DllInfo *info){
  // Free resources required for single subject solve.
  SEXP _mv = PROTECT(_rxGetModelLib("m1__model_vars"));
  if (!isNull(_mv)){
    _rxRmModelLib("m1__model_vars");
  }
  UNPROTECT(1);
}
