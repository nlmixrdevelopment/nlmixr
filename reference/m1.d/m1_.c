#include <RxODE_model_shared.h>
#define __MAX_PROD__ 0
#define _CMT CMT
#define _SYNC_simeps_ for (int _svari=_solveData->neps; _svari--;){  if (_solveData->svar[_svari] == 0) {KA = _PP[0];};   if (_solveData->svar[_svari] == 1) {KE = _PP[1];}; }
#define _SYNC_simeta_ for (int _ovari=_solveData->neta; _ovari--;){  if (_solveData->ovar[_ovari] == 0) {KA = _PP[0];};   if (_solveData->ovar[_ovari] == 1) {KE = _PP[1];}; }
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
#define __doBuf__  sprintf(buf, "un]\"BAAA@QBtHACAAAAAAAv(rAAAv7#aT):ST,KACeV(ftyxTANR*5cO28nu*)S](X2juOra5QEdEGeY;~T_av<gD$/0Kwk9bvBDeM#isY|5M=tQ`$<_{WnEOBaKQ@k__g>n{SvJ6PU}(6G0}$8V},H|=C[V}NA$bL=d\"5>}\?D8of_dSiRJqCbk^oq[z;!hWGym>[v]rZ68nS]FQ:qy^oz*ITx*6CdPe2.}V,T`ogu\?mKsFe5HkhMH86!vnU\"vJX]CsIcRa9`O&.8[Pq[5gpbl<jLW|!(T#u[9L3ti!d!cLl:Rg]xDoz%%Z(Fzt}76uZa9K^lnMGU%%(~lNyPJ#ZFBxvRH^F90AH2;)kdl;.~et}`)zDua3\?Jt3Rmjf7yI=T7cz@m(=Ptm3O<N492(xwy,d5^6,{:k@k\?Lh7CA]9SKF|qvBI;|TaBx{c[l:O:>*N_arc,6|agt^4r!GTip[d\"cSx8R;eNUtBozUwOJiQOoQM=*S]PP>Qul/_ju!+d[.PnqYK43E/N6T$j,r>e.h^|g1/]+,\?]im|_GW&Q/oI~6a,Nisv81MC)`7@pL&#*0AL5T,}!3VP\"fd|~]TmXN}l~WAbig.1Nv,oh;lY$!xb4nz>2]|U5ej\?&uh}uYo~4$vi|T^{c`;@A#S5{~ipqQRrl66jK1JV7SHRyv\"a@j*z*xYYMAjA+O&gGDu{fv}eKC+thZnnSySz#}NYHMc40B(>5=e{4Uvz`4[WAMqR^FB$B8rV0S9*7F&JBooYGyvn9ptp=oC5Ai&,RO}lH&Ib9~u.#Rm>IaR[&.\"3_wTP;gTI=5)Dq5xq~6b~bzs{q1_D:0#O\?Pbb+)FZ;r8*hq{xGO`qv,(mW>DtWw!9cIH[H\?\"S{UX:CHu#/9<>\?OlLw{A");
    char buf[892];
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
