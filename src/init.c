#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void hermite_ek_compute_(void *, void *, void *);
extern void parse_ode(void *, void *, void *, void *);
extern void parse_pars(void *, void *, void *, void *);

/* Internal C calls, should not be called outside of C code. */
typedef void (*S_fp) (double *, double *);
extern void nelder_fn(S_fp func, int n, double *start, double *step,
	       int itmax, double ftol_rel, double rcoef, double ecoef, double ccoef,
	       int *iconv, int *it, int *nfcall, double *ynewlo, double *xmin,
	       int *iprint);

/* .Call calls */
extern SEXP neldermead_wrap(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
/* extern SEXP n1qn1_wrap(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); */
extern SEXP nlmixr_lin_cmt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _nlmixr_llik_binomial_c(SEXP, SEXP, SEXP);
extern SEXP slice_wrap(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP _nlmixr_llik_poisson(SEXP, SEXP);
extern SEXP _nlmixr_llik_normal(SEXP, SEXP);
extern SEXP _nlmixr_llik_betabinomial(SEXP, SEXP, SEXP);
extern SEXP _nlmixr_llik_student_t(SEXP, SEXP);
extern SEXP _nlmixr_llik_beta(SEXP, SEXP);
extern SEXP _nlmixr_lin_cmt_stan(SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP);
extern SEXP _nlmixr_llik_neg_binomial(SEXP, SEXP);
extern SEXP _nlmixr_grFOCEi(SEXP, SEXP);
extern SEXP _nlmixr_sFOCEi(SEXP, SEXP);

// FOCEi
extern SEXP _nlmixr_rxGrad(SEXP rhoSEXP);
extern SEXP _nlmixr_rxInner(SEXP etanewsSEXP, SEXP rhoSEXP);
extern SEXP _nlmixr_rxInnerNum(SEXP etanewsSEXP, SEXP rhoSEXP);
extern SEXP _nlmixr_rxHessian(SEXP rhoSEXP);
extern SEXP _nlmixr_RxODE_focei_eta_lik(SEXP sexp_etaSEXP, SEXP sexp_rhoSEXP);
extern SEXP _nlmixr_RxODE_focei_eta_lp(SEXP sexp_etaSEXP, SEXP sexp_rhoSEXP);
extern SEXP _nlmixr_RxODE_focei_eta(SEXP fstrSEXP);
extern SEXP _nlmixr_RxODE_focei_finalize_llik(SEXP rhoSEXP);
extern SEXP _nlmixr_RxODE_finalize_log_det_OMGAinv_5(SEXP rhoSEXP);
extern SEXP _nlmixr_rxDetaDomega(SEXP rhoSEXP);
extern SEXP _nlmixr_rxOuter_(SEXP rhoSEXP);
extern SEXP _nlmixr_rxDetaDtheta(SEXP rhoSEXP);
extern SEXP _nlmixr_rxOuter(SEXP rhoSEXP);
extern SEXP _nlmixr_rxUpdateEtas(SEXP DnDhSSEXP, SEXP DhSSEXP, SEXP initSSEXP, SEXP acceptNSSEXP);
extern SEXP _nlmixr_chkSolvedInf(SEXP, SEXP);
extern SEXP _nlmixr_chkSortIDTime(SEXP _id,SEXP _time);

static const R_CMethodDef CEntries[] = {
    {"hermite_ek_compute_",     (DL_FUNC) &hermite_ek_compute_,      3},
    {"parse_ode",               (DL_FUNC) &parse_ode,                4},
    {"parse_pars",              (DL_FUNC) &parse_pars,               4},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
  {"_nlmixr_chkSolvedInf",(DL_FUNC) &_nlmixr_chkSolvedInf, 2},
  {"_nlmixr_rxInner", (DL_FUNC) &_nlmixr_rxInner, 2},
  {"_nlmixr_rxInnerNum", (DL_FUNC) &_nlmixr_rxInnerNum, 2},
  {"_nlmixr_rxGrad", (DL_FUNC) &_nlmixr_rxGrad, 1},
  {"_nlmixr_rxHessian", (DL_FUNC) &_nlmixr_rxHessian, 1},
  {"_nlmixr_RxODE_focei_eta_lik", (DL_FUNC) &_nlmixr_RxODE_focei_eta_lik, 2},
  {"_nlmixr_RxODE_focei_eta_lp", (DL_FUNC) &_nlmixr_RxODE_focei_eta_lp, 2},
  {"_nlmixr_RxODE_focei_eta", (DL_FUNC) &_nlmixr_RxODE_focei_eta, 1},
  {"_nlmixr_RxODE_focei_finalize_llik", (DL_FUNC) &_nlmixr_RxODE_focei_finalize_llik, 1},
  {"_nlmixr_rxDetaDomega", (DL_FUNC) &_nlmixr_rxDetaDomega, 1},
  {"_nlmixr_rxOuter_", (DL_FUNC) &_nlmixr_rxOuter_, 1},
  {"_nlmixr_rxDetaDtheta", (DL_FUNC) &_nlmixr_rxDetaDtheta, 1},
  {"_nlmixr_rxOuter", (DL_FUNC) &_nlmixr_rxOuter, 1},
  {"_nlmixr_rxUpdateEtas", (DL_FUNC) &_nlmixr_rxUpdateEtas, 4},
  {"neldermead_wrap",      (DL_FUNC) &neldermead_wrap,      11},
  /* {"n1qn1_wrap",           (DL_FUNC) &n1qn1_wrap,           13}, */
  {"nlmixr_lin_cmt",       (DL_FUNC) &nlmixr_lin_cmt,        9},
  {"_nlmixr_lin_cmt_stan",  (DL_FUNC) &_nlmixr_lin_cmt_stan,   9},
  {"_nlmixr_llik_binomial_c", (DL_FUNC) &_nlmixr_llik_binomial_c,  3},
  {"_nlmixr_llik_poisson",  (DL_FUNC) &_nlmixr_llik_poisson,   2},
  {"_nlmixr_llik_normal",   (DL_FUNC) &_nlmixr_llik_normal,    2},
  {"_nlmixr_llik_betabinomial", (DL_FUNC) &_nlmixr_llik_betabinomial, 3},
  {"_nlmixr_llik_student_t",  (DL_FUNC) &_nlmixr_llik_student_t, 2},
  {"_nlmixr_llik_beta",     (DL_FUNC) &_nlmixr_llik_beta, 2},
  {"_nlmixr_llik_neg_binomial", (DL_FUNC) &_nlmixr_llik_neg_binomial, 2},
  {"_nlmixr_grFOCEi",(DL_FUNC) &_nlmixr_grFOCEi, 2},
  {"_nlmixr_sFOCEi",(DL_FUNC) &_nlmixr_sFOCEi, 2},
  {"_nlmixr_chkSortIDTime",(DL_FUNC) &_nlmixr_chkSortIDTime, 2},
  {"slice_wrap",           (DL_FUNC) &slice_wrap,            7},
  {NULL, NULL, 0}
};

void R_init_nlmixr(DllInfo *dll)
{
  R_RegisterCCallable("nlmixr","nelder_fn", (DL_FUNC) &nelder_fn);
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, TRUE);
  R_forceSymbols(dll,FALSE);
}
