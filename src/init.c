#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void hermite_ek_compute_(void *, void *, void *);
extern void parse_ode(void *, void *, void *, void *);
extern void parse_pars(void *, void *, void *, void *);

/* .Call calls */
extern SEXP neldermead_wrap(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP nlmixr_lin_cmt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP nlmixr_llik_binomial_c(SEXP, SEXP, SEXP);
extern SEXP slice_wrap(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP nlmixr_llik_poisson(SEXP, SEXP);
extern SEXP nlmixr_llik_normal(SEXP, SEXP);
extern SEXP nlmixr_llik_betabinomial(SEXP, SEXP, SEXP);
extern SEXP nlmixr_llik_student_t(SEXP, SEXP);
extern SEXP nlmixr_llik_beta(SEXP, SEXP);
extern SEXP nlmixr_lin_cmt_stan(SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP);
extern SEXP nlmixr_llik_neg_binomial(SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"hermite_ek_compute_",     (DL_FUNC) &hermite_ek_compute_,      3},
    {"parse_ode",               (DL_FUNC) &parse_ode,                4},
    {"parse_pars",              (DL_FUNC) &parse_pars,               4},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"neldermead_wrap",      (DL_FUNC) &neldermead_wrap,      11},
    {"nlmixr_lin_cmt",       (DL_FUNC) &nlmixr_lin_cmt,        9},
    {"nlmixr_lin_cmt_stan",  (DL_FUNC) &nlmixr_lin_cmt_stan,   9},
    {"nlmixr_llik_binomial_c", (DL_FUNC) &nlmixr_llik_binomial_c,  3},
    {"nlmixr_llik_poisson",  (DL_FUNC) &nlmixr_llik_poisson,   2},
    {"nlmixr_llik_normal",   (DL_FUNC) &nlmixr_llik_normal,    2},
    {"nlmixr_llik_betabinomial", (DL_FUNC) &nlmixr_llik_betabinomial, 3},
    {"nlmixr_llik_student_t",  (DL_FUNC) &nlmixr_llik_student_t, 2},
    {"nlmixr_llik_beta",     (DL_FUNC) &nlmixr_llik_beta, 2},
    {"nlmixr_llik_neg_binomial", (DL_FUNC) &nlmixr_llik_neg_binomial, 2},
    {"slice_wrap",           (DL_FUNC) &slice_wrap,            7},
    {NULL, NULL, 0}
};

void R_init_nlmixr(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll,TRUE);
}
