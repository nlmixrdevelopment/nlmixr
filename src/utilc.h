#ifndef __UTILC_H__
#define __UTILC_H__

#if defined(__cplusplus)
extern "C" {
#endif

  void setSilentErr(int silent);
  SEXP _nlmixr_setSilentErr(SEXP in);
  void RSprintf(const char *format, ...);
  SEXP _nlmixr_powerL(SEXP xS, SEXP lambdaS, SEXP yjS, SEXP lowS, SEXP hiS);
  SEXP _nlmixr_powerD(SEXP xS, SEXP lambdaS, SEXP yjS, SEXP lowS, SEXP hiS);
  SEXP getDfSubsetVars(SEXP ipred, SEXP lhs);

#if defined(__cplusplus)
}
#endif 


#endif
