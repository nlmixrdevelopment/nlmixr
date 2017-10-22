// Adopted by Wenping Wang from Rcpp.
// written originally by Dirk Eddelbuettel <edd@debian.org>


#ifndef Rcpp_DE_evaluate_h_
#define Rcpp_DE_evaluate_h_

#include <RcppArmadillo.h>

namespace Rcpp {

  class EvalBase {
  public:
    EvalBase() : neval(0) {};
    virtual Rcpp::NumericVector eval(SEXP par) = 0;
    unsigned long getNbEvals() { return neval; }
  protected:
    unsigned long int neval;
  };

  class EvalStandard : public EvalBase {
  public:
    EvalStandard(SEXP fcall_, SEXP env_) : fcall(fcall_), env(env_) {} 
    Rcpp::NumericVector eval(SEXP par) {
      neval++;
      return defaultfun(par);
    }
  private:
    SEXP fcall, env;
    Rcpp::NumericVector defaultfun(SEXP par) { 
      SEXP fn = ::Rf_lang2(fcall, par);
      SEXP sexp_fvec = ::Rf_eval(fn, env);
      Rcpp::NumericVector f_result = (Rcpp::NumericVector) Rcpp::as<Rcpp::NumericVector>(sexp_fvec);
      return(f_result); 
    }
  };

  typedef Rcpp::NumericVector (*funcPtr)(SEXP, SEXP);

  class EvalCompiled : public EvalBase {
  public:
    EvalCompiled(Rcpp::XPtr<funcPtr> xptr, SEXP __env) {
      funptr = *(xptr);
      env = __env;
    };
    EvalCompiled(SEXP xps, SEXP __env) {
      Rcpp::XPtr<funcPtr> xptr(xps);
      funptr = *(xptr);
      env = __env;
    };
    Rcpp::NumericVector eval(SEXP par) {
      neval++;
      return funptr(par, env);
    }
  private:
    funcPtr funptr;
    SEXP env;
  };

}

#endif
