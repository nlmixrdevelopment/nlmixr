#ifndef __RxODE_types__
#define __RxODE_types__
#define BOOST_DISABLE_ASSERTS true
#define NDEBUG
#include <RcppArmadillo.h>
typedef Rcpp::NumericVector (*rxFn2)(SEXP,SEXP);
#endif // __RxODE_types__

