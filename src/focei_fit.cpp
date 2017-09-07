#include <Rcpp.h>
using namespace Rcpp;


typedef double (*rxode_sum_fn)(double *input, int n);
// Returns -2*grad
//[[Rcpp::export]]
NumericVector foceiGrad(List objf){
  static rxode_sum_fn RxODE_sum = NULL;
  if (RxODE_sum == NULL) RxODE_sum = (rxode_sum_fn) R_GetCCallable("RxODE","RxODE_sum");
  
  int nsub = objf.size();
  int i, j;
  int n = (as<NumericVector>((as<NumericVector>(objf[0])).attr("grad"))).size();
  NumericVector g, ret(n);
  double *d = Calloc(nsub*n, double);
  for (i = 0; i < nsub; i++){
    g = (as<NumericVector>(as<NumericVector>(objf[i])).attr("grad"));
    for (j = 0; j < n; j++){
      d[i+j*nsub] = g[j];
    }
  }
  for (j = 0; j < n; j++){
    ret[j] = RxODE_sum(d+j*nsub,nsub);
  }
  Free(d);
  return ret;
}
