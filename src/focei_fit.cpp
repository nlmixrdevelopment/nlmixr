// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

extern "C" double RxODE_sum (double *input, int n);

//[[Rcpp::export]]
arma::mat sFOCEi(NumericVector par, Environment e){
  Function objf = as<Function>(e["ofv.FOCEi.ind"]);
  List llikSubj = objf(par);
  int nsub = llikSubj.size();
  int n = (as<NumericVector>((as<NumericVector>(llikSubj[0])).attr("grad"))).size();
  int i, j, k, l = 0;
  double *d = Calloc(n*(n+1)/2*nsub,double);
  mat m1(1, n), m3(n, n);
  NumericVector g;
  for (i =0; i < nsub; i++){
    g = (as<NumericVector>(as<NumericVector>(llikSubj[i])).attr("grad"));
    for (j = 0; j < n; j++){
      m1(0, j) = g[j];
    }
    // Cross product matrix
    m3 = m1.t() * m1;
    l = 0;
    for (j = 0; j < n; j++){
      for (k = 0; k <= j; k++){
	// Defered sum vector
	d[i+l*nsub] = m3(j,k);
	l++;
      }
    }
  }
  l = 0;
  // Summed matrix
  l = 0;
  for (j = 0; j < n; j++){
    for (k = 0; k <= j; k++){
      m3(j,k) = RxODE_sum(d+l*nsub,nsub);
      m3(k,j) = m3(j, k);
      l++;
    }
  }
  Free(d);
  // m3 is S matrix
  return m3;
}

//[[Rcpp::export]]
NumericVector grFOCEi(NumericVector par, Environment e){
  Function objf = as<Function>(e["ofv.FOCEi.ind"]);
  Function optimObj = as<Function>(e["optim.obj"]);
  List con = as<List>(e["con"]);
  double ridgeDecay = as<double>(con["ridge.decay"]);
  List llikSubj = objf(par);
  int i, j;
  int nsub = llikSubj.size();
  int n = (as<NumericVector>((as<NumericVector>(llikSubj[0])).attr("grad"))).size();
  double *d = Calloc(nsub, double);
  double *dg = Calloc(nsub*n, double);
  NumericVector g, ret(n);
  for (i =0; i < nsub; i++){
    g = llikSubj[i];
    d[i] = as<double>(g);
    g = (as<NumericVector>(as<NumericVector>(g)).attr("grad"));
    for (j = 0; j < n; j++){
      dg[i+j*nsub] = g[j];
    }
  }
  double llik = -2 * RxODE_sum(d, nsub);
  Free(d);
  double extra = 1.0;
  double precision = as<double>(con["precision"]);
  Nullable<NumericVector> curDiff = e["cur.diff"];
  if (ridgeDecay != 0.0 && curDiff.isNull()){
    curDiff = 0;
  } else if (!curDiff.isNull()){
    // extra <- exp(-con$ridge.decay * cur.diff);
    extra = REAL(as<SEXP>(curDiff))[0];
    extra = exp(-ridgeDecay * extra);
  }
  //gr <- -2 * foceiGrad(llik.subj) + pars * con$precision  * extra;
  for (j = 0; j < n; j++){
    ret[j] = -2 * RxODE_sum(dg+j*nsub,nsub) + par[j]*extra*precision;
  }
  Free(dg);
  // assign(optim.obj(llik, "l"), gr, envir=ofv.cache, inherits=FALSE);
  Environment ofvCache = as<Environment>(e["ofv.cache"]);
  std::string lStr = "l";
  std::string llikO =as<std::string>(optimObj(llik,lStr));
  ofvCache[llikO] = ret;
  return ret;
}
