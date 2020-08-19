#include <stdio.h>
#include <stdarg.h>
#include <thread>
#include <chrono>
#include <R_ext/Rdynload.h>
#include <RcppArmadillo.h>
#include <RxODE.h>
#include "saem_class_rcpp.h"
using namespace Rcpp;

typedef rx_solve *(*getRxSolve_t)();
getRxSolve_t getRx_ = (getRxSolve_t) R_GetCCallable("RxODE","getRxSolve_");

typedef t_calc_lhs (*getRxLhs_t)();

getRxLhs_t getRxLhs = (getRxLhs_t) R_GetCCallable("RxODE","getRxLhs");


typedef t_update_inis (*getUpdateInis_t)();
getUpdateInis_t getUpdateInis = (getUpdateInis_t) R_getCCallable("RxODE", "getUpdateInis");

t_calc_lhs saem_lhs = NULL;
t_update_inis saem_inis = NULL;

typedef void (*par_solve_t)(rx_solve *rx);

par_solve_t saem_solve = (par_solve_t) R_GetCCallable("RxODE","par_solve");

typedef int (*iniSubjectE_t)(int solveid, int inLhs, rx_solving_options_ind *ind, rx_solving_options *op, rx_solve *rx,
			     t_update_inis u_inis);

iniSubjectE_t iniSubjectE = (iniSubjectE_t) R_GetCCallable("RxODE","iniSubjectE");

rx_solve* _rx = NULL;

RObject mat2NumMat(const mat &m) {
  RObject x = wrap( m.memptr() , m.memptr() + m.n_elem ) ;
  x.attr( "dim" ) = Dimension( m.n_rows, m.n_cols ) ;
  return x;
}

// FIXME use C to calculate sd
Function ff("sd");

vec Ruser_function(const mat &phi_, const mat &evt_, const List &opt) {
  RObject phi, evt;
  phi = mat2NumMat(phi_);
  evt = mat2NumMat(evt_);
  NumericVector g;
  g = ff(phi, evt);
  vec yp(g);
  return yp;
}

vec user_function(const mat &_phi, const mat &_evt, const List &_opt){
  // yp has all the observations in the dataset
  rx_solving_options_ind *ind;
  rx_solving_options *op = _rx->op;
  saem_solve(_rx); // Solve the complete system (possibly in parallel)
  vec g(_rx->nobs2); // nobs EXCLUDING EVID=2
  int elt = 0;
  for (int id = 0; id < _rx->nsub; ++id) {
    ind = &(_rx->subjects[id]);
    iniSubjectE(op->neq, 1, ind, op, _rx, saem_inis);
    for (int j = 0; j < ind->n_all_times; ++j){
      ind->idx=j;
      if (isDose(ind->evid[j])){
	ind->tlast = ind->all_times[j];
	// Need to calculate for advan sensitivities
	saem_lhs((int)id, ind->all_times[j],
		 getSolve(j), ind->lhs);
      } else if (ind->evid[j] == 0) {
	saem_lhs((int)id, ind->all_times[j],
		 getSolve(j), ind->lhs);
	g[elt++] = ind->lhs[0];
      } // evid=2 does not need to be calculated
    }
  }
  return g;
}

vec _do_pred(SEXP in_phi, SEXP in_evt, SEXP in_opt) {
  saem_lhs = getRxLhs();
  saem_inis = getUpdateInis();
  mat phi = as<mat>(in_phi);
  mat evt = as<mat>(in_evt);
  List opt= as<List>(in_opt);
  int distribution = as<int>(opt["distribution"]);
  vec g = user_function(phi, evt, opt);
  if (distribution == 4) g = log(g);
  return g;
}

SEXP _saem_fit(SEXP xSEXP) {
  if (getRx_ == NULL) getRx_ = (getRxSolve_t) R_GetCCallable("RxODE","getRxSolve_");
  // if (rxSingleSolve == NULL) rxSingleSolve = (rxSingleSolve_t) R_GetCCallable("RxODE","rxSingleSolve");
  saem_lhs = getRxLhs();
  saem_inis = getUpdateInis();
  _rx=getRx_();
  List x(xSEXP);
  getRxLhs

  SAEM saem;
  saem.inits(x);

  if(x.containsElementNamed("Rfn")) {
    ff = as<Function>(x["Rfn"]);
    saem.set_fn(Ruser_function);
  } else {
    saem.set_fn(user_function);
  }

  saem.saem_fit();

  List out = List::create(
    Named("resMat") = saem.get_resMat(),
    Named("mprior_phi") = saem.get_mprior_phi(),
    Named("mpost_phi") = saem.get_mpost_phi(),
    Named("Gamma2_phi1") = saem.get_Gamma2_phi1(),
    Named("Plambda") = saem.get_Plambda(),
    Named("Ha") = saem.get_Ha(),
    Named("sig2") = saem.get_sig2(),
    Named("eta") = saem.get_eta(),
    Named("par_hist") = saem.get_par_hist()
  );
  out.attr("saem.cfg") = x;
  out.attr("class") = "saemFit";
  return out;
}

// This is a thin layer to the dll.  This is to try to avoid corruption

typedef SEXP (*t_dopred)(SEXP in_phi, SEXP in_evt, SEXP in_opt);
typedef SEXP (*t_saem_fit)(SEXP xSEXP);

t_dopred dopred = NULL;
t_saem_fit saem_fit = NULL;

void setCur(std::string cur, std::string dll){
  if (cur == "") stop("Need to specify a saem dll string.");
  Function isLoaded("is.loaded", R_BaseNamespace);
  std::string saemPredStr = "_" + cur + "_dopred";
  std::string saemFitStr = "_" + cur + "_saem_fit";
  std::string saemCallStr = "_" + cur + "_call";
  if (!isLoaded(wrap(saemCallStr))){
    for (int j = 4; j--;){
      RxODE::dynLoad(dll);
      if (isLoaded(wrap(saemCallStr))){
	break;
      } else {
	std::this_thread::sleep_for(std::chrono::milliseconds(250));
      }
    }
    if (!isLoaded(wrap(saemCallStr))){
      stop("saem dll not loaded.");
    }
  }
  dopred =(t_dopred) R_GetCCallable(cur.c_str(), saemPredStr.c_str());
  saem_fit =(t_saem_fit) R_GetCCallable(cur.c_str(), saemFitStr.c_str());
}

//[[Rcpp::export]]
RObject saemDoPred(RObject in_phi, RObject in_evt, RObject in_opt,
		   CharacterVector cur, RObject model,
		   std::string dll){
  if (!RxODE::rxIs(model, "NULL")){
    if (!RxODE::rxDynLoad(model)){
      stop("Cannot load RxODE dlls for the SAEM prediction");
    }
  }
  setCur(as<std::string>(cur[0]), dll);
  return dopred(in_phi, clone(in_evt), in_opt);
}

//[[Rcpp::export]]
RObject saemFit(RObject xSEXP, CharacterVector cur, RObject model, std::string dll){
  if (!RxODE::rxIs(model, "NULL")){
    if (!RxODE::rxDynLoad(model)){
      stop("Cannot load RxODE dlls for the SAEM fit");
    }
  }
  setCur(as<std::string>(cur[0]),dll);
  return saem_fit(xSEXP);
}

