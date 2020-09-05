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

typedef void (*sortIds_t)(rx_solve* rx, int ini);

sortIds_t sortIds = (sortIds_t) R_GetCCallable("RxODE", "sortIds");

typedef t_update_inis (*getUpdateInis_t)();
getUpdateInis_t getUpdateInis = (getUpdateInis_t) R_GetCCallable("RxODE", "getUpdateInis");

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

mat user_function(const mat &_phi, const mat &_evt, const List &_opt) {
  // yp has all the observations in the dataset
  rx_solving_options_ind *ind;
  rx_solving_options *op = _rx->op;
  vec _id = _evt.col(0);
  int _N=_id.max()+1;
  SEXP paramUpdate = _opt["paramUpdate"];
  int *doParam = INTEGER(paramUpdate);
  int nPar = Rf_length(paramUpdate);
  // Fill in subject parameter information
  for (int _i = 0; _i < _N; ++_i) {
    ind = &(_rx->subjects[_i]);
    ind->solved = -1;
    // ind->par_ptr
    int k=0;
    for (int _j = 0; _j < nPar; _j++){
      if (doParam[_j] == 1) {
	ind->par_ptr[_j] = _phi(_i, k++);
      }
    }
  }
  saem_solve(_rx); // Solve the complete system (possibly in parallel)
  mat g(_rx->nobs2, 3); // nobs EXCLUDING EVID=2
  int elt=0;
  for (int id = 0; id < _N; ++id) {
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
	g(elt, 0) = ind->lhs[0];
	if (_rx->cens) {
	  g(elt, 1) = ind->cens[j];
	} else {
	  g(elt, 1) = 0;
	}
	if (_rx->limit) {
	  g(elt, 2) = ind->limit[j];
	} else {
	  g(elt, 2) = R_NegInf;
	}
	elt++;
      } // evid=2 does not need to be calculated
    }
  }
  if (op->stiff == 2) { // liblsoda
    // Order by the overall solve time
    // Should it be done every time? Every x times?
    sortIds(_rx, 0);
  }
  return g;
}

void setupRx(List &opt, SEXP evt, SEXP evtM) {
  RObject obj = opt[".rx"];
  if (!Rf_isNull(obj)){
    // Now need to get the largest item to setup the solving space
    // foceiSetupEta_(etaMat0);
    RObject pars = opt[".pars"];
    List odeO = opt["ODEopt"];
    // SEXP evt = x["evt"];
    // SEXP evtM = x["evtM"];
    int nEvt = INTEGER(Rf_getAttrib(evt, R_DimSymbol))[0];
    int nEvtM = INTEGER(Rf_getAttrib(evtM, R_DimSymbol))[0];
    SEXP ev;
    if (nEvt > nEvtM) {
      ev = evt;
    } else {
      ev = evtM;
    }
    RxODE::rxSolve_(obj, odeO,
     		    R_NilValue,//const Nullable<CharacterVector> &specParams =
     		    R_NilValue,//const Nullable<List> &extraArgs =
     		    pars,//const RObject &params =
     		    ev,//const RObject &events =
     		    R_NilValue, // inits
     		    1);//const int setupOnly = 0
  } else {
    stop("cannot find RxODE model");
  }
}

//[[Rcpp::export]]
SEXP saem_do_pred(SEXP in_phi, SEXP in_evt, SEXP in_opt) {
  List opt = List(in_opt);
  setupRx(opt, in_evt, in_evt);
  saem_lhs = getRxLhs();
  saem_inis = getUpdateInis();
  mat phi = as<mat>(in_phi);
  mat evt = as<mat>(in_evt);
  int distribution = as<int>(opt["distribution"]);
  mat gMat = user_function(phi, evt, opt);
  vec g = gMat.col(0);
  if (distribution == 4) g = log(g);
  return wrap(g);
}


//[[Rcpp::export]]
SEXP saem_fit(SEXP xSEXP) {
  if (getRx_ == NULL) getRx_ = (getRxSolve_t) R_GetCCallable("RxODE","getRxSolve_");
  List x(xSEXP);
  List opt = x["opt"];
  setupRx(opt,x["evt"],x["evtM"]);
  // if (rxSingleSolve == NULL) rxSingleSolve = (rxSingleSolve_t) R_GetCCallable("RxODE","rxSingleSolve");
  saem_lhs = getRxLhs();
  saem_inis = getUpdateInis();
  _rx=getRx_();

  SAEM saem;
  saem.inits(x);
  saem.set_fn(user_function);
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
    Named("par_hist") = saem.get_par_hist(),
    Named("res_info") = saem.get_resInfo()
  );
  out.attr("saem.cfg") = x;
  out.attr("class") = "saemFit";
  return out;
}
