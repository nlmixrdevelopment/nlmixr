## saem_fit.R: population PK/PD modeling library
##
## Copyright (C) 2014 - 2016  Wenping Wang
##
## This file is part of nlmixr.
##
## nlmixr is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## nlmixr is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with nlmixr.  If not, see <http://www.gnu.org/licenses/>.

saem_ode_str <- '#define ARMA_DONT_PRINT_ERRORS
#define ARMA_DONT_USE_OPENMP // Known to cause speed problems
#include <RcppArmadillo.h>
#include <RxODE.h>
#include "saem_class_rcpp.hpp"
//#include <omp.h>


using namespace std;
using namespace arma;

extern "C" {

typedef void (*rxSingleSolve_t)(int subid, double *_theta, double *timep,
			  int *evidp, int *ntime,
			  double *initsp, double *dosep,
			  double *ii, double *retp,
			  double *lhsp, int *rc,
			  double *newTime, int *newEvid,
			  int *on, int *ix,
			  int *slvr_counter, int *dadt_counter, int *jac_counter,
			  double *InfusionRate, int *BadDose, int *idose,
			  double *scale, int *stateIgnore, double *mtime);

rxSingleSolve_t rxSingleSolve = (rxSingleSolve_t) R_GetCCallable("RxODE","rxSingleSolve");

typedef rx_solve *(*getRxSolve_t)();
getRxSolve_t getRx = (getRxSolve_t) R_GetCCallable("RxODE","getRxSolve_");

rx_solve* _rx = NULL;

}

Function ff("sd");

RObject mat2NumMat(const mat &m) {
	RObject x = wrap( m.memptr() , m.memptr() + m.n_elem ) ;
	x.attr( "dim" ) = Dimension( m.n_rows, m.n_cols ) ;
	return x;
}

vec Ruser_function(const mat &phi_, const mat &evt_, const List &opt) {
  RObject phi, evt;
  phi = mat2NumMat(phi_);
  evt = mat2NumMat(evt_);
  NumericVector g;
  g = ff(phi, evt);
  vec yp(g);
  return yp;
}


vec user_function(const mat &_phi, const mat &_evt, const List &_opt) {
  rx_solving_options* _op = _rx->op;
  vec _id = _evt.col(0);
  int _N=_id.max()+1;
  int _cores = 1;//_op->cores;
  uvec _ix;
  _ix = find(_evt.col(2) == 0);
  vec _yp(_ix.n_elem);
  double *_p=_yp.memptr();
  vec _id0 = _id(_ix);
  int _DEBUG = _opt["DEBUG"];
  uvec _cmt_endpnt = _opt["cmt_endpnt"];
  for (int _i=0; _i<_N; _i++) {
     int _nlhs = _op->nlhs;
     vec _inits(_op->neq, fill::zeros);
     vec _scale(_op->neq, fill::ones);
     ivec _stateIgnore(_op->neq, fill::zeros);

     mat _wm;
     vec _wv, _wv2;

     //////////////////////////////////////////////////////////////////////
     // declPars
<%=declPars%>
     int _slvr_counter=0, _dadt_counter=0, _jac_counter=0;

     //////////////////////////////////////////////////////////////////////
     // assgnPars
<%=assgnPars%>

     vec _InfusionRate(_op->neq, fill::zeros);
     ivec _BadDose(_op->neq, fill::zeros);
     ivec _on(_op->neq, fill::ones);
    _wm = _evt.rows( find(_id == _i) );
    if(_wm.n_rows==0) {
        Rcout << "ID = " << _i+1 << " has no data. Please check." << endl;
        arma_stop_runtime_error("");
    }
    vec _time__;
    _time__ = _wm.col(1);
    int _ntime = _time__.n_elem;
    vec _newTime(_ntime);
    _wv = _wm.col(2);
    ivec _evid(_ntime);
    ivec _evid2(_ntime);
    _wv2 = _wm.col(5);
    ivec _cmt(_ntime);
    for (int _k=_ntime; _k--;){
      _evid(_k) = _wv(_k);
      _cmt(_k) = _wv2(_k);
    }
    _wv = _wm.col(3);
    uvec _ds  = find(_evid > 99 || _evid == 3);
    vec _amt;
    _amt = _wv(_ds);
    ivec _idose(_amt.n_elem);
    _wv = _wm.col(4);
    vec _ii;
    _ii = _wv(_ds);
    ivec _ix2(_ntime);
    for (int _jj = _ntime; _jj--;) _ix2[_jj] = _jj;
    //std::iota(_ix2.memptr(),_ix2.memptr()+_ntime, 0); // 0, 1, 2, 3...
    vec _mtime(<%=nmtime%>, fill::zeros);

    // _inits.zeros();
    //std::copy(&_op->inits[0], &_op->inits[0]+_op->neq, &_inits[0]);
    //_inits.zeros(); //as<vec>(_opt["inits"]);	//FIXME

<%=foo%>

<%=pars%>
<%=inits%>


    int _rc=0;

    mat _ret(_op->neq, _ntime);
    mat _lhs(_nlhs, _ntime);

    rxSingleSolve(_i, _params.memptr(), _time__.memptr(),
	    _evid.memptr(), &_ntime, _inits.memptr(), _amt.memptr(), _ii.memptr(),
            _ret.memptr(), _lhs.memptr(), &_rc, _newTime.memptr(), _evid2.memptr(),
            _on.memptr(), _ix2.memptr(), &_slvr_counter,&_dadt_counter, &_jac_counter,
            _InfusionRate.memptr(), _BadDose.memptr(), _idose.memptr(), _scale.memptr(),
            _stateIgnore.memptr(), _mtime.memptr());

    if ( _DEBUG > 4 && _rc != 0 ) {

        Rcout << "pars: " << _params.t();
        Rcout << "_inits: " << _inits.t();
        Rcout << "LSODA return code: " << _rc << endl;
        Rcout << _wm << endl;
    }
	_ret = join_cols(join_cols(_newTime.t(), _ret), _lhs).t();
	uvec _r  = find(_evid2 == 0);
	_ret = _ret.rows(_r);
	ivec _cmtObs = _cmt(find(_evid==0));

<%=model_vars_decl%>


mat _g(time.n_elem, <%=nendpnt%>);
<%=pred_expr%>

if (_g.has_nan()) {
	Rcout << "NaN in prediction. Consider to: relax atol & rtol; change initials; change seed; change structure model." << endl;
    if ( _DEBUG > 4) {
	Rcout << "pars: " << _params.t();
	Rcout << "inits: " << _inits.t();
	Rcout << "LSODA code: " << _rc << endl;
	Rcout << "input data:" << endl;
	Rcout << _wm;
	Rcout << "LSODA solutions:" << endl;
	Rcout << _ret << endl;
	}
	_g.replace(datum::nan, 1.0e99);
}

uvec _b0(1), _b1(1); _b0(0) = 0;

for (int _b=1; _b< <%=nendpnt%>; ++_b) {
  _b1(0) = _b;
  uvec _r;
  _r = find( _cmtObs==_cmt_endpnt(_b) );
  _g.submat(_r, _b0) = _g.submat(_r, _b1);
}


    int _no = _cmtObs.n_elem;
    memcpy(_p, _g.memptr(), _no*sizeof(double));
    _p += _no;
  }
  return _yp;
}

// definition
SEXP _<%=saem.base%>_dopred( SEXP in_phi, SEXP in_evt, SEXP in_opt){
    if (getRx == NULL) getRx = (getRxSolve_t) R_GetCCallable("RxODE","getRxSolve_");
    if (rxSingleSolve == NULL) rxSingleSolve = (rxSingleSolve_t) R_GetCCallable("RxODE","rxSingleSolve");
    _rx=getRx();
    mat phi = as<mat>(in_phi);
    mat evt = as<mat>(in_evt);
    List opt= as<List>(in_opt);
    int distribution = as<int>(opt["distribution"]);
    vec g = user_function(phi, evt, opt);
    if (distribution == 4) g = log(g);
    return Rcpp::wrap(g);
}

SEXP _<%=saem.base%>_saem_fit(SEXP xSEXP) {
  if (getRx == NULL) getRx = (getRxSolve_t) R_GetCCallable("RxODE","getRxSolve_");
  if (rxSingleSolve == NULL) rxSingleSolve = (rxSingleSolve_t) R_GetCCallable("RxODE","rxSingleSolve");
  _rx=getRx();
  List x(xSEXP);

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
extern "C" {
SEXP _<%=saem.base%>_call(){
  return R_NilValue;
}
void R_init_<%=saem.base%>(DllInfo *dll)
{
  R_RegisterCCallable("<%=saem.base%>","_<%=saem.base%>_saem_fit",(DL_FUNC) &_<%=saem.base%>_saem_fit);
  R_RegisterCCallable("<%=saem.base%>","_<%=saem.base%>_dopred",(DL_FUNC) &_<%=saem.base%>_dopred);
  static const R_CallMethodDef callMethods[]  = {
    {"_<%=saem.base%>_call", (DL_FUNC) &_<%=saem.base%>_call, 0},
    {NULL, NULL, 0}
  };
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
}
'


saem_cmt_str <- '#include <RcppArmadillo.h>
#include <R_ext/Rdynload.h>
#include <Eigen/Dense>
#include "saem_class_rcpp.hpp"
#include "lin_cmt.hpp"


using namespace std;
using namespace arma;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

Function ff("sd");

RObject mat2NumMat(const mat &m) {
	RObject x = wrap( m.memptr() , m.memptr() + m.n_elem ) ;
	x.attr( "dim" ) = Dimension( m.n_rows, m.n_cols ) ;
	return x;
}

vec Ruser_function(const mat &phi_, const mat &evt_, const List &opt) {
  RObject phi, evt;
  phi = mat2NumMat(phi_);
  evt = mat2NumMat(evt_);
  NumericVector g;
  g = ff(phi, evt);
  vec yp(g);

  return yp;
}


vec user_function(const mat &_phi, const mat &_evt, const List &_opt) {
  uvec _ix;
  vec _id = _evt.col(0);
  mat _wm;
  vec _obs_time, _dose_time, _dose, _wv;

  _ix = find(_evt.col(2) == 0);
  vec _yp(_ix.n_elem);
  double *_p=_yp.memptr();
  int _N=_id.max()+1;

  for (int _i=0; _i<_N; _i++) {
    _ix = find(_id == _i);
    _wm = _evt.rows(_ix);

    _ix = find(_wm.col(2) == 0);
    _wv = _wm.col(1);
    _wv = _wv(_ix);
    const Map<MatrixXd> __obs_time(_wv.memptr(), _wv.n_elem, 1);
    const VectorXd _obs_time(__obs_time);

    _ix = find(_wm.col(2) > 0);
    _wv = _wm.col(1);
    _wv = _wv(_ix);
    const Map<MatrixXd> __dose_time(_wv.memptr(), _wv.n_elem, 1);
    const VectorXd _dose_time(__dose_time);

    _wv = _wm.col(3);
    _wv = _wv(_ix);
    const Map<MatrixXd> __dose(_wv.memptr(), _wv.n_elem, 1);
    const VectorXd _dose(__dose);

    _wv = _wm.col(4);
    _wv = _wv(_ix);
    const Map<MatrixXd> __Tinf(_wv.memptr(), _wv.n_elem, 1);
    const VectorXd _Tinf(__Tinf);

<%=foo%>
<%=pars%>

    int _no=_obs_time.size();
    VectorXd _g(_obs_time.size());
    int _ncmt=<%=ncmt%>, _oral=<%=oral%>, _infusion=<%=infusion%>, _parameterization=<%=parameterization%>;

	_g = generic_cmt_interface(
      _obs_time,
      _dose_time,
      _dose,
      _Tinf,
      _params,
      _ncmt,
      _oral,
      _infusion,
      _parameterization);

    memcpy(_p, _g.data(), _no*sizeof(double));
    _p += _no;
	//cout << "ok " << _i <<endl;
  }

  return _yp;
}

// definition
extern "C"  SEXP _<%=saem.base%>_dopred( SEXP in_phi, SEXP in_evt, SEXP in_opt ) {
    Rcpp::traits::input_parameter< mat& >::type phi(in_phi);
    Rcpp::traits::input_parameter< mat& >::type evt(in_evt);
    List opt(in_opt);
    int distribution = as<int>(opt["distribution"]);
    vec g = user_function(phi, evt, opt);
    if (distribution == 4) g = log(g);
    return Rcpp::wrap(g);
}

extern "C" SEXP _<%=saem.base%>_saem_fit(SEXP xSEXP) {
  List x(xSEXP);
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

extern "C" {
SEXP _<%=saem.base%>_call(){
  return R_NilValue;
}
void R_init_<%=saem.base%>(DllInfo *dll)
{
  R_RegisterCCallable("<%=saem.base%>","_<%=saem.base%>_saem_fit",(DL_FUNC) &_<%=saem.base%>_saem_fit);
R_RegisterCCallable("<%=saem.base%>","_<%=saem.base%>_dopred",(DL_FUNC) &_<%=saem.base%>_dopred);
  static const R_CallMethodDef callMethods[]  = {
    {"_<%=saem.base%>_call", (DL_FUNC) &_<%=saem.base%>_call, 0},
    {NULL, NULL, 0}
  };
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
}
'


##' Get a list of directories for inclusion
##'
##' @param pkg a string or list of string for packages to be included.
##' @return An inclusion string for Makevars
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
nmxInclude <- function(pkg = "nlmixr") {
  if (length(pkg) == 1) {
    x <- system.file("", package = pkg)
    if (.Platform$OS.type == "windows") x <- gsub("\\\\", "/", utils::shortPathName(x))
    x <- paste0("-I", x, "/include")
    return(x)
  } else {
    paste(sapply(pkg, nmxInclude), collapse = " ")
  }
}

..saemCountDll <- list()
.saemCountDll <- function(dll, inc = NULL) {
  .what <- ..saemCountDll[[dll]]
  if (is.null(.what)) {
    if (is.integer(inc)) {
      .lst <- ..saemCountDll
      .lst[[dll]] <- inc
      assignInMyNamespace("..saemCountDll", .lst)
      return(inc)
    } else {
      return(0L)
    }
  } else if (is.null(inc)) {
    return(.what)
  } else {
    .lst <- ..saemCountDll
    .what <- .what + inc
    .lst[[dll]] <- .what
    assignInMyNamespace("..saemCountDll", .lst)
    return(.what)
  }
}

#' Generate an SAEM model
#'
#' Generate an SAEM model using either closed-form solutions or ODE-based model definitions
#'
#' @param model a compiled SAEM model by gen_saem_user_fn()
#' @param PKpars PKpars function
#' @param pred  pred function
#' @param inPars a character vector of parameters to be read from the input dataset
#' @details
#'    Fit a generalized nonlinear mixed-effect model using the Stochastic
#'    Approximation Expectation-Maximization (SAEM) algorithm
#'
#' @author Wenping Wang
#' @export
gen_saem_user_fn <- function(model, PKpars = attr(model, "default.pars"), pred = NULL, inPars = NULL) {
  is.ode <- class(model) == "RxODE"
  is.win <- .Platform$OS.type == "windows"
  env <- environment()
  lwd <- getwd()
  .md5 <- digest::digest(list(
    ifelse(is.ode, RxODE::rxModelVars(model)$trans["lib.name"],
      ## lib.name includes RxODE version MD5
      deparse(model)
    ),
    ## Should give different models for different nlmixr versions
    sessionInfo()$otherPkgs$nlmixr$Version,
    deparse(body(PKpars)), deparse(body(pred)),
    deparse(inPars)
  ))
  if (getOption("RxODE.tempfiles", TRUE)) {
    ## .wd <- tempfile()
    ## dir.create(.wd, recursive = TRUE)
    ## setwd(.wd)
    .wd <- RxODE::rxTempDir()
    if (.wd == "") {
      .wd <- tempfile()
      dir.create(.wd, recursive = TRUE)
      setwd(.wd)
      on.exit({
        setwd(lwd)
        unlink(.wd, recursive = TRUE, force = TRUE)
      })
    } else {
      ## This makes all the parsing files, cpp and so in their own
      ## directory.  No collisions.
      .wd <- file.path(.wd, paste0("saem-", .md5, ".saemd"))
      suppressWarnings({
        dir.create(.wd, recursive = TRUE)
      })
      setwd(.wd)
      on.exit({
        setwd(lwd)
      })
    }
  }

  saem.cpp <- paste0("saem", .md5, .Platform$r_arch)
  ## }
  saem.base <- saem.cpp
  saem.dll <- paste0(saem.cpp, .Platform$dynlib.ext)
  saem.cpp <- paste0(saem.cpp, ".cpp")
  saem.lock <- paste0(saem.cpp, ".lock")

  if (!file.exists("saem.cpp")) {
    if (is.ode) {
      modelVars <- model$cmpMgr$get.modelVars()
      pars <- modelVars$params
      pars <- gsub("[.]", "_DoT_", pars)
      npar <- length(pars)
      pars <- paste(c(
        sprintf("    vec _params(%d);\n", npar),
        sprintf("    _params(%d) = %s;\n", 1:npar - 1, pars)
      ), collapse = "")

      model_vars <- names <- c("time", modelVars$state, modelVars$lhs)
      s <- lapply(1:length(model_vars), function(k) {
        sprintf("vec %s;\n%s=_ret.col(%d);\n", model_vars[k], model_vars[k], k - 1)
      })
      model_vars_decl <- paste0(s, collapse = "")

      ode_solver <- model$cmpMgr$ode_solver
      dll <- sub("[.].*", "", basename(RxODE::rxDll(model)))

      neq <- length(modelVars$state)
      nlhs <- length(modelVars$lhs)

      x <- deparse(body(pred))
      len <- length(x)
      x <- if (x[1] == "{") x[2:(len - 1)] else x
      len <- length(x)
      nendpnt <- len
      pred_expr <- paste(paste("_g.col(", 1:len - 1, ") = ", x, ";", sep = ""), collapse = "\n")
    } else {
      neq <- nlhs <- 0
      pars <- model
      list2env(attr(pars, "calls"), envir = env)
    }
    # str(ls(envir=env))

    ### deal with explicit initCondition statement
    x <- deparse(body(PKpars))
    ix <- grep(reINITS, x, perl = T)
    if (length(ix) > 0) {
      inits <- getInits(x[ix], reINITS)
      x <- x[-ix]
    } else {
      inits <- ""
    }

    x <- gsub("[.]", "_DoT_", x)

    len <- length(x)
    .tmp <- sprintf("%s;\n", x[2:(len - 1)])
    .tmp <- gsub("([{}]);", "\\1", .tmp)
    cat(.tmp, file = "eqn__.txt")

    nrhs <- integer(1)

    if (is.null(inPars)) {
      offset <- 0L
      nignore <- 0L
      ignore_vars <- ""
    } else {
      offset <- cumsum(c(0L, nchar(inPars) + 1L))
      nignore <- length(inPars)
      ignore_vars <- paste(c(inPars, ""), collapse = ",")
    }

    RxODE::rxReq("dparser")
    x <- .C("parse_pars", "eqn__.txt", "foo__.txt", nrhs, as.integer(FALSE), ignore_vars, offset, nignore)
    nrhs <- x[[3]]
    foo <- paste(readLines("foo__.txt"), collapse = "\n")

    nm <- system.file("", package = "nlmixr")
    if (is.null(inPars)) {
      assgnPars <- declPars <- ""
    } else {
      s <- sprintf("      %s = _mPars(_i,%d);", inPars, 1:length(inPars) - 1)
      assgnPars <- paste0(s, collapse = "\n")
      s <- "mat _mPars=as<mat>(_opt[\"mPars\"]);"
      declPars <- sprintf("\tdouble %s;\n\t%s", paste0(inPars, collapse = ", "), s)
    }
    nmtime <- 0
    if (is.ode) {
      nmtime <- RxODE::rxModelVars(model)$nMtime
    }

    brew(text = c(saem_cmt_str, saem_ode_str)[1 + is.ode], output = saem.cpp)
  }
  .i <- 0
  while (file.exists(saem.lock) && .i < 100) {
    if (.i == 0) message(sprintf("Waiting for %s to be removed", saem.lock))
    .i <- .i + 1
    Sys.sleep(1)
    message(".", appendLF = FALSE)
  }
  if (.i > 0) message("")
  if (!file.exists(file.path(getwd(), saem.dll))) {
    sink(saem.lock)
    cat("")
    sink()
    on.exit(if (file.exists(saem.lock)) {
      unlink(saem.lock)
    }, add = TRUE)
    ## unlink(c("eqn__.txt", "foo__.txt"))
    ## if (inPars == "") inPars = NULL

    ## gen Markevars
    ## if(is.win) x = gsub("\\\\", "/", utils::shortPathName(x))
    ## x = sub("/nlmixr", "", x)
    ## .lib=  if(is.ode) model$cmpMgr$dllfile else ""
    ## if (is.ode && .Platform$OS.type=="windows") .lib <- gsub("\\\\", "/", utils::shortPathName(.lib));

    make_str <- "PKG_CXXFLAGS=%s\nPKG_LIBS=%s $(BLAS_LIBS) $(LAPACK_LIBS)\n"
    make_str <- sprintf(make_str, nmxInclude(c("nlmixr", "StanHeaders", "Rcpp", "RcppArmadillo", "RcppEigen", "BH", "RxODE")), "")

    cat(paste0(make_str, "\n"), file = file.path(getwd(), "Makevars"))
    ## cat(make_str)

    rexec <- paste(R.home(component = "bin"), .Platform$file.sep, "R", sep = "")

    args <- c("CMD", "SHLIB", saem.cpp, "-o", saem.dll)
    .rxBinpref <- Sys.getenv("rxBINPREF")
    if (.rxBinpref != "") {
      .oldBinpref <- Sys.getenv("BINPREF")
      Sys.setenv("BINPREF" = .rxBinpref)
      on.exit(Sys.setenv("BINPREF" = .oldBinpref), add = TRUE)
    }
    ## do.call("system", list(shlib))
    .badBuild <- function(msg, stop = TRUE) {
      message(msg)
      message(cli::rule(left = "stdout output"))
      message(paste(rawToChar(.out$stdout), sep = "\n"))
      message(cli::rule(left = "stderr output"))
      message(paste(rawToChar(.out$stderr), sep = "\n"))
      if (stop) stop(msg, call. = FALSE)
    }
    message("Building SAEM model...", appendLF = FALSE)
    .out <- sys::exec_internal(cmd = rexec, args = args, error = FALSE)
    if (!(.out$status == 0 & file.exists(saem.dll))) {
      message("error")
      unlink(saem.lock)
      .badBuild("Error building SAEM model")
    }
    message("done")
    unlink(saem.lock)
  }
  ## file.copy(file.path(.wd, saem.dll), file.path(lwd, saem.dll));
  ## file.copy(file.path(.wd, saem.cpp), file.path(lwd, saem.cpp));
  ## setwd(lwd);
  saem.dll <- file.path(getwd(), saem.dll)

  if (is.ode) {
    RxODE::rxLoad(model)
    ## .saemCountDll(RxODE::rxModelVars(model)$trans["lib.name"], 1L);
  }
  ## .saemCountDll(saem.dll, 1L);
  `.DLL` <- dyn.load(saem.dll)
  assignInMyNamespace(".protectSaemDll", saem.dll)
  on.exit(
    {
      assignInMyNamespace(".protectSaemDll", "")
    },
    add = TRUE
  )
  .mod <- model
  if (!is.ode) {
    .mod <- NULL
  }
  fn.pred <- eval(bquote(function(a, b, c) {
    if (!file.exists(.(saem.dll))) stop(sprintf("Stopping since '%s' does not exist", .(saem.dll)))
    dyn.load(.(saem.dll))
    nlmixr::.protectSaem(.(saem.dll))
    on.exit(
      {
        nlmixr::.unprotectSaem()
      },
      add = TRUE
    )
    .Call(
      `_nlmixr_saemDoPred`, a, b, c, .(saem.base),
      .(.mod), .(saem.dll)
    )
  }))
  fn1 <- eval(bquote(function(a) {
    dyn.load(.(saem.dll))
    nlmixr::.protectSaem(.(saem.dll))
    on.exit(
      {
        nlmixr::.unprotectSaem()
      },
      add = TRUE
    )
    if (.(is.ode)) {
      RxODE::rxLoad(.(model))
      RxODE::rxLock(.(model))
      RxODE::rxAllowUnload(FALSE)
      on.exit({
        RxODE::rxUnlock(.(model))
      })
      .l1 <- length(unique(a$evt[, "ID"]))
      .l2 <- length(unique(a$evtM[, "ID"]))
      if (.l2 > .l1) {
        suppressWarnings(do.call(
          RxODE::rxSolve,
          c(
            list(
              object = .(model), params = a$opt$.pars,
              events = a$evtM, .setupOnly = 1L
            ),
            a$optM
          )
        ))
      } else {
        suppressWarnings(do.call(
          RxODE::rxSolve,
          c(
            list(
              object = .(model), params = a$opt$.pars,
              events = a$evt, .setupOnly = 1L
            ),
            a$optM
          )
        ))
      }
    }
    .Call(
      `_nlmixr_saemFit`, a, .(saem.base),
      .(.mod), .(saem.dll)
    )
  }))
  if (is.ode) {
    fn <- eval(bquote(function(a, b, c) {
      RxODE::rxLoad(.(model))
      RxODE::rxLock(.(model))
      on.exit({
        RxODE::rxUnlock(.(model))
        RxODE::rxAllowUnload(TRUE)
      })
      if (missing(b) && missing(c)) {
        cur.fn <- .(fn1)
        ret <- cur.fn(a)
        attr(ret, "dopred") <- .(fn.pred)
        attr(ret, "env") <- .(env)
        return(ret)
      } else {
        cur.fn <- .(fn.pred)
        return(cur.fn(a, b, c))
      }
    }))
  } else {
    fn <- eval(bquote(function(a, b, c) {
      if (missing(b) && missing(c)) {
        cur.fn <- .(fn1)
        ret <- cur.fn(a)
        attr(ret, "dopred") <- .(fn.pred)
        attr(ret, "env") <- .(env)
        return(ret)
      } else {
        cur.fn <- .(fn.pred)
        return(cur.fn(a, b, c))
      }
    }))
  }

  attr(fn, "form") <- if (is.ode) "ode" else "cls"
  attr(fn, "neq") <- neq
  attr(fn, "nlhs") <- nlhs
  attr(fn, "nrhs") <- nrhs
  attr(fn, "saem.dll") <- saem.dll
  attr(fn, "saem.cpp") <- saem.cpp
  attr(fn, "rx") <- if (is.ode) model else ""
  attr(fn, "inPars") <- inPars
  if (is.ode) attr(fn, "nendpnt") <- nendpnt
  reg.finalizer(env, saem.cleanup, onexit = TRUE) ## remove dlls on gc or proper exit of R.
  fn
}
.protectSaemDll <- ""
##' @title SAEM dll prodection from garbage collection
##'
##' @description
##' This protects the saem dll from being prematurely unloaded while
##' running it due to garbage collection.
##'
##' @param dll dll to protect
##'
##' @return nothing
##'
##' @export
.protectSaem <- function(dll) {
  assignInMyNamespace(".protectSaemDll", dll)
}
##' @rdname .protectSaem
##' @export
.unprotectSaem <- function(dll) {
  assignInMyNamespace(".protectSaemDll", "")
}
##' Cleanup saem_fit environment by removing dll after the object is no longer used by R.
##'
##' @param env Environment where cleanup needs to occur.
##' @author Matthew L. Fidler
##' @export
saem.cleanup <- function(env) {
  if (is(env, "nlmixr.ui.saem")) env <- as.saem(env)
  if (is(env, "saemFit")) env <- attr(env, "env")
  ## if (!any(.protectSaemDll== env$saem.dll)){
  ##     try({dyn.unload(env$saem.dll)}, silent=TRUE);
  ##     if (env$is.ode){
  ##         try({RxODE::rxUnload(env$model)}, silent=TRUE)
  ##     }
  ##     if (file.exists(env$saem.cpp))
  ##         unlink(env$saem.cpp);
  ## }
}

parfn.list <- c(
  par.1cmt.CL,
  par.1cmt.CL.oral,
  par.1cmt.CL.oral.tlag,
  par.2cmt.CL,
  par.2cmt.CL.oral,
  par.2cmt.CL.oral.tlag,
  par.3cmt.CL,
  par.3cmt.CL.oral,
  par.3cmt.CL.oral.tlag,
  par.1cmt.micro,
  par.1cmt.micro.oral,
  par.1cmt.micro.oral.tlag,
  par.2cmt.micro,
  par.2cmt.micro.oral,
  par.2cmt.micro.oral.tlag,
  par.3cmt.micro,
  par.3cmt.micro.oral,
  par.3cmt.micro.oral.tlag
)

#' Parameters for a linear compartment model for SAEM
#'
#' Parameters for a linear compartment model using closed-form solutions and the SAEM algorithm.
#'
#' @param ncmt number of compartments
#' @param oral logical, whether oral absorption is true
#' @param tlag logical, whether lag time is present
#' @param infusion logical, whether infusion is true
#' @param parameterization type of parameterization, 1=clearance/volume, 2=micro-constants
#' @return parameters for a linear compartment model
#' @author Wenping Wang
#' @export
lincmt <- function(ncmt, oral = T, tlag = F, infusion = F, parameterization = 1) {
  # ncmt=1; oral=T; tlag=F; parameterization=1
  # print(as.list(match.call()))
  if (infusion) {
    oral <- F
    tlag <- F
  }

  # master par list
  pm <- list(
    c("CL", "V", "KA", "TLAG"),
    c("CL", "V", "CLD", "VT", "KA", "TLAG"),
    c("CL", "V", "CLD", "VT", "CLD2", "VT2", "KA", "TLAG"),
    c("KE", "V", "KA", "TLAG"),
    c("KE", "V", "K12", "K21", "KA", "TLAG"),
    c("KE", "V", "K12", "K21", "K13", "K31", "KA", "TLAG")
  )
  dim(pm) <- c(3, 2)

  pars <- pm[[ncmt, parameterization]]
  if (!tlag) pars[2 * ncmt + 2] <- "0"
  if (!oral) pars <- pars[1:(2 * ncmt)]
  npar <- length(pars)
  s <- sprintf("_params(%d) = %s;", 1:npar - 1, pars)
  pars <- paste(c(sprintf("VectorXd _params(%d);", 2 * ncmt + 2), s), collapse = "\n")
  attr(pars, "calls") <- list(ncmt = ncmt, oral = oral, tlag = tlag, infusion = infusion, parameterization = parameterization)
  ix <- (parameterization - 1) * 9 + (ncmt - 1) * 3 + oral + tlag + 1
  attr(pars, "default.pars") <- parfn.list[[ix]]
  pars
}

#' Configure an SAEM model
#'
#' Configure an SAEM model by generating an input list to the SAEM model function
#'
#' @param model a compiled saem model by gen_saem_user_fn()
#' @param data input data
#' @param inits initial values
#' @param mcmc a list of various mcmc options
#' @param ODEopt optional ODE solving options
#' @param seed seed for random number generator
#' @param distribution one of c("normal","poisson","binomial", "lnorm")
#' @param fixed a character vector of fixed effect only parameters (no random effects attached) to be fixed
#' @param DEBUG Integer determining if debugging is enabled.
#' @details
#'    Fit a generalized nonlinear mixed-effect model by he Stochastic
#'    Approximation Expectation-Maximization (SAEM) algorithm
#'
#' @author Wenping Wang
#' @examples
#' \dontrun{
#' library(nlmixr)
#'
#' # ode <- "d/dt(depot) =-KA*depot;
#' #        d/dt(centr) = KA*depot - KE*centr;"
#' # m1 = RxODE(ode, modName="m1")
#' # ode <- "C2 = centr/V;
#' #      d/dt(depot) =-KA*depot;
#' #      d/dt(centr) = KA*depot - KE*centr;"
#' # m2 = RxODE(ode, modName="m2")
#'
#' # Note: only use the '=' assignment, not the '<-' at this point
#'
#' PKpars <- function() {
#'   CL <- exp(lCL)
#'   V <- exp(lV)
#'   KA <- exp(lKA)
#'   KE <- CL / V
#'   # initCondition = c(0, KE - CL/V)
#' }
#' PRED <- function() centr / V
#' PRED2 <- function() C2
#'
#' saem_fit <- gen_saem_user_fn(model = lincmt(ncmt = 1, oral = T))
#' # saem_fit <- gen_saem_user_fn(model=m1, PKpars, pred=PRED)
#' # saem_fit <- gen_saem_user_fn(model=m2, PKpars, pred=PRED2)
#'
#'
#' #--- saem cfg
#' nmdat <- theo_sd
#' model <- list(saem_mod = saem_fit, covars = "WT")
#' inits <- list(theta = c(.05, .5, 2))
#' cfg <- configsaem(model, nmdat, inits)
#' cfg$print <- 50
#'
#' # cfg$Rfn = nlmixr:::Ruser_function_cmt
#' # dyn.load("m1.d/m1.so");cfg$Rfn = nlmixr:::Ruser_function_ode
#' fit <- saem_fit(cfg)
#' df <- simple.gof(fit)
#' xyplot(DV ~ TIME | ID, df, type = c("p", "l"), lwd = c(NA, 1), pch = c(1, NA), groups = grp)
#' fit
#' }
#' @export
configsaem <- function(model, data, inits,
                       mcmc = list(niter = c(200, 300), nmc = 3, nu = c(2, 2, 2)),
                       ODEopt = list(atol = 1e-6, rtol = 1e-4, method = "lsoda", transitAbs = FALSE, maxeval = 100000),
                       distribution = c("normal", "poisson", "binomial", "lnorm"),
                       seed = 99, fixed = NULL, DEBUG = 0) {
  names(ODEopt) <- gsub("transit_abs", "transitAbs", names(ODEopt))
  ODEopt <- do.call(RxODE::rxControl, ODEopt)
  # mcmc=list(niter=c(200,300), nmc=3, nu=c(2,2,2));ODEopt = list(atol=1e-6, rtol=1e-4, stiff=1, transit_abs=0);distribution=c("normal","poisson","binomial");seed=99;data=dat;distribution=1;fixed=NULL
  set.seed(seed)
  distribution.idx <- c("normal" = 1, "poisson" = 2, "binomial" = 3, "lnorm" = 4)
  distribution <- match.arg(distribution)
  distribution <- distribution.idx[distribution]
  .extraLL <- 0
  if (distribution == 4) {
    data$DV[data$EVID == 0] <- log(data$DV[data$EVID == 0])
    if (any(is.na(data$DV[data$EVID == 0]))) stop("Data must be positive for log-normally distributed data.")
    .extraLL <- sum(data$DV[data$EVID == 0])
  }
  .data <- data
  ## RxODE::rxTrans(data, model)
  data <- list(nmdat = data)

  neq <- attr(model$saem_mod, "neq")
  nlhs <- attr(model$saem_mod, "nlhs")
  inPars <- attr(model$saem_mod, "inPars")
  ninputpars <- length(inPars)
  opt <- optM <- c(list(neq = neq, nlhs = nlhs, inits = numeric(neq)), ODEopt,
    ninputpars = ninputpars, inPars = inPars
  )

  model$N.eta <- attr(model$saem_mod, "nrhs")
  model$nendpnt <- attr(model$saem_mod, "nendpnt")
  if (is.null(model$nendpnt)) model$nendpnt <- 1

  if (is.null(model$log.eta)) model$log.eta <- rep(TRUE, model$N.eta)
  if (is.null(model$omega)) model$omega <- diag(model$N.eta)
  if (is.null(model$res.mod)) model$res.mod <- rep(1, model$nendpnt)
  if (is.null(inits$omega)) inits$omega <- rep(1, model$N.eta) * 4
  if (is.null(inits$ares)) inits$ares <- 10
  if (is.null(inits$bres)) inits$bres <- 1
  if (is.null(mcmc$print)) mcmc$print <- 1
  if (is.null(names(inits$theta))) names(inits$theta) <- rep("", length(inits$theta))
  inits.save <- inits
  inits$theta.fix <- matrix(names(inits$theta), byrow = T, ncol = model$N.eta)
  inits$theta <- matrix(inits$theta, byrow = T, ncol = model$N.eta)
  model$cov.mod <- 1 - is.na(inits$theta)
  data$N.covar <- nrow(inits$theta) - 1
  inits$theta[is.na(inits$theta)] <- 0


  ###  FIXME
  mcmc$stepsize <- 0:1
  mcmc$burn.in <- 300

  ###  FIXME: chk covars as char vec
  wh <- setdiff(c(model$covars, inPars), names(data$nmdat))
  if (length(wh)) {
    msg <- paste0("covariate(s) not found: ", paste(wh, collapse = ", "))
    stop(msg)
  }
  s <- subset(data$nmdat, EVID == 0)
  data$data <- as.matrix(s[, c("ID", "TIME", "DV", c(model$covars, inPars))])

  ###  chk for no obs records
  wh <- setdiff(unique(data$nmdat$ID), unique(data$data[, "ID"]))
  if (length(wh)) {
    msg <- paste0("No data with ID: ", paste(wh, collapse = ", "))
    stop(msg)
  }

  nphi <- model$N.eta
  mcov <- model$cov.mod
  covstruct <- model$omega

  check <- sum((covstruct - t(covstruct)) != 0)
  if (check) stop("illegal covstruct")
  check <- nphi - dim(covstruct)[1]
  if (check) stop("nphi and covstruct dim mismatch")

  check <- prod(mcov[1, ])
  if (check == 0) stop("structural par(s) absent")
  check <- nphi - dim(mcov)[2]
  if (check) stop("nphi and ncol(mcov) mismatch")
  check <- sum(dim(inits$theta) - dim(mcov) != 0)
  if (check) stop("initial theta's and mcov dim mismatch")
  check <- data$N.covar + 1 - dim(mcov)[1]
  if (check) stop("dim mcov and N.covar mismatch")

  check <- length(model$log.eta) - nphi
  if (check) stop("jlog length and nphi mismatch")

  check <- length(inits$omega) - nphi
  if (check) stop("length of omega inits and nphi mismatch")

  # check = mcmc$burn.in>sum(mcmc$niter)
  # if (check) stop("#burn-in exceeds niter")

  check <- prod(is.element(covstruct, c(0, 1)))
  if (check == 0) warning("non-zero value(s) in covstruct set to 1")
  covstruct[covstruct != 0] <- 1

  check <- prod(is.element(mcov, c(0, 1)))
  if (check == 0) warning("non-zero value(s) in mcov set to 1")
  mcov[mcov != 0] <- 1

  check <- sum(inits$theta[1, model$log.eta] <= 0)
  if (check) stop("illegal initial theta's")
  check <- sum(inits$omega <= 0)
  if (check) stop("illegal initial omega")
  # check = inits$sigma2<=0
  # if (check) stop("illegal initial sigma2")
  check <- sum(diag(covstruct) == 1)
  if (!check) stop("0 ETA's")
  y <- data$data[, "DV"]
  id <- data$data[, "ID"]
  check <- any(diff(unique(id)) != 1)
  if (check) stop("saem classic UI needs sequential ID. check your data")
  ntotal <- length(id)
  N <- length(unique(id))
  covariables <- if (is.null(model$covars)) NULL else unlist(stats::aggregate(.as.data.frame(data$data[, model$covars, drop = FALSE]), list(id), unique)[, -1, drop = FALSE])
  if (!is.null(covariables)) dim(covariables) <- c(N, data$N.covar)
  nb_measures <- table(id)
  ncov <- data$N.covar + 1
  nmc <- mcmc$nmc
  nM <- mcmc$nmc * N
  yM <- rep(y, nmc)
  mlen <- max(nb_measures)
  io <- t(sapply(nb_measures, function(x) rep(1:0, c(x, mlen - x))))
  ix <- rep(1:dim(io)[1], nmc)
  ioM <- io[ix, ]
  indioM <- grep(1, t(ioM)) - 1
  mPars <- if (ninputpars == 0) NULL else unlist(stats::aggregate(.as.data.frame(data$data[, inPars]), list(id), unique)[, -1])
  if (!is.null(mPars)) {
    dim(mPars) <- c(N, ninputpars)
    opt$mPars <- mPars
    ix <- rep(1:dim(mPars)[1], nmc)
    optM$mPars <- mPars[ix, ]
    dim(optM$mPars) <- c(nmc * N, ninputpars)
  }

  if (is.null(data$nmdat$CMT)) data$nmdat$CMT <- 1 ## CHECKME
  if (any(is.na(data$nmdat$CMT))) {
    stop("'CMT' has NA(s)")
  }
  ## CHECKME
  form <- attr(model$saem_mod, "form")
  .nobs <- 0
  if (form != "cls") {
    dat <- RxODE::etTrans(data$nmdat, attr(model$saem_mod, "rx"), TRUE, TRUE)
    .nobs <- attr(class(dat), ".RxODE.lst")$nobs
    ## if(length(dat) !=7) stop("SAEM doesn't support time varying covariates yet.");
    .rx <- attr(model$saem_mod, "rx")
    .pars <- .rx$params
    .pars <- setNames(rep(1.1, length(.pars)), .pars)
    ## opt$.rx <- .rx;
    opt$.pars <- .pars
    ## opt$.dat <- dat;
    dat <- .as.data.frame(dat[, -6])
    names(dat) <- toupper(names(dat))
    dat$ID <- as.integer(dat$ID)
  } else {
    dat <- data$nmdat[, c("ID", "TIME", "EVID", "AMT", "CMT")]
    infusion <- max(dat$EVID) > 10000
    if (infusion) {
      dat <- .fmtInfusionData(dat)
    } else {
      dat$DUR <- -1
    }
  }

  evt <- dat
  evt$ID <- evt$ID - 1
  evtM <- evt[rep(1:dim(evt)[1], nmc), ]
  evtM$ID <- cumsum(c(FALSE, diff(evtM$ID) != 0))

  i1 <- grep(1, diag(covstruct))
  i0 <- grep(0, diag(covstruct))
  nphi1 <- sum(diag(covstruct))
  nphi0 <- nphi - nphi1
  na <- length(mcmc$stepsize)
  nlambda1 <- sum(mcov[, i1])
  nlambda0 <- sum(mcov[, i0])
  nlambda <- nlambda1 + nlambda0
  nd1 <- nphi1 + nlambda1 + 1
  nd2 <- nphi1 + nlambda1 + nlambda0
  nb_param <- nd2 + 1
  Mcovariables <- cbind(rep(1, N), covariables)[, 1:nrow(mcov)]
  dim(Mcovariables) <- c(length(Mcovariables) / nrow(mcov), nrow(mcov)) # FIXME

  # get fixed ix
  fixed <- inits$theta.fix != ""
  wh <- fixed[, i1][mcov[, i1] == 1]
  len <- length(wh)
  fixed.i1 <- (1:len)[wh] - 1
  wh <- fixed[, i0][mcov[, i0] == 1]
  len <- length(wh)
  fixed.i0 <- (1:len)[wh] - 1

  jlog1 <- grep(T, model$log.eta)
  jcov <- grep(T, apply(mcov, 1, sum) > 0)
  covstruct1 <- covstruct[i1, i1]
  dim(covstruct1) <- c(nphi1, nphi1)
  ind_cov <- grep(1, mcov[mcov > 0])

  mcov1 <- matrix(mcov[, i1], ncol = length(i1))
  mcov0 <- matrix(mcov[, i0], nrow = nrow(mcov), ncol = length(i0))
  ind_cov1 <- grep(1, mcov1[mcov1 > 0]) - 1
  ind_cov0 <- grep(1, mcov0[mcov0 > 0]) - 1

  pc <- apply(mcov, 2, sum)
  ipc <- cumsum(c(0, pc[1:(nphi - 1)])) + 1
  ipcl1 <- ipc[jlog1]
  for (x in jlog1) inits$theta[1, x] <- log(inits$theta[1, x])

  idx <- as.vector(mcov1 > 0)
  COV1 <- Mcovariables[, row(mcov1)[idx]]
  dim(COV1) <- c(N, sum(idx))
  COV21 <- crossprod(COV1)

  x <- mcov1 * col(mcov1)
  x <- sapply(x[idx], function(x) {
    ret <- rep(0, nphi1)
    ret[x] <- 1
    ret
  })
  LCOV1 <- t(x)
  dim(LCOV1) <- c(nlambda1, nphi1)
  pc1 <- apply(LCOV1, 2, sum)

  x1 <- diag(sum(idx))
  diag(x1) <- inits$theta[, i1][idx]
  MCOV1 <- x1 %*% LCOV1
  jcov1 <- grep(1, LCOV1) - 1

  idx <- as.vector(mcov0 > 0)
  COV0 <- Mcovariables[, row(mcov0)[idx]]
  dim(COV0) <- c(N, sum(idx))
  COV20 <- crossprod(COV0)

  x <- mcov0 * col(mcov0)
  x <- sapply(x[idx], function(x) {
    ret <- rep(0, nphi0)
    ret[x] <- 1
    ret
  })
  LCOV0 <- t(x)
  dim(LCOV0) <- c(nlambda0, nphi0)

  x1 <- diag(sum(idx))
  diag(x1) <- inits$theta[, i0][idx]
  if (dim(x1)[1] > 0) {
    MCOV0 <- x1 %*% LCOV0
  } else {
    MCOV0 <- matrix(x1, nrow = 0, ncol = dim(LCOV0)[2])
  }
  jcov0 <- grep(1, LCOV0) - 1

  mprior_phi1 <- Mcovariables %*% inits$theta[, i1]
  mprior_phi0 <- Mcovariables %*% inits$theta[, i0]

  Gamma2_phi1 <- diag(nphi1)
  diag(Gamma2_phi1) <- inits$omega[i1]
  Gamma2_phi0 <- diag(nphi0)
  diag(Gamma2_phi0) <- inits$omega[i0]

  phiM <- matrix(0, N, nphi)
  phiM[, i1] <- mprior_phi1
  phiM[, i0] <- mprior_phi0
  phiM <- phiM[rep(1:N, nmc), , drop = FALSE]
  .tmp <- diag(sqrt(inits$omega))
  if (model$N.eta == 1) .tmp <- matrix(sqrt(inits$omega))
  phiM <- phiM + matrix(rnorm(phiM), dim(phiM)) %*% .tmp


  mc.idx <- rep(1:N, nmc)
  statphi <- sapply(1:nphi, function(x) tapply(phiM[, x], mc.idx, mean))
  statphi11 <- statphi[, i1]
  dim(statphi11) <- c(N, length(i1))
  statphi01 <- statphi[, i0]
  dim(statphi01) <- c(N, length(i0))
  statphi12 <- crossprod(phiM[, i1])
  statphi02 <- crossprod(phiM[, i0])

  # x = mcov *cumsum(mcov)
  # x1 = cbind(x[,i1], x[,i0])
  # indiphi = order(x1[x1>0])	#FINDME

  niter <- sum(mcmc$niter)
  niter_phi0 <- niter * .5
  nb_sa <- mcmc$niter[1] * .75
  nb_correl <- mcmc$niter[1] * .75
  va <- mcmc$stepsize
  vna <- mcmc$niter
  na <- length(va)
  pas <- 1 / (1:vna[1])^va[1]
  for (ia in 2:na)
  {
    end <- length(pas)
    k1 <- pas[end]^(-1 / va[ia])
    pas <- c(pas, 1 / ((k1 + 1):(k1 + vna[ia]))^va[ia])
  }
  pash <- c(rep(1, mcmc$burn.in), 1 / (1:niter))
  minv <- rep(1e-20, nphi)

  # preserve par order when printing iter history
  mcov[mcov == 1] <- 1:nlambda
  ilambda1 <- mcov[, i1]
  ilambda1 <- ilambda1[ilambda1 > 0] - 1
  ilambda0 <- mcov[, i0]
  ilambda0 <- ilambda0[ilambda0 > 0] - 1
  # print(mcov)
  # print(mcov[,i1]); print(ilambda1)
  # print(mcov[,i0]); print(ilambda0)

  i1 <- i1 - 1
  i0 <- i0 - 1
  opt$distribution <- distribution
  cfg <- list(
    inits = inits.save,
    nu = mcmc$nu,
    niter = niter,
    nb_sa = nb_sa,
    nb_correl = nb_correl,
    niter_phi0 = niter_phi0,
    nmc = nmc,
    coef_phi0 = .9638, # FIXME
    rmcmc = .5,
    coef_sa = .95,
    pas = pas,
    pash = pash,
    minv = minv,

    N = N,
    ntotal = ntotal,
    y = y,
    yM = yM,
    phiM = phiM,
    evt = as.matrix(evt),
    evtM = as.matrix(evtM),
    mlen = mlen,
    indioM = indioM,

    pc1 = pc1,
    covstruct1 = covstruct1,
    Mcovariables = Mcovariables,

    i1 = i1,
    i0 = i0,
    nphi1 = nphi1,
    nphi0 = nphi0,
    nlambda1 = nlambda1,
    nlambda0 = nlambda0,
    COV0 = COV0,
    COV1 = COV1,
    COV20 = COV20,
    COV21 = COV21,
    LCOV0 = LCOV0,
    LCOV1 = LCOV1,
    MCOV0 = MCOV0,
    MCOV1 = MCOV1,
    Gamma2_phi0 = Gamma2_phi0,
    Gamma2_phi1 = Gamma2_phi1,
    mprior_phi0 = mprior_phi0,
    mprior_phi1 = mprior_phi1,
    jcov0 = jcov0,
    jcov1 = jcov1,
    ind_cov0 = ind_cov0,
    ind_cov1 = ind_cov1,
    statphi11 = statphi11,
    statphi01 = statphi01,
    statphi02 = statphi02,
    statphi12 = statphi12,
    res.mod = model$res.mod,
    ares = inits$ares,
    bres = inits$bres,
    opt = opt,
    optM = optM,
    print = mcmc$print,
    distribution = distribution,
    par.hist = matrix(0, sum(niter), nlambda1 + nlambda0 + nphi1 + 1 + (model$res.mod > 2)),
    seed = seed,
    fixed.i1 = fixed.i1,
    fixed.i0 = fixed.i0,
    ilambda1 = as.integer(ilambda1),
    ilambda0 = as.integer(ilambda0),
    extraLL = .extraLL,
    nobs = .nobs
  )


  ## CHECKME
  s <- cfg$evt[cfg$evt[, "EVID"] == 0, "CMT"]
  cfg$opt$cmt_endpnt <- cfg$optM$cmt_endpnt <- sort(unique(s))
  cfg$nendpnt <- length(unique(s))
  if (model$nendpnt != cfg$nendpnt) {
    msg <- sprintf("mis-match in nbr endpoints in model & in data")
    stop(msg)
  }
  t <- unlist(split(1L:length(s), s))
  cfg$ysM <- rep(cfg$y[t], cfg$nmc)
  cfg$ix_sorting <- t - 1 # c-index for sorting by endpnt
  cfg$y_offset <- c(0, cumsum(table(s)))
  s <- cfg$evtM[cfg$evtM[, "EVID"] == 0, "CMT"]
  cfg$ix_endpnt <- as.integer(as.factor(s)) - 1 # to derive vecares & vecbres
  s <- cfg$evtM[cfg$evtM[, "EVID"] == 0, "ID"]
  t <- cumsum(c(0, table(s)))
  cfg$ix_idM <- cbind(t[-length(t)], t[-1] - 1) # c-index of obs records of each subject

  cfg$ares <- rep(10, cfg$nendpnt)
  cfg$bres <- rep(1, cfg$nendpnt)
  cfg$ares[cfg$res.mod == 2] <- 0
  cfg$bres[cfg$res.mod == 1] <- 0

  nres <- (1:2)[(cfg$res.mod == 3) + 1]
  cfg$res_offset <- cumsum(c(0, nres))
  cfg$par.hist <- matrix(0, cfg$niter, nlambda1 + nlambda0 + nphi1 + sum(nres))

  cfg$DEBUG <- cfg$opt$DEBUG <- cfg$optM$DEBUG <- DEBUG
  cfg$phiMFile <- tempfile()
  cfg
}



reINITS <- "^\\s*initCondition\\s*=\\s*c\\((?<inits>.+)\\)\\s*$"
reDATAPAR <- "^\\s*ParamFromData\\s*=\\s*c\\((?<inits>.+)\\)\\s*$"

getInits <- function(x, re, collapse = TRUE) {
  # x = "initCondition = c(1,2,3,4)"
  # re = "^\\s*initCondition\\s*=\\s*c\\((?<inits>.+)\\)\\s*$"
  m <- regexpr(re, x, perl = T)
  if (m < 0) stop("invalid initCondition input")
  start <- attr(m, "capture.start")
  len <- attr(m, "capture.length")
  inits <- substring(x, start, start + len - 1)
  s <- strsplit(inits, ",")[[1]]

  if (collapse) {
    paste0(sprintf("\t_inits[%d] = %s;", 1:length(s) - 1, s), collapse = "\n")
  } else {
    s
  }
}

#' Print an SAEM model fit summary
#'
#' Print an SAEM model fit summary
#'
#' @param object a saemFit object
#' @param ... others
#' @return a list
#' @export
summary.saemFit <- function(object, ...) {
  fit <- object ## Rcheck hack

  th <- fit$Plambda
  nth <- length(th)
  H <- solve(fit$Ha[1:nth, 1:nth])
  se <- sqrt(diag(H))

  m <- cbind(exp(th), th, se) # FIXME
  ## lhsVars = scan("LHS_VARS.txt", what="", quiet=T)
  ## if (length(lhsVars)==nth) dimnames(m)[[1]] = lhsVars
  dimnames(m)[[2]] <- c("th", "log(th)", "se(log_th)")
  cat("THETA:\n")
  print(m)
  cat("\nOMEGA:\n")
  print(fit$Gamma2_phi1)
  if (any(fit$sig2 == 0)) {
    cat("\nSIGMA:\n")
    print(max(fit$sig2^2))
  } else {
    cat("\nARES & BRES:\n")
    print(fit$sig2)
  }

  invisible(list(theta = th, se = se, H = H, omega = fit$Gamma2_phi1, eta = fit$mpost_phi))
}

#' Print an SAEM model fit summary
#'
#' Print an SAEM model fit summary
#'
#' @param x a saemFit object
#' @param ... others
#' @return a list
#' @export
print.saemFit <- function(x, ...) {
  fit <- x ## Rcheck hack

  th <- fit$Plambda
  nth <- length(th)
  H <- solve(fit$Ha[1:nth, 1:nth])
  se <- sqrt(diag(H))

  m <- cbind(exp(th), th, se) # FIXME
  ## lhsVars = scan("LHS_VARS.txt", what="", quiet=T)
  ## if (length(lhsVars)==nth) dimnames(m)[[1]] = lhsVars
  dimnames(m)[[2]] <- c("th", "log(th)", "se(log_th)")
  cat("THETA:\n")
  print(m)
  cat("\nOMEGA:\n")
  print(fit$Gamma2_phi1)
  if (any(fit$sig2 == 0)) {
    cat("\nSIGMA:\n")
    print(max(fit$sig2^2))
  } else {
    cat("\nARES & BRES:\n")
    print(fit$sig2)
  }

  invisible(list(theta = th, se = se, H = H, omega = fit$Gamma2_phi1, eta = fit$mpost_phi))
}


#' Fit an SAEM model
#'
#' Fit an SAEM model using either closed-form solutions or ODE-based model definitions
#'
#' @param model an RxODE model or lincmt()
#' @param data input data
#' @param inits initial values
#' @param PKpars PKpars function
#' @param pred  pred function
#' @param covars Covariates in data
#' @param mcmc a list of various mcmc options
#' @param ODEopt optional ODE solving options
#' @inheritParams configsaem
#' @param seed seed for random number generator
#' @details
#'    Fit a generalized nonlinear mixed-effect model using the Stochastic
#'    Approximation Expectation-Maximization (SAEM) algorithm
#'
#' @author Matthew Fidler & Wenping Wang
#' @export
saem.fit <- function(model, data, inits,
                     PKpars = NULL, pred = NULL,
                     covars = NULL,
                     mcmc = list(niter = c(200, 300), nmc = 3, nu = c(2, 2, 2)),
                     ODEopt = list(atol = 1e-06, rtol = 1e-04, method = "lsoda", transitAbs = FALSE),
                     distribution = c("normal", "poisson", "binomial", "lnorm"),
                     seed = 99) {
  UseMethod("saem.fit")
}
##' @rdname saem.fit
##' @export
saem <- saem.fit

##' @rdname saem.fit
##' @export
saem.fit.nlmixr.ui.nlme <- function(model, data, inits,
                                    PKpars = NULL, pred = NULL,
                                    covars = NULL,
                                    mcmc = list(niter = c(200, 300), nmc = 3, nu = c(2, 2, 2)),
                                    ODEopt = list(atol = 1e-06, rtol = 1e-04, method = "lsoda", transitAbs = FALSE),
                                    distribution = c("normal", "poisson", "binomial", "lnorm"),
                                    seed = 99) {
  call <- as.list(match.call(expand.dots = TRUE))[-1]
  names(call)[1] <- "object"
  call$est <- "saem"
  return(do.call(getFromNamespace("nlmixr", "nlmixr"), call, envir = parent.frame(1)))
}

##' @rdname saem.fit
##' @export
saem.fit.function <- saem.fit.nlmixr.ui.nlme

##' @rdname saem.fit
##' @export
saem.fit.nlmixrUI <- saem.fit.nlmixr.ui.nlme

##' @rdname saem.fit
##' @export
saem.fit.RxODE <- function(model, data, inits,
                           PKpars = NULL, pred = NULL,
                           covars = NULL,
                           mcmc = list(niter = c(200, 300), nmc = 3, nu = c(2, 2, 2)),
                           ODEopt = list(atol = 1e-06, rtol = 1e-04, method = "lsoda", transitAbs = FALSE),
                           distribution = c("normal", "poisson", "binomial", "lnorm"),
                           seed = 99) {
  saem_fit <- gen_saem_user_fn(model, PKpars, pred)
  model <- list(saem_mod = saem_fit, covars = covars)
  cfg <- configsaem(model, data, inits, mcmc, ODEopt, distribution, seed)
  fit <- saem_fit(cfg)
  ## dyn.unload("saem_main.dll")
  fit
}

##' @rdname saem.fit
##' @export
saem.fit.default <- function(model, data, inits,
                             PKpars = NULL, pred = NULL,
                             covars = NULL,
                             mcmc = list(niter = c(200, 300), nmc = 3, nu = c(2, 2, 2)),
                             ODEopt = list(atol = 1e-06, rtol = 1e-04, method = "lsoda", transitAbs = FALSE),
                             distribution = c("normal", "poisson", "binomial", "lnorm"),
                             seed = 99) {
  saem_fit <- gen_saem_user_fn(model)
  model <- list(saem_mod = saem_fit, covars = covars)
  cfg <- configsaem(model, data, inits, mcmc, ODEopt, distribution, seed)
  fit <- saem_fit(cfg)
  # dyn.unload("saem_main.dll")
  fit
}

##' @export
ranef.saemFit <- function(object, ...) {
  return(object$eta)
}

##' @export
fixef.saemFit <- function(object, ...) {
  return(object$Plambda)
}

focei.theta.saemFit <- function(object, uif, ...) {
  ## Get the thetas needed for FOCEi fit.
  this.env <- environment()
  if (class(uif) == "function") {
    uif <- nlmixr(uif)
  }
  n <- uif$focei.names
  thetas <- structure(rep(NA, length(n)), .Names = n)
  sf <- structure(as.vector(fixed.effects(object)), .Names = uif$saem.theta.name)
  for (n in names(sf)) {
    thetas[n] <- sf[n]
  }
  .predDf <- uif$predDf
  .ini <- .as.data.frame(uif$ini)
  .resMat <- object$resMat
  sapply(seq_along(.predDf$cond), function(i) {
    x <- paste(.predDf$cond[i])
    .tmp <- .ini[which(.ini$condition == x), ]
    .w <- which(sapply(.tmp$err, function(x) any(x == "prop")))
    if (length(.w) == 1) {
      thetas[paste(.tmp$name[.w])] <- .resMat[i, 2]
    }
    .w <- which(sapply(.tmp$err, function(x) {
      any(x == c(
        "add", "norm", "dnorm", "lnorm", "dlnorm",
        "dlogn", "logn"
      ))
    }))
    if (length(.w) == 1) {
      thetas[paste(.tmp$name[.w])] <- .resMat[i, 1]
    }
    assign("thetas", thetas, this.env)
  })
  return(thetas)
}

focei.eta.saemFit <- function(object, uif, ...) {
  if (class(uif) == "function") {
    uif <- nlmixr(uif)
  }
  ## Reorder based on translation
  eta.trans <- uif$saem.eta.trans
  for (i in seq(1, max(eta.trans))) {
    while (!(any(i == eta.trans)) && max(eta.trans) > i) {
      eta.trans[eta.trans >= i] <- eta.trans[eta.trans >= i] - 1
    }
  }
  ## orig eta ->  new eta
  df <- .as.data.frame(uif$ini)
  eta <- df[!is.na(df$neta1), ]
  len <- length(eta$name)
  cur.lhs <- character()
  cur.rhs <- numeric()
  ome <- character()
  cur.ome <- object$Gamma2_phi1
  for (i in seq_along(eta$name)) {
    last.block <- FALSE
    if (i == len) {
      last.block <- TRUE
    } else if (eta$neta1[i + 1] == eta$neta2[i + 1]) {
      last.block <- TRUE
    }
    if (eta$neta1[i] == eta$neta2[i]) {
      cur.lhs <- c(cur.lhs, sprintf("ETA[%d]", eta$neta1[i]))
      cur.rhs <- c(cur.rhs, cur.ome[eta.trans[eta$neta1[i]], eta.trans[eta$neta2[i]]])
      if (last.block) {
        ome[length(ome) + 1] <- sprintf(
          "%s ~ %s", paste(cur.lhs, collapse = " + "),
          paste(deparse(cur.rhs), collapse = " ")
        )
        cur.lhs <- character()
        cur.rhs <- numeric()
      }
    } else {
      cur.rhs <- c(cur.rhs, cur.ome[eta.trans[eta$neta1[i]], eta.trans[eta$neta2[i]]])
    }
  }
  ome <- eval(parse(text = sprintf("list(%s)", paste(ome, collapse = ","))))
  return(ome)
}

as.focei.saemFit <- function(object, uif, pt = proc.time(), ..., data, calcResid = TRUE, obf = NULL,
                             nnodes.gq = 1, nsd.gq = 3, adjObf = TRUE,
                             calcCov = TRUE, covMethod = NULL,
                             calcCovTime = NULL, calcTables = TRUE) {
  .saemCfg <- attr(object, "saem.cfg")
  .saemTime <- proc.time() - pt
  if (class(uif) == "function") {
    uif <- nlmixr(uif)
  }
  .dist <- ""
  if (any(uif$saem.distribution == c("poisson", "binomial"))) {
    calcResid <- NA
    .dist <- uif$saem.distribution
  }
  uif.new <- uif
  fit <- object
  mat <- random.effects(fit)
  ## Reorder based on translation
  eta.trans <- uif$saem.eta.trans
  for (i in seq(1, max(eta.trans))) {
    while (!(any(i == eta.trans)) && max(eta.trans) > i) {
      eta.trans[eta.trans >= i] <- eta.trans[eta.trans >= i] - 1
    }
  }
  mat2 <- mat[, eta.trans, drop = FALSE]
  th <- focei.theta(fit, uif)
  for (n in names(th)) {
    uif.new$est[uif.new$name == n] <- th[n]
  }
  ome <- focei.eta(fit, uif)
  init <- list(
    THTA = as.vector(th),
    OMGA = ome
  )
  saem.time <- proc.time() - pt
  if (missing(data)) {
    stop("Requires Data...")
  } else {
    dat <- data
  }
  .tn <- uif$saem.theta.name
  .ini <- .as.data.frame(uif$ini)
  .ini <- .ini[uif$ini$name %in% .tn, ]
  if (any(.ini$fix)) {
    .fixed <- paste(.ini$name[.ini$fix])
    .tn <- .tn[!(.tn %in% .fixed)]
  }
  .calcCov <- calcCov
  if (uif$env$covMethod == "") {
    .cov <- NULL
    .addCov <- FALSE
  } else if (inherits(calcCov, "matrix")) {
    .cov <- calcCov
    .addCov <- TRUE
  } else {
    .nth <- length(.tn)

    .ini <- .as.data.frame(uif$ini)
    .ini <- .ini[is.na(.ini$err), ]
    .ini <- .ini[!is.na(.ini$ntheta), ]
    .ini <- .ini[!.ini$fix, ]
    .ini <- paste(.ini$name)
    .calcCovTime <- proc.time()
    if (calcCov) {
      .covm <- object$Ha[1:.nth, 1:.nth]
      .covm <- try(calc.COV(object))
      .doIt <- !inherits(.covm, "try-error")
      if (.doIt && dim(.covm)[1] != .nth) .doIt <- FALSE
      if (.doIt) {
        .tmp <- try(chol(.covm), silent = TRUE)
        .addCov <- TRUE
        .sqrtm <- FALSE
        if (inherits(.tmp, "try-error")) {
          .tmp <- .covm
          .tmp <- try(sqrtm(.tmp %*% t(.tmp)), silent = FALSE)
          if (inherits(.tmp, "try-error")) {
            .calcCov <- FALSE
            .covm <- object$Ha[1:.nth, 1:.nth]
            .tmp <- try(chol(.covm), silent = TRUE)
            .addCov <- TRUE
            .sqrtm <- FALSE
            if (inherits(.tmp, "try-error")) {
              .tmp <- object$Ha[1:.nth, 1:.nth]
              .tmp <- try(sqrtm(.tmp %*% t(.tmp)), silent = FALSE)
              if (inherits(.tmp, "try-error")) {
                .addCov <- FALSE
              } else {
                .sqrtm <- TRUE
              }
            } else {
              .tmp <- object$Ha[1:.nth, 1:.nth]
            }
          } else {
            .sqrtm <- TRUE
          }
        } else {
          .tmp <- .covm
        }
      } else {
        .tmp <- object$Ha[1:.nth, 1:.nth]
        .tmp <- try(chol(.covm), silent = TRUE)
        .calcCov <- FALSE
        .addCov <- TRUE
        .sqrtm <- FALSE
        if (inherits(.tmp, "try-error")) {
          .tmp <- object$Ha[1:.nth, 1:.nth]
          .tmp <- try(sqrtm(.tmp %*% t(.tmp)), silent = FALSE)
          if (inherits(.tmp, "try-error")) {
            .addCov <- FALSE
          } else {
            .sqrtm <- TRUE
          }
        } else {
          .tmp <- object$Ha[1:.nth, 1:.nth]
          .calcCov <- FALSE
        }
      }
    } else {
      .tmp <- try(chol(.covm), silent = TRUE)
      .addCov <- TRUE
      .sqrtm <- FALSE
      if (inherits(.tmp, "try-error")) {
        .tmp <- object$Ha[1:.nth, 1:.nth]
        .tmp <- try(sqrtm(.tmp %*% t(.tmp)), silent = FALSE)
        if (inherits(.tmp, "try-error")) {
          .addCov <- FALSE
        } else {
          .sqrtm <- TRUE
        }
      } else {
        .tmp <- object$Ha[1:.nth, 1:.nth]
        .calcCov <- FALSE
      }
    }
    if (.addCov) {
      if (!.calcCov) {
        .cov <- RxODE::rxInv(.tmp)
      } else {
        .cov <- .tmp
      }
      attr(.cov, "dimnames") <- list(.tn, .tn)
      .cov <- .cov[.ini, .ini, drop = FALSE]
    }
    .calcCovTime <- proc.time() - .calcCovTime
    .calcCovTime <- .calcCovTime["elapsed"]
  }
  .ini <- .as.data.frame(uif$ini)
  .ini <- .ini[!is.na(.ini$ntheta), ]
  .skipCov <- !is.na(.ini$err)
  .fixed <- uif$focei.fixed
  .skipCov <- .skipCov | .fixed[seq_along(.skipCov)]
  .covMethod <- uif$env$covMethod
  if (!any(.covMethod == c("", "r", "s", "r,s"))) {
    .covMethod <- ""
  }
  if (is.na(calcResid)) .covMethod <- ""
  .allThetaNames <- c(uif$saem.theta.name, uif$saem.omega.name, uif$saem.res.name)
  .m <- object$par_hist
  if (ncol(.m) > length(.allThetaNames)) {
    .m <- .m[, seq_along(.allThetaNames)]
  }
  if (.dist == "binomial") {
    .dist <- "bernoulli"
  }
  dimnames(.m) <- list(NULL, .allThetaNames)
  .fixedNames <- paste(uif$ini$name[which(uif$ini$fix)])
  .rn <- ""
  .likTime <- 0
  if (is.na(obf)) {
    .saemObf <- NA
  } else if (is.null(obf)) {
    .likTime <- proc.time()
    .saemObf <- calc.2LL(object, nnodes.gq = nnodes.gq, nsd.gq = nsd.gq)
    .likTime <- proc.time() - .likTime
    .likTime <- .likTime["elapsed"]
    if (nnodes.gq == 1) {
      .rn <- paste0("laplace", nsd.gq)
    } else {
      .rn <- paste0("gauss", nnodes.gq, "_", nsd.gq)
    }
  } else if (is(obf, "logical")) {
    if (is.na(obf)) {
      .saemObf <- NA
    } else if (obf) {
      .likTime <- proc.time()
      .saemObf <- calc.2LL(object, nnodes.gq = nnodes.gq, nsd.gq = nsd.gq)
      .likTime <- proc.time() - .likTime
      .likTime <- .likTime["elapsed"]
      if (nnodes.gq == 1) {
        .rn <- paste0("laplace", nsd.gq)
      } else {
        .rn <- paste0("gauss", nnodes.gq, "_", nsd.gq)
      }
    } else {
      .saemObf <- NA
    }
  } else if (is(object, "numeric")) {
    .saemObf <- obf
  }
  .notCalced <- TRUE
  .cwresTime <- proc.time()
  while (.notCalced) {
    .env <- new.env(parent = emptyenv())
    .env$nobs2 <- .saemCfg$nobs
    .env$nnodes.gq <- nnodes.gq
    .env$nsd.gq <- nsd.gq
    .env$adjObf <- adjObf
    .env$method <- "SAEM"
    .env$uif <- uif
    .env$saem <- object
    if (.addCov) {
      .env$cov <- .cov
    }
    .env$parHistStacked <- .data.frame(
      val = as.vector(.m),
      par = rep(.allThetaNames, each = nrow(.m)),
      iter = rep(1:nrow(.m), ncol(.m))
    )
    .env$parHist <- .data.frame(iter = rep(1:nrow(.m)), .as.data.frame(.m))
    if (length(.fixedNames) > 0) {
      .env$parHistStacked <- .env$parHistStacked[!(.env$parHistStacked$par %in% .fixedNames), , drop = FALSE]
      .env$parHist <- .env$parHist[, !(names(.env$parHist) %in% .fixedNames), drop = FALSE]
    }
    if (is.na(calcResid)) {
      .setSaemExtra(.env, .rn)
      .env$theta <- .data.frame(
        lower = -Inf, theta = init$THTA, upper = Inf, fixed = .fixed[seq_along(init$THTA)],
        row.names = uif$focei.names
      )
      .env$fullTheta <- setNames(init$THTA, uif$focei.names)
      .om0 <- .genOM(.parseOM(init$OMGA))
      attr(.om0, "dimnames") <- list(uif$eta.names, uif$eta.names)
      .env$omega <- .om0
      .env$etaObf <- .data.frame(ID = seq_along(mat2[, 1]), setNames(.as.data.frame(mat2), uif$eta.names), OBJI = NA)
      .env$noLik <- FALSE
      .env$objective <- .saemObf
    } else if (calcResid) {
      .setSaemExtra(.env, .rn)
    } else {
      .setSaemExtra(.env, .rn)
      .env$theta <- .data.frame(
        lower = -Inf, theta = init$THTA, upper = Inf, fixed = .fixed[seq_along(init$THTA)],
        row.names = uif$focei.names
      )
      .env$fullTheta <- setNames(init$THTA, uif$focei.names)
      .om0 <- .genOM(.parseOM(init$OMGA))
      attr(.om0, "dimnames") <- list(uif$eta.names, uif$eta.names)
      .env$omega <- .om0
      .env$etaObf <- .data.frame(ID = seq_along(mat2[, 1]), setNames(.as.data.frame(mat2), uif$eta.names), OBJI = NA)
      .env$noLik <- TRUE
      .env$objective <- .saemObf
    }
    .ctl <- uif$env$ODEopt
    names(.ctl) <- sub("maxsteps", "maxstepsOde", names(.ctl))
    .ctl <- .ctl[names(.ctl) != "scale"]
    .ctl$maxOuterIterations <- 0
    .ctl$maxInnerIterations <- 0
    .ctl$covMethod <- .covMethod
    .ctl$sumProd <- uif$env$sum.prod
    .ctl$optExpression <- uif$env$optExpression
    .ctl$scaleTo <- 0
    .ctl$calcTables <- calcTables
    .ctl <- do.call(foceiControl, .ctl)
    fit.f <- try(foceiFit.data.frame(
      data = dat,
      inits = init,
      PKpars = uif$theta.pars,
      ## par_trans=fun,
      model = uif$rxode.pred,
      pred = function() {
        return(nlmixr_pred)
      },
      err = uif$error,
      lower = uif$focei.lower,
      upper = uif$focei.upper,
      thetaNames = uif$focei.names,
      etaNames = uif$eta.names,
      etaMat = mat2,
      env = .env,
      fixed = .fixed,
      skipCov = .skipCov,
      control = .ctl
    ), silent = FALSE)
    if (inherits(fit.f, "try-error")) {
      if (is.na(calcResid)) {
        warning("Error calculating nlmixr object, return classic object")
        .notCalced <- FALSE
        return(object)
      } else if (calcResid) {
        calcResid <- FALSE
      } else {
        calcResid <- NA
      }
    } else {
      .notCalced <- FALSE
    }
  }
  .cwresTime <- proc.time() - .cwresTime
  if (is.na(calcResid)) {
    .cwresTime <- 0
  } else if (!calcResid) {
    .cwresTime <- 0
  }
  .env <- fit.f$env
  if (uif$env$covMethod == "") {
  } else if (inherits(.calcCov, "matrix")) {
    if (!is.null(covMethod)) {
      .env$covMethod <- covMethod
    }
    .calcCovTime <- calcCovTime
  } else if (.calcCov) {
    .env$covMethod <- "linFim"
    if (.addCov & .sqrtm) {
      .env$covMethod <- "|linFim|"
      warning("Covariance matrix non-positive definite, corrected by sqrtm(linFim %*% linFim)")
    }
  } else {
    if (calcCov) {
      warning("Linearization of FIM could not be used to calculate covariance.")
    }
    if (.addCov & .sqrtm) {
      .env$covMethod <- "|fim|"
      warning("Covariance matrix non-positive definite, corrected by sqrtm(fim %*% fim)")
    } else if (!.addCov) {
      warning("FIM non-positive definite and cannot be used to calculate the covariance")
    }
  }

  if (is.null(.env$time)) {
    .env$time <- .data.frame(saem = .saemTime["elapsed"], check.names = FALSE, row.names = c(""))
  } else {
    .time <- .env$time
    .time <- .time[, !(names(.time) %in% c("optimize", "covariance"))]
    .saemTime <- .saemTime["elapsed"]
    if (calcResid) {
      .saemTime <- .saemTime - .cwresTime["elapsed"]
      .time <- .data.frame(.time, cwres = .cwresTime["elapsed"], check.names = FALSE)
    }
    if (.likTime > 0) {
      .time <- .data.frame(.time, logLik = .likTime, check.names = FALSE)
      .saemTime <- .saemTime - .likTime
    }
    if (uif$env$covMethod != "") {
      .saemTime <- .saemTime - .calcCovTime
      .time <- .data.frame(.time, covariance = .calcCovTime, check.names = FALSE)
    }
    .env$time <- .data.frame(saem = .saemTime, .time, check.names = FALSE, row.names = c(""))
  }
  .env$message <- ""
  if (is.na(calcResid)) {
    row.names(.env$objDf) <- .rn
  } else if (calcResid) {
    if (!is.na(.saemObf)) {
      .llik <- -.saemObf / 2
      .nobs <- .env$nobs
      attr(.llik, "df") <- attr(get("logLik", .env), "df")
      .objf <- ifelse(.env$adjObf, .saemObf - .nobs * log(2 * pi), .saemObf)
      ## for (.t in c("OBJF","objective", "objf")){
      ##   assign(.t,.objf,.env);
      ## }
      .tmp <- .data.frame(
        OBJF = .objf,
        AIC = .saemObf + 2 * attr(get("logLik", .env), "df"),
        BIC = .saemObf + log(.env$nobs) * attr(get("logLik", .env), "df"),
        "Log-likelihood" = as.numeric(.llik), check.names = FALSE
      )
      if (any(names(.env$objDf) == "Condition Number")) {
        .tmp <- .data.frame(.tmp,
          "Condition Number" = .env$objDf[, "Condition Number"],
          check.names = FALSE
        )
      }
      .env$objDf <- rbind(.env$objDf, .tmp)
      row.names(.env$objDf) <- c("FOCEi", .rn)
    }
    .setSaemExtra(.env, "FOCEi")
  } else {
    row.names(.env$objDf) <- .rn
  }
  if (inherits(fit.f, "nlmixrFitData")) {
    .cls <- class(fit.f)
    .env <- attr(.cls, ".foceiEnv")
    .cls[1] <- "nlmixrSaem"
    class(fit.f) <- .cls
  }
  return(fit.f)
}

## FIXME: coef_phi0, rmcmc, coef_sa
## FIXME: Klog, rho, sa, nmc
## FIXME: N.design
## FIXME: g = gc = 1
## FIXME: ODE inits
## FIXME: Tinf for ODE
## FIXME: chk infusion poor fit

## Local Variables:
## ess-indent-level: 2
## indent-tabs-mode: nil
## End:
