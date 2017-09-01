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
## along with nlmixr.  If not, see <http:##www.gnu.org/licenses/>.

saem_ode_str = '#include <RcppArmadillo.h>
#include "saem_class_rcpp.hpp"


using namespace std;
using namespace arma;

extern "C" {

void <%=ode_solver%>(
	int *neq,
	double *theta,	//order:
	double *time,
	int *evid,
	int *ntime,
	double *inits,
	double *dose,
	double *ret,
	double *atol,
	double *rtol,
	int *stiff,
	int *transit_abs,
	int *nlhs,
	double *lhs,
    int *rc
);

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


vec user_function(const mat &phi, const mat &evt, const List &opt) {
  uvec ix;
  vec id = evt.col(0);
  mat wm;
  vec wv;

  ix = find(evt.col(2) == 0);
  vec yp(ix.n_elem);
  double *p=yp.memptr();
  int N=id.max()+1;

  for (int i=0; i<N; i++) {
    wm = evt.rows( find(id == i) );

    vec time__;
    time__ = wm.col(1);
    int ntime = time__.n_elem;
    wv = wm.col(2);
    ivec evid(ntime);
    for (int k=0; k<ntime; ++k) evid(k) = wv(k);
    vec amt;
    amt = wm.col(3);
    amt = amt( find(evid > 0) );

    int neq=as<int>(opt["neq"]);
    vec inits(neq);
    inits.zeros(); //as<vec>(opt["inits"]);	//FIXME

<%=foo%>

<%=pars%>
<%=inits%>

    int stiff=as<int>(opt["stiff"]);
    int transit_abs=as<int>(opt["transit_abs"]);
    int nlhs=as<int>(opt["nlhs"]);
    double atol=as<double>(opt["atol"]);
    double rtol=as<double>(opt["rtol"]);
    int rc=0;

    mat ret(neq, ntime);
    mat lhs(nlhs, ntime);
	<%=ode_solver%>(&neq, params.memptr(), time__.memptr(),
	    evid.memptr(), &ntime, inits.memptr(), amt.memptr(), ret.memptr(),
	    &atol, &rtol, &stiff, &transit_abs, &nlhs, lhs.memptr(), &rc);
	ret = join_cols(join_cols(time__.t(), ret), lhs).t();
	ret = ret.rows( find(evid == 0) );

<%=model_vars_decl%>

vec g;
g = <%=pred_expr%>;

    int no = g.n_elem;
    memcpy(p, g.memptr(), no*sizeof(double));
    p += no;
  }

  return yp;
}

// definition
RcppExport SEXP dopred( SEXP in_phi, SEXP in_evt, SEXP in_opt ) {
BEGIN_RCPP
    Rcpp::traits::input_parameter< mat& >::type phi(in_phi);
    Rcpp::traits::input_parameter< mat& >::type evt(in_evt);
    List opt(in_opt);

    vec g = user_function(phi, evt, opt);
    return Rcpp::wrap(g);
END_RCPP
}

RcppExport SEXP saem_fit(SEXP xSEXP) {
BEGIN_RCPP
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
END_RCPP
}
'


saem_cmt_str = '#include <RcppArmadillo.h>
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


vec user_function(const mat &phi, const mat &evt, const List &opt) {
  uvec ix;
  vec id = evt.col(0);
  mat wm;
  vec obs_time, dose_time, dose, wv;

  ix = find(evt.col(2) == 0);
  vec yp(ix.n_elem);
  double *p=yp.memptr();
  int N=id.max()+1;

  for (int i=0; i<N; i++) {
    ix = find(id == i);
    wm = evt.rows(ix);

    ix = find(wm.col(2) == 0);
    wv = wm.col(1);
    wv = wv(ix);
    const Map<MatrixXd> _obs_time(wv.memptr(), wv.n_elem, 1);
    const VectorXd obs_time(_obs_time);

    ix = find(wm.col(2) > 0);
    wv = wm.col(1);
    wv = wv(ix);
    const Map<MatrixXd> _dose_time(wv.memptr(), wv.n_elem, 1);
    const VectorXd dose_time(_dose_time);

    wv = wm.col(3);
    wv = wv(ix);
    const Map<MatrixXd> _dose(wv.memptr(), wv.n_elem, 1);
    const VectorXd dose(_dose);

    wv = wm.col(4);
    wv = wv(ix);
    const Map<MatrixXd> _Tinf(wv.memptr(), wv.n_elem, 1);
    const VectorXd Tinf(_Tinf);

<%=foo%>
<%=pars%>

    int no=obs_time.size();
    VectorXd g(obs_time.size());
    int ncmt=<%=ncmt%>, oral=<%=oral%>, infusion=<%=infusion%>, parameterization=<%=parameterization%>;

	g = generic_cmt_interface(
      obs_time,
      dose_time,
      dose,
      Tinf,
      params,
      ncmt,
      oral,
      infusion,
      parameterization);

    memcpy(p, g.data(), no*sizeof(double));
    p += no;
	//cout << "ok " << i <<endl;
  }

  return yp;
}

// definition
RcppExport SEXP dopred( SEXP in_phi, SEXP in_evt, SEXP in_opt ) {
BEGIN_RCPP
    Rcpp::traits::input_parameter< mat& >::type phi(in_phi);
    Rcpp::traits::input_parameter< mat& >::type evt(in_evt);
    List opt(in_opt);
    vec g = user_function(phi, evt, opt);
    return Rcpp::wrap(g);
END_RCPP
}

RcppExport SEXP saem_fit(SEXP xSEXP) {
BEGIN_RCPP
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
END_RCPP
}
'


##' Get a list of directories for inclusion
##'
##' @param pkg a string or list of string for packages to be included.
##' @return An inclusion string for Makevars
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
nmxInclude <- function(pkg="nlmixr"){
    if  (length(pkg) == 1){
        x <- system.file("", package = pkg)
        if(.Platform$OS.type=="windows") x <- gsub("\\\\", "/", utils::shortPathName(x))
        x <- paste0("-I", x, "/include");
        return(x);
    } else {
        paste(sapply(pkg, nmxInclude), collapse=" ");
    }
}

#' Generate an SAEM model
#'
#' Generate an SAEM model using either closed-form solutions or ODE-based model definitions
#'
#' @param model a compiled SAEM model by gen_saem_user_fn()
#' @param PKpars PKpars function
#' @param pred  pred function
#' @details
#'    Fit a generalized nonlinear mixed-effect model using the Stochastic
#'    Approximation Expectation-Maximization (SAEM) algorithm
#'
#' @author Wenping Wang
#' @examples
#' \dontrun{
#' #ode <- "d/dt(depot) =-KA*depot;
#' #        d/dt(centr) = KA*depot - KE*centr;"
#' #m1 = RxODE(ode, modName="m1")
#' #ode <- "C2 = centr/V;
#' #      d/dt(depot) =-KA*depot;
#' #      d/dt(centr) = KA*depot - KE*centr;"
#' #m2 = RxODE(ode, modName="m2")
#'
#' PKpars = function()
#' {
#'   CL = exp(lCL)
#'   V  = exp(lV)
#'   KA = exp(lKA)
#'   KE = CL / V
#'   #initCondition = c(0, KE - CL/V)
#' }
#' PRED = function() centr / V
#' PRED2 = function() C2
#'
#' saem_fit <- gen_saem_user_fn(model=lincmt(ncmt=1, oral=T))
#' #saem_fit <- gen_saem_user_fn(model=m1, PKpars, pred=PRED)
#' #saem_fit <- gen_saem_user_fn(model=m2, PKpars, pred=PRED2)
#'
#'
#' #--- saem cfg
#' nmdat = read.table("theo_sd.dat",  head=T)
#' model = list(saem_mod=saem_fit, covars="WT")
#' inits = list(theta=c(.05, .5, 2))
#' cfg   = configsaem(model, nmdat, inits)
#' cfg$print = 50
#'
#' #cfg$Rfn = nlmixr:::Ruser_function_cmt
#' #dyn.load("m1.d/m1.so");cfg$Rfn = nlmixr:::Ruser_function_ode
#' fit = saem_fit(cfg)
#' df = simple.gof(fit)
#' xyplot(DV~TIME|ID, df, type=c("p","l"), lwd=c(NA,1), pch=c(1,NA), groups=grp)
#' fit
#' }
#' @export
gen_saem_user_fn = function(model, PKpars=attr(model, "default.pars"), pred=NULL) {
  is.ode = class(model) == "RxODE"
  is.win <- .Platform$OS.type=="windows"
  env = environment()
  saem.cpp <- paste0(basename(tempfile(pattern="saem", getwd())), .Platform$r_arch);
  saem.base <- saem.cpp
  saem.dll <- paste0(saem.cpp, .Platform$dynlib.ext)
  saem.cpp <- paste0(saem.cpp, ".cpp");
  lwd <- getwd();
  .wd <- tempfile()
  dir.create(.wd, recursive = TRUE)
  setwd(.wd)
  on.exit({setwd(lwd);unlink(.wd, recursive=TRUE, force=TRUE)});

  if (is.ode) {
    modelVars = model$cmpMgr$get.modelVars()

    pars = modelVars$params
    npar = length(pars)
    pars = paste(c(
    sprintf("    vec params(%d);\n", npar),
    sprintf("    params(%d) = %s;\n", 1:npar-1, pars)
    ), collapse="")

	model_vars = names=c("time", modelVars$state, modelVars$lhs)
	s = lapply(1:length(model_vars), function(k) {
		sprintf("vec %s;\n%s=ret.col(%d);\n", model_vars[k], model_vars[k], k-1)
	})
	model_vars_decl = paste0(s, collapse="")

	ode_solver = model$cmpMgr$ode_solver

	neq = length(modelVars$state)
	nlhs = length(modelVars$lhs)

    x = deparse(body(pred))
    len = length(x)
    pred_expr = if(len>2) x[2:(len-1)] else x
  } else {
	neq = nlhs = 0
	pars = model
	list2env(attr(pars, "calls"), envir=env)
  }
  #str(ls(envir=env))

  ### deal with explicit initCondition statement
  x = deparse(body(PKpars))
  ix = grep(reINITS, x, perl=T)
  if (length(ix)>0) {
	inits = getInits(x[ix], reINITS)
	x = x[-ix]
  } else {
	inits = ""
  }
##print(inits);print(x)


  len = length(x)
  cat(sprintf("%s;\n", x[2:(len-1)]), file="eqn__.txt")
  nrhs = integer(1)
  RxODE::rxReq("dparser");
  x = .C("parse_pars", "eqn__.txt", "foo__.txt", nrhs, as.integer(FALSE))
  nrhs = x[[3]]
  foo = paste(readLines("foo__.txt"), collapse="\n")
  nm = system.file("", package = "nlmixr");
  brew(text=c(saem_cmt_str, saem_ode_str)[1+is.ode], output=saem.cpp)
  unlink(c("eqn__.txt", "foo__.txt"))

  ##gen Markevars
  ## if(is.win) x = gsub("\\\\", "/", utils::shortPathName(x))
  ## x = sub("/nlmixr", "", x)
  .lib=  if(is.ode) model$cmpMgr$dllfile else ""

  make_str = 'PKG_CXXFLAGS=%s\nPKG_LIBS=%s $(BLAS_LIBS) $(LAPACK_LIBS) $(FLIBS)\n'
  make_str = sprintf(make_str, nmxInclude(c("nlmixr","StanHeaders","Rcpp","RcppArmadillo","RcppEigen","BH")), .lib)
  cat(make_str, file="Makevars")

  shlib = sprintf('R CMD SHLIB %s -o %s', saem.cpp, saem.dll)
  ## shlib = sprintf(shlib, system.file("include/neldermead.cpp", package = "nlmixr"))
  system(shlib)
  file.copy(file.path(.wd, saem.dll), file.path(lwd, saem.dll));
  file.copy(file.path(.wd, saem.cpp), file.path(lwd, saem.cpp));
  setwd(lwd);
  saem.dll <- file.path(lwd, saem.dll);

  if(is.ode) dyn.load(model$cmpMgr$dllfile)
  `.DLL` <- dyn.load(saem.dll);
  fn.pred <- sourceCppFunction(function(a,b,c) {}, FALSE, `.DLL`, 'dopred')
  fn1 <- sourceCppFunction(function(a) {}, FALSE, `.DLL`, 'saem_fit')
  fn <- eval(bquote(function(a, b, c){
      if (missing(b) && missing(c)){
          cur.fn <- .(fn1)
          ret <- cur.fn(a)
          attr(ret, "dopred") <- .(fn.pred);
          return(ret);
      } else {
          cur.fn <- .(fn.pred)
          return(cur.fn(a, b, c));
      }
  }))
  attr(fn, "form") = if (is.ode) "ode" else "cls"
  attr(fn, "neq") = neq
  attr(fn, "nlhs") = nlhs
  attr(fn, "nrhs") = nrhs
  attr(fn, "saem.dll") = saem.dll
  reg.finalizer(env, saem.cleanup, onexit=TRUE); ## remove dlls on gc or proper exit of R.
  fn
}

##' Cleanup saem_fit environment by removing dll after the object is no logner used by R.
##'
##' @param env Environment where cleanup needs to occur.
##' @author Matthew L. Fidler
##' @export
saem.cleanup <- function(env){
    try({dyn.unload(env$saem.dll)}, silent=TRUE);
    unlink(env$saem.dll);
    unlink(env$saem.cpp);
}

parfn.list = c(
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
lincmt = function(ncmt, oral=T, tlag=F, infusion=F, parameterization=1) {
#ncmt=1; oral=T; tlag=F; parameterization=1
	#print(as.list(match.call()))
	if (infusion) {
		oral = F
		tlag = F
	}

	#master par list
	pm = list(
		c("CL", "V", "KA", "TLAG"),
		c("CL", "V", "CLD", "VT", "KA", "TLAG"),
		c("CL", "V", "CLD", "VT", "CLD2", "VT2", "KA", "TLAG"),
		c("KE", "V", "KA", "TLAG"),
		c("KE", "V", "K12", "K21", "KA", "TLAG"),
		c("KE", "V", "K12", "K21", "K13", "K31", "KA", "TLAG")
	)
	dim(pm)=c(3,2)

	pars = pm[[ncmt, parameterization]]
	if (!tlag) pars[2*ncmt+2] = "0"
	if (!oral) pars = pars[1:(2*ncmt)]
	npar = length(pars)
	s = sprintf("params(%d) = %s;", 1:npar-1, pars)
	pars = paste(c(sprintf("VectorXd params(%d);", 2*ncmt+2), s), collapse="\n")
	attr(pars, "calls") = list(ncmt=ncmt, oral=oral, tlag=tlag, infusion=infusion, parameterization=parameterization)
	ix = (parameterization-1)*9+(ncmt-1)*3+oral+tlag+1
	attr(pars, "default.pars") = parfn.list[[ix]]
	pars
}


#' Plot an SAEM model fit
#'
#' Plot an SAEM model fit
#'
#' @param x a saemFit object
#' @param ... others
#' @return a list
#' @export
plot.saemFit = function(x, ...)
{
    fit = x
    saem.cfg = attr(fit, "saem.cfg")
    dat = as.data.frame(saem.cfg$evt)
    dat = cbind(dat[dat$EVID == 0, ], DV = saem.cfg$y)

    df = rbind(cbind(dat, grp = 1), cbind(dat, grp = 2))
    dopred <- attr(x, "dopred");
    yp = dopred(fit$mpost_phi, saem.cfg$evt, saem.cfg$opt)
    df$DV[df$grp == 2] = yp

    p4 = ggplot(subset(df, grp==1), aes(TIME, DV)) +
        geom_point() +
        facet_wrap(~ID) +
        geom_line(aes(TIME, DV), subset(df, grp==2), col='red')

    df = cbind(dat, IPRED=yp)
    df$IRES = df$DV - df$IPRED

    p2 = ggplot(df, aes(IPRED, DV)) +
        geom_point() +
        geom_abline(intercept = 0, slope=1, col='red')

    p3 = ggplot(df, aes(IPRED, IRES)) +
        geom_point() +
        geom_abline(intercept = 0, slope=0, col='red')

    m = x$par_hist
    df = data.frame(
      val=as.vector(m),
      par=rep(1:ncol(m), each=nrow(m)),
      iter=rep(1:nrow(m), ncol(m))
    )

    p1 = ggplot2::ggplot(df, aes(iter, val)) +
        ggplot2::geom_line() +
        ggplot2::facet_wrap(~par, scales = "free_y")

    print(p1)
    print(p2)
    print(p3)
    print(p4)
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
#' @details
#'    Fit a generalized nonlinear mixed-effect model by he Stochastic
#'    Approximation Expectation-Maximization (SAEM) algorithm
#'
#' @author Wenping Wang
#' @examples
#' \dontrun{
#' library(nlmixr)
#'
#' #ode <- "d/dt(depot) =-KA*depot;
#' #        d/dt(centr) = KA*depot - KE*centr;"
#' #m1 = RxODE(ode, modName="m1")
#' #ode <- "C2 = centr/V;
#' #      d/dt(depot) =-KA*depot;
#' #      d/dt(centr) = KA*depot - KE*centr;"
#' #m2 = RxODE(ode, modName="m2")
#'
#' #Note: only use the '=' assignment, not the '<-' at this point
#'
#' PKpars = function()
#' {
#'   CL = exp(lCL)
#'   V  = exp(lV)
#'   KA = exp(lKA)
#'   KE = CL / V
#'   #initCondition = c(0, KE - CL/V)
#' }
#' PRED = function() centr / V
#' PRED2 = function() C2
#'
#' saem_fit <- gen_saem_user_fn(model=lincmt(ncmt=1, oral=T))
#' #saem_fit <- gen_saem_user_fn(model=m1, PKpars, pred=PRED)
#' #saem_fit <- gen_saem_user_fn(model=m2, PKpars, pred=PRED2)
#'
#'
#' #--- saem cfg
#' nmdat = read.table("theo_sd.dat",  head=T)
#' model = list(saem_mod=saem_fit, covars="WT")
#' inits = list(theta=c(.05, .5, 2))
#' cfg   = configsaem(model, nmdat, inits)
#' cfg$print = 50
#'
#' #cfg$Rfn = nlmixr:::Ruser_function_cmt
#' #dyn.load("m1.d/m1.so");cfg$Rfn = nlmixr:::Ruser_function_ode
#' fit = saem_fit(cfg)
#' df = simple.gof(fit)
#' xyplot(DV~TIME|ID, df, type=c("p","l"), lwd=c(NA,1), pch=c(1,NA), groups=grp)
#' fit
#' }
#' @export
configsaem = function(model, data, inits,
	mcmc=list(niter=c(200,300), nmc=3, nu=c(2,2,2)),
	ODEopt = list(atol=1e-8, rtol=1e-6, stiff=1, transit_abs=0),
	seed=99)
#data =  saem.data; mcmc = saem.mcmc; model = saem.model; inits = saem.inits
{
  set.seed(seed)
  data = list(nmdat=data)

  neq  = attr(model$saem_mod, "neq")
  nlhs = attr(model$saem_mod, "nlhs")
  opt = optM = c(list(neq=neq, nlhs=nlhs, inits=numeric(neq)), ODEopt)
  model$N.eta = attr(model$saem_mod, "nrhs")

  if (is.null(model$log.eta)) model$log.eta = rep(TRUE, model$N.eta)
  if (is.null(model$omega)) model$omega = diag(model$N.eta)
  if (is.null(model$res.mod)) model$res.mod = 1
  if (is.null(inits$omega)) inits$omega = rep(1, model$N.eta)*4
  if (is.null(inits$ares)) inits$ares = 10
  if (is.null(inits$bres)) inits$bres = 1
  if (is.null(mcmc$print)) mcmc$print=1
  inits$theta = matrix(inits$theta, byrow=T, ncol=model$N.eta)
  model$cov.mod=1-is.na(inits$theta)
  data$N.covar=nrow(inits$theta)-1
  inits$theta[is.na(inits$theta)]=0


  ###  FIXME
  mcmc$stepsize=0:1;
  mcmc$burn.in=300

  ###  FIXME: chk covars as char vec
  wh = setdiff(model$covars, names(data$nmdat))
  if (length(wh)) {
  	msg = paste0("covariate(s) not found: ", paste(wh, collapse=", "))
  	stop(msg)
  }
  s = subset(data$nmdat, EVID==0)
  data$data = as.matrix(s[,c("ID", "TIME", "DV", model$covars)])

  nphi = model$N.eta
  mcov = model$cov.mod
  covstruct = model$omega

  check = sum((covstruct-t(covstruct)) != 0)
  if (check) stop("illegal covstruct")
  check = nphi-dim(covstruct)[1]
  if (check) stop("nphi and covstruct dim mismatch")

  check = prod(mcov[1,])
  if (check==0) stop("structural par(s) absent")
  check = nphi-dim(mcov)[2]
  if (check) stop("nphi and ncol(mcov) mismatch")
  check = sum(dim(inits$theta)-dim(mcov) != 0)
  if (check) stop("initial theta's and mcov dim mismatch")
  check = data$N.covar+1 - dim(mcov)[1]
  if (check) stop("dim mcov and N.covar mismatch")

  check = length(model$log.eta)-nphi
  if (check) stop("jlog length and nphi mismatch")

  check = length(inits$omega)-nphi
  if (check) stop("length of omega inits and nphi mismatch")

  #check = mcmc$burn.in>sum(mcmc$niter)
  #if (check) stop("#burn-in exceeds niter")

  check = prod(is.element(covstruct, c(0,1)))
  if (check==0) warning("non-zero value(s) in covstruct set to 1")
  covstruct[covstruct!=0] = 1

  check = prod(is.element(mcov, c(0,1)))
  if (check==0) warning("non-zero value(s) in mcov set to 1")
  mcov[mcov!=0] = 1

  check = sum(inits$theta[1, model$log.eta] <= 0)
  if (check) stop("illegal initial theta's")
  check = sum(inits$omega<=0)
  if (check) stop("illegal initial omega")
  #check = inits$sigma2<=0
  #if (check) stop("illegal initial sigma2")
  check = sum(diag(covstruct)==1)
  if (!check) stop("0 ETA's")


  y = data$data[,"DV"]
  id = data$data[,"ID"]
  ntotal = length(id)
  N = length(unique(id))
  covariables = if(is.null(model$covars)) NULL else tapply(data$data[, model$covars], id, unique)
  nb_measures = table(id)
  ncov = data$N.covar + 1
  nmc = mcmc$nmc
  nM  = mcmc$nmc*N
  yM = rep(y, nmc)
  mlen = max(nb_measures)
  io = t(sapply(nb_measures, function(x) rep(1:0, c(x,mlen-x))))
  ix = rep(1:dim(io)[1], nmc)
  ioM = io[ix,]
  indioM = grep(1, t(ioM)) - 1

  dat = data$nmdat[,c("ID", "TIME", "EVID", "AMT")]		## CHECKME
  form = attr(model$saem_mod, "form")
  infusion = max(dat$EVID)>10000
  if (form=="cls" && infusion) {
	  dat = fmt_infusion_data(dat)
  } else {
	  dat$DUR = -1
  }
  evt = dat
  evt$ID = evt$ID -1
  evtM = evt[rep(1:dim(evt)[1], nmc),]
  evtM$ID = cumsum(c(FALSE, diff(evtM$ID) != 0))


  i1 = grep(1, diag(covstruct))
  i0 = grep(0, diag(covstruct))
  nphi1 = sum(diag(covstruct))
  nphi0 = nphi - nphi1
  na=length(mcmc$stepsize)
  nlambda1 = sum(mcov[,i1])
  nlambda0 = sum(mcov[,i0])
  nlambda=nlambda1+nlambda0
  nd1 =nphi1+nlambda1+1;
  nd2 =nphi1+nlambda1+nlambda0;
  nb_param = nd2+1
  Mcovariables = cbind(rep(1, N), covariables)[,1:nrow(mcov)]
  dim(Mcovariables) = c(length(Mcovariables)/nrow(mcov), nrow(mcov))	#FIXME

  jlog1 = grep(T, model$log.eta)
  jcov = grep(T, apply(mcov, 1, sum)>0)
  covstruct1=covstruct[i1,i1]; dim(covstruct1) = c(nphi1, nphi1)
  ind_cov = grep(1, mcov[mcov>0])

  mcov1 = matrix(mcov[, i1], ncol=length(i1))
  mcov0 = matrix(mcov[, i0], nrow=nrow(mcov), ncol=length(i0))
  ind_cov1 = grep(1, mcov1[mcov1>0]) - 1
  ind_cov0 = grep(1, mcov0[mcov0>0]) - 1

  pc=apply(mcov,2,sum)
  ipc= cumsum(c(0, pc[1:(nphi-1)]))+1
  ipcl1=ipc[jlog1]
  for (x in jlog1) inits$theta[1,x] = log(inits$theta[1,x])


  idx = as.vector(mcov1>0)
  COV1 = Mcovariables[, row(mcov1)[idx]]
  dim(COV1) = c(N, sum(idx))
  COV21 = crossprod(COV1)

  x = mcov1 * col(mcov1)
  x = sapply(x[idx], function(x)
  {
    ret = rep(0, nphi1)
    ret[x] = 1
    ret
  }
  )
  LCOV1 = t(x)
  dim(LCOV1) = c(nlambda1, nphi1)
  pc1 = apply(LCOV1, 2, sum)

  x1 = diag(sum(idx))
  diag(x1) = inits$theta[,i1][idx]
  MCOV1 = x1 %*% LCOV1
  jcov1 = grep(1,LCOV1) - 1


  idx = as.vector(mcov0>0)
  COV0 = Mcovariables[, row(mcov0)[idx]]
  dim(COV0) = c(N, sum(idx))
  COV20 = crossprod(COV0)

  x = mcov0 * col(mcov0)
  x = sapply(x[idx], function(x)
  {
    ret = rep(0, nphi0)
    ret[x] = 1
    ret
  }
  )
  LCOV0 = t(x)
  dim(LCOV0) = c(nlambda0, nphi0)

  x1 = diag(sum(idx))
  diag(x1) = inits$theta[,i0][idx]
  if (dim(x1)[1]>0) {MCOV0 = x1 %*% LCOV0} else {MCOV0 = matrix(x1, nrow=0, ncol=dim(LCOV0)[2])}
  jcov0 = grep(1,LCOV0) - 1

  mprior_phi1 = Mcovariables %*% inits$theta[,i1]
  mprior_phi0 = Mcovariables %*% inits$theta[,i0]

  Gamma2_phi1=diag(nphi1); diag(Gamma2_phi1) = inits$omega[i1]
  Gamma2_phi0=diag(nphi0); diag(Gamma2_phi0) = inits$omega[i0]

  phiM = matrix(0, N, nphi)
  phiM[,i1] = mprior_phi1
  phiM[,i0] = mprior_phi0
  phiM = phiM[rep(1:N, nmc),]
  phiM = phiM + matrix(rnorm(phiM), dim(phiM)) %*% diag(sqrt(inits$omega))

  mc.idx = rep(1:N, nmc)
  statphi = sapply(1:nphi, function(x) tapply(phiM[,x], mc.idx, mean))
  statphi11 = statphi[,i1]; dim(statphi11) = c(N, length(i1))
  statphi01 = statphi[,i0]; dim(statphi01) = c(N, length(i0))
  statphi12 = crossprod(phiM[,i1])
  statphi02 = crossprod(phiM[,i0])

  #x = mcov *cumsum(mcov)
  #x1 = cbind(x[,i1], x[,i0])
  #indiphi = order(x1[x1>0])	#FINDME

  niter = sum(mcmc$niter)
  niter_phi0 = niter*.5
  nb_sa = mcmc$niter[1]*.75
  nb_correl = mcmc$niter[1]*.75
  va=mcmc$stepsize;
  vna=mcmc$niter;
  na=length(va);
  pas=1/(1:vna[1])^va[1] ;
  for (ia in 2:na)
  {
  	end=length(pas);
  	k1=pas[end]^(-1/va[ia]);
  	pas=c(pas, 1/((k1+1):(k1+vna[ia]))^va[ia] );
  }
  pash=c(rep(1,mcmc$burn.in), 1/(1:niter));
  minv = rep(1e-20, nphi)

  i1 = i1 - 1
  i0 = i0 - 1

  cfg=list(
    nu=mcmc$nu,
    niter=niter,
    nb_sa=nb_sa,
    nb_correl=nb_correl,
    niter_phi0=niter_phi0,
    nmc=nmc,
    coef_phi0=.9638,    #FIXME
    rmcmc=.5,
    coef_sa=.95,
    pas=pas,
    pash=pash,
    minv=minv,

    N=N,
    ntotal=ntotal,
    y=y,
    yM=yM,
    phiM=phiM,
    evt=as.matrix(evt),
    evtM=as.matrix(evtM),
    mlen=mlen,
    indioM=indioM,

    pc1=pc1,
    covstruct1=covstruct1,
    Mcovariables=Mcovariables,

    i1=i1,
    i0=i0,
    nphi1=nphi1,
    nphi0=nphi0,
    nlambda1=nlambda1,
    nlambda0=nlambda0,
    COV0=COV0,
    COV1=COV1,
    COV20=COV20,
    COV21=COV21,
    LCOV0=LCOV0,
    LCOV1=LCOV1,
    MCOV0=MCOV0,
    MCOV1=MCOV1,
    Gamma2_phi0=Gamma2_phi0,
    Gamma2_phi1=Gamma2_phi1,
    mprior_phi0=mprior_phi0,
    mprior_phi1=mprior_phi1,
    jcov0=jcov0,
    jcov1=jcov1,
    ind_cov0=ind_cov0,
    ind_cov1=ind_cov1,
    statphi11=statphi11,
    statphi01=statphi01,
    statphi02=statphi02,
    statphi12=statphi12,
    res.mod = model$res.mod,
    ares=inits$ares,
    bres=inits$bres,
    opt=opt,
    optM=optM,
    print=mcmc$print,
    par.hist = matrix(0, sum(niter), nlambda1+nlambda0+nphi1+1+(model$res.mod>2))
  )
}

reINITS = "^\\s*initCondition\\s*=\\s*c\\((?<inits>.+)\\)\\s*$"
reDATAPAR = "^\\s*ParamFromData\\s*=\\s*c\\((?<inits>.+)\\)\\s*$"

getInits = function(x, re, collapse=TRUE) {
	#x = "initCondition = c(1,2,3,4)"
	#re = "^\\s*initCondition\\s*=\\s*c\\((?<inits>.+)\\)\\s*$"
	m = regexpr(re, x, perl=T)
	if (m<0) stop("invalid initCondition input")
	start = attr(m,"capture.start")
	len = attr(m,"capture.length")
	inits = substring(x, start, start+len-1)
	s = strsplit(inits, ",")[[1]]

	if (collapse)
		paste0(sprintf("\tinits[%d] = %s;", 1:length(s)-1, s), collapse="\n")
	else s
}

#' Print an SAEM model fit summary
#'
#' Print an SAEM model fit summary
#'
#' @param object a saemFit object
#' @param ... others
#' @return a list
#' @export
summary.saemFit = function(object, ...)
{
	fit = object ##Rcheck hack

	th = fit$Plambda
	nth = length(th)
	H = solve(fit$Ha[1:nth,1:nth])
	se = sqrt(diag(H))

	m =  cbind(exp(th), th, se)	#FIXME
	## lhsVars = scan("LHS_VARS.txt", what="", quiet=T)
	## if (length(lhsVars)==nth) dimnames(m)[[1]] = lhsVars
	dimnames(m)[[2]] = c("th", "log(th)", "se(log_th)")
	cat("THETA:\n")
	print(m)
	cat("\nOMEGA:\n")
	print(fit$Gamma2_phi1)
	if (any(fit$sig2==0)) {
		cat("\nSIGMA:\n")
		print(max(fit$sig2^2))
	} else {
		cat("\nARES & BRES:\n")
		print(fit$sig2)
	}

	invisible(list(theta=th, se=se, H=H, omega=fit$Gamma2_phi1, eta=fit$mpost_phi))
}

#' Print an SAEM model fit summary
#'
#' Print an SAEM model fit summary
#'
#' @param x a saemFit object
#' @param ... others
#' @return a list
#' @export
print.saemFit = function(x, ...)
{
	fit = x ##Rcheck hack

	th = fit$Plambda
	nth = length(th)
	H = solve(fit$Ha[1:nth,1:nth])
	se = sqrt(diag(H))

	m =  cbind(exp(th), th, se)	#FIXME
	## lhsVars = scan("LHS_VARS.txt", what="", quiet=T)
	## if (length(lhsVars)==nth) dimnames(m)[[1]] = lhsVars
	dimnames(m)[[2]] = c("th", "log(th)", "se(log_th)")
	cat("THETA:\n")
	print(m)
	cat("\nOMEGA:\n")
	print(fit$Gamma2_phi1)
	if (any(fit$sig2==0)) {
		cat("\nSIGMA:\n")
		print(max(fit$sig2^2))
	} else {
		cat("\nARES & BRES:\n")
		print(fit$sig2)
	}

	invisible(list(theta=th, se=se, H=H, omega=fit$Gamma2_phi1, eta=fit$mpost_phi))
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
#' @param seed seed for random number generator
#' @details
#'    Fit a generalized nonlinear mixed-effect model using the Stochastic
#'    Approximation Expectation-Maximization (SAEM) algorithm
#'
#' @author Matthew Fidler & Wenping Wang
#' @export
#' @examples
#' \dontrun{
#' #ode <- "C2 = centr/V;
#' #      d/dt(depot) =-KA*depot;
#' #      d/dt(centr) = KA*depot - KE*centr;"
#' #m1 = RxODE(ode, modName="m2")
#'
#' PKpars = function()
#' {
#'   CL = exp(lCL)
#'   V  = exp(lV)
#'   KA = exp(lKA)
#'   KE = CL / V
#'   #initCondition = c(0, KE - CL/V)
#' }
#' PRED = function() centr / V
#' PRED2 = function() C2
#'
#' dat = theo_sd
#' inits = list(theta=c(.05, .5, 2))
#' fit = saem.fit(model=lincmt(ncmt=1, oral=T), dat, inits)
#' #fit = saem.fit(model=m1, dat, inits, PKpars, pred=PRED)
#' #fit = saem.fit(model=m1, dat, inits, PKpars, pred=PRED2)
#'
#' fit
#' df = plot(fit)
#' }
saem.fit <- function(model, data, inits,
                     PKpars=NULL, pred=NULL,
                     covars=NULL,
                     mcmc = list(niter = c(200, 300), nmc = 3, nu = c(2, 2, 2)),
                     ODEopt = list(atol = 1e-08, rtol = 1e-06, stiff = 1, transit_abs = 0),
                     seed = 99)
{
    UseMethod("saem.fit");
}

##' @rdname saem.fit
##' @export
saem.fit.nlmixr.ui.nlme <- function(model, data, inits,
                                    PKpars=NULL, pred=NULL,
                                    covars=NULL,
                                    mcmc = list(niter = c(200, 300), nmc = 3, nu = c(2, 2, 2)),
                                    ODEopt = list(atol = 1e-08, rtol = 1e-06, stiff = 1, transit_abs = 0),
                                    seed = 99){
    call <- as.list(match.call(expand.dots=TRUE))[-1];
    names(call)[1] <- "object"
    call$est <- "saem";
    return(do.call(getFromNamespace("nlmixr","nlmixr"), call, envir = parent.frame(1)));
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
                           PKpars=NULL, pred=NULL,
                           covars=NULL,
                           mcmc = list(niter = c(200, 300), nmc = 3, nu = c(2, 2, 2)),
                           ODEopt = list(atol = 1e-08, rtol = 1e-06, stiff = 1, transit_abs = 0),
                           seed = 99){
    saem_fit = gen_saem_user_fn(model, PKpars, pred)
    model = list(saem_mod=saem_fit, covars=covars)
    cfg   = configsaem(model, data, inits, mcmc, ODEopt, seed)
    fit = saem_fit(cfg)
    ##dyn.unload("saem_main.dll")
    fit
}

##' @rdname saem.fit
##' @export
saem.fit.default <- function(model, data, inits,
                             PKpars=NULL, pred=NULL,
                             covars=NULL,
                             mcmc = list(niter = c(200, 300), nmc = 3, nu = c(2, 2, 2)),
                             ODEopt = list(atol = 1e-08, rtol = 1e-06, stiff = 1, transit_abs = 0),
                             seed = 99){
    saem_fit = gen_saem_user_fn(model)
    model = list(saem_mod=saem_fit, covars=covars)
    cfg   = configsaem(model, data, inits, mcmc, ODEopt, seed)
    fit = saem_fit(cfg)
                                        #dyn.unload("saem_main.dll")
    fit
}

##' @export
ranef.saemFit <- function(object, ...){
    return(object$eta)
}

##' @export
fixef.saemFit <- function(object, ...){
    return(object$Plambda)
}

focei.theta.saemFit <- function(object, uif, ...){
    ## Get the thetas needed for FOCEi fit.
    if (class(uif) == "function"){
        uif <- nlmixr(uif);
    }
    n <- uif$focei.names
    thetas <- rep(NA, length(n));
    names(thetas) <- n;
    sf <- as.vector(fixed.effects(object))
    theta.trans <- uif$saem.theta.trans
    theta.trans <- theta.trans[!is.na(theta.trans)];
    f <- rep(NA, length(sf))
    for (i in seq_along(theta.trans)){
        ## i = old theta->  theta.trans[i] = saem.theta
        f[i] <- sf[theta.trans[i]];
    }
    ## Now get the names
    n <- uif$theta.names;
    n <- n[is.na(uif$focei.err.type)];
    names(f) <- n;
    for (n in names(f)){
        thetas[n] <- f[n];
    }
    err <- abs(as.vector(object$sig2)) ## abs?
    err.type <- uif$focei.err.type;
    add <- which(err.type == "add")
    prop <- which(err.type == "prop")
    if (length(add) > 0){
        thetas[add] <- sqrt(err[1]); ## Check is this variance (others are sd)
    }
    if (length(prop) > 0){
        thetas[prop] <- sqrt(err[2]); ## Check is this variance (others are sd)
    }
    return(thetas)
}

focei.eta.saemFit <- function(object, uif, ...){
    if (class(uif) == "function"){
        uif <- nlmixr(uif);
    }
    ## Reorder based on translation
    eta.trans <- uif$saem.eta.trans
    ## orig eta ->  new eta
    df <- as.data.frame(uif$ini);
    eta <- df[!is.na(df$neta1), ];
    len <- length(eta$name)
    cur.lhs <- character()
    cur.rhs <- numeric()
    ome <- character()
    cur.ome <- object$Gamma2_phi1;
    for (i in seq_along(eta$name)){
        last.block <- FALSE;
        if (i == len){
            last.block <- TRUE
        } else if (eta$neta1[i + 1] == eta$neta2[i + 1]){
            last.block <- TRUE
        }
        if (eta$neta1[i] == eta$neta2[i]){
            cur.lhs <- c(cur.lhs, sprintf("ETA[%d]", eta$neta1[i]));
            cur.rhs <- c(cur.rhs, cur.ome[eta.trans[eta$neta1[i]], eta.trans[eta$neta2[i]]]);
            if (last.block){
                ome[length(ome) + 1] <- sprintf("%s ~ %s", paste(cur.lhs, collapse=" + "),
                                                paste(deparse(cur.rhs), collapse=" "));
                cur.lhs <- character();
                cur.rhs <- numeric()
            }
        } else {
            cur.rhs <- c(cur.rhs, cur.ome[eta.trans[eta$neta1[i]], eta.trans[eta$neta2[i]]]);
        }
    }
    ome <- eval(parse(text=sprintf("list(%s)", paste(ome, collapse=","))))
    return(ome)
}

as.focei.saemFit <- function(object, uif, pt=proc.time(), ..., data){
    if (class(uif) == "function"){
        uif <- nlmixr(uif);
    }
    uif.new <- uif;
    fit <- object;
    mat <- random.effects(fit);
    ## Reorder based on translation
    eta.trans <- uif$saem.eta.trans
    ## orig eta ->  new eta
    mat2 <- mat;
    for (i in seq_along(eta.trans)){
        ## i = old eta->  eta.trans[i] = saem.eta
        mat2[,i] <- mat[, eta.trans[i]];
    }
    th <- focei.theta(fit, uif)
    for (n in names(th)){
        uif.new$est[uif.new$name == n] <- th[n];
    }
    ome <- focei.eta(fit, uif);
    init <- list(THTA=th,
                 OMGA=ome)
    saem.time <- proc.time() - pt;
    if (missing(data)){
        stop("Requires Data...")
    } else {
        dat <- data;
    }
    fit.f <- focei.fit(data=dat,
                       inits=init,
                       PKpars=uif$theta.pars,
                       ## par_trans=fun,
                       model=uif$rxode.pred,
                       pred=function(){return(nlmixr_pred)},
                       err=uif$error,
                       lower=uif$focei.lower,
                       upper=uif$focei.upper,
                       theta.names=uif$focei.names,
                       eta.names=uif$eta.names,
                       control=list(NOTRUN=TRUE,
                                    inits.mat=mat2,
                                    cores=1,
                                    find.best.eta=FALSE));
    ome <- fit.f$omega;
    w <- which(!is.na(uif.new$ini$neta1))
    for (i in w){
        uif.new$ini$est[i] <- ome[uif.new$ini$neta1[i], uif.new$ini$neta2[i]];
    }
    ## enclose the nlme fit in the .focei.env
    env <- attr(fit.f, ".focei.env");
    env$fit$saem <- fit
    tmp <- cbind(data.frame(saem=saem.time["elapsed"]), env$fit$time);
    names(tmp) <- gsub("optimize", "FOCEi Evaulate", names(tmp))
    env$fit$time <- tmp;
    env$uif <- uif;
    env$uif.new <- uif.new;
    class(fit.f) <- c("nlmixr.ui.saem", class(fit.f))
    return(fit.f)
}

#FIXME: coef_phi0, rmcmc, coef_sa
#FIXME: Klog, rho, sa, nmc
#FIXME: N.design
#FIXME: g = gc = 1
#FIXME: ODE inits
#FIXME: Tinf for ODE
                                        #FIXME: chk infusion poor fit

