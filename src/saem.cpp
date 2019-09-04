#include <stdio.h>
#include <stdarg.h>
#include <thread>
#include <chrono>
#include <R_ext/Rdynload.h>
#include <RcppArmadillo.h>
#include <RxODE.h>
using namespace Rcpp;

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

