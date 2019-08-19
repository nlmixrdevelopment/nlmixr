#include <stdio.h>
#include <stdarg.h>
#include <R_ext/Rdynload.h>
#include <Rcpp.h>
using namespace Rcpp;

// This is a thin layer to the dll.  This is to try to avoid corruption

typedef SEXP (*t_dopred)(SEXP in_phi, SEXP in_evt, SEXP in_opt);
typedef SEXP (*t_saem_fit)(SEXP xSEXP);

t_dopred dopred;
t_saem_fit saem_fit;

void setCur(std::string cur){
  if (cur == "") stop("Need to specify a saem dll string.");
  std::string saemPredStr = "_" + cur + "_dopred";
  std::string saemFitStr = "_" + cur + "_saem_fit";
  dopred =(t_dopred) R_GetCCallable(cur.c_str(), saemPredStr.c_str());
  saem_fit =(t_saem_fit) R_GetCCallable(cur.c_str(), saemFitStr.c_str());
}

//[[Rcpp::export]]
RObject saemDoPred(RObject in_phi, RObject in_evt, RObject in_opt, CharacterVector cur = ""){
  setCur(as<std::string>(cur[0]));
  return dopred(in_phi, in_evt, in_opt);
}

//[[Rcpp::export]]
RObject saemFit(RObject xSEXP, CharacterVector cur = ""){
  setCur(as<std::string>(cur[0]));
  return saem_fit(xSEXP);
}

