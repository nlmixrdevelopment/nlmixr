#ifndef __NPDE_H__
#define __NPDE_H__

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("nlmixr", String)
/* replace pkg as appropriate */
#else
#define _(String) (String)
#endif

#if defined(__cplusplus)
#include "armahead.h"
#include "censResid.h"
#include "res.h"

#define NPDE_CHOL 0
#define NPDE_CHOL_PINV 1
#define NPDE_DECORRELATE_EIGEN  2
#define NPDE_DECORRELATE_EIGEN_PINV  3
#define NPDE_CHOLSE 4
#define NPDE_CHOLSE_PINV 5
#define NPDE_NPD 99

typedef const char *(*rxGetId2_t)(int id);
extern rxGetId2_t rxGetId2;


typedef struct {
  arma::mat matsim;
  arma::umat namat;
  arma::mat epredt;
  arma::mat epred;
  arma::umat obs;
  arma::mat varsim;
  arma::mat ymat;
  arma::mat ymat2;
  arma::mat ydsim;
  arma::mat ydsim2;
  arma::mat yobst;
  arma::mat yobst2;
  arma::mat yobs;
  arma::mat ydobs;
  arma::mat ydobs2;
  arma::mat tcomp;
  arma::mat tcomp2;
  arma::mat pd;
  arma::mat pd2;
  arma::mat npde;
  arma::mat npd;
  arma::mat eres;
  unsigned int warn = 0;
} calcNpdeInfoId;

extern "C" {
#endif

  SEXP _nlmixr_npdeCalc(SEXP npdeSim, SEXP dvIn, SEXP evidIn, SEXP censIn, SEXP limitIn, SEXP npdeOpt);

#if defined(__cplusplus)
}
#endif

#endif
