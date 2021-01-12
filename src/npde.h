#ifndef __NPDE_H__
#define __NPDE_H__

#include "armahead.h"
#define CENS_OMIT 1
#define CENS_EPRED 2
#define CENS_CDF 3

typedef struct {
  arma::mat matsim;
  arma::umat namat;
  arma::mat epredt;
  arma::mat epred;
  arma::mat varsim;
  arma::mat ymat;
  arma::mat ydsim;
  arma::mat yobst;
  arma::mat yobs;
  arma::mat ydobs;
  arma::mat tcomp;
  arma::mat pd;
  arma::mat npde;
  arma::mat eres;
  unsigned int warn = 0;
} calcNpdeInfoId;

#endif

