// saem_class_rcpp.hpp: population PK/PD modeling library
//
// Copyright (C) 2014 - 2017  Wenping Wang
//
// This file is part of nlmixr.
//
// nlmixr is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// nlmixr is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with nlmixr.  If not, see <http://www.gnu.org/licenses/>.


#ifndef __SAEM_CLASS_RCPP_HPP__
#define __SAEM_CLASS_RCPP_HPP__
#include <RcppArmadillo.h>
#include <neldermead.hpp>
#define MAXENDPNT 40
using namespace std;
using namespace arma;
using namespace Rcpp;


struct mcmcphi {
  int nphi;
  uvec i;
  mat Gamma_phi;
  mat Gdiag_phi;
  mat IGamma2_phi;
  mat mprior_phiM;
};

struct mcmcaux {
  int nM;
  uvec indioM;
  //double sigma2;  //not needed?
  vec yM;
  mat evtM;
  List optM;
};


uvec getObsIdx(umat m) {
    uvec x;
	x.set_size(0);

	for (unsigned int b=0; b<m.n_rows; ++b) {
	    uvec i=linspace<uvec>(m(b,0), m(b,1), m(b,1) - m(b,0) + 1);
    	x = join_cols(x, i);
    }

	return x;
}


// class def starts
class SAEM {
typedef vec (*user_funct) (const mat&, const mat&, const List&);

public:

SAEM() {
  user_fn = NULL;
}

~SAEM() {}

void set_fn(user_funct f) {
  user_fn = f;
}

mat get_resMat(){
  mat m(nendpnt,2);
  m.col(0) = ares;
  m.col(1) = bres;
  return m;
}
mat get_mprior_phi(){
  mat m = mpost_phi;
  m.cols(i1) = mprior_phi1;
  return m;
}
mat get_mpost_phi(){
  return mpost_phi;
}
mat get_Plambda(){
  return Plambda;
}
mat get_Gamma2_phi1(){
  return Gamma2_phi1;
}
mat get_Ha(){
  return Ha;
}
vec get_sig2(){
  return vcsig2;                                       //FIXME: regression due to multiple endpnts?
}
mat get_par_hist(){
  return par_hist;
}

mat get_eta(){
  mat eta = mpost_phi.cols(i1);
  eta -= mprior_phi1;
  return eta;
}

void inits(List x) {
  nmc = as<int>(x["nmc"]);
  nu = as<uvec>(x["nu"]);
  niter = as<int>(x["niter"]);
  nb_correl = as<int>(x["nb_correl"]);
  niter_phi0 = as<int>(x["niter_phi0"]);
  coef_phi0 = as<double>(x["coef_phi0"]);
  nb_sa = as<int>(x["nb_sa"]);
  coef_sa = as<double>(x["coef_sa"]);
  rmcmc = as<double>(x["rmcmc"]);
  pas = as<vec>(x["pas"]);
  pash = as<vec>(x["pash"]);
  minv = as<vec>(x["minv"]);

  N = as<int>(x["N"]);
  ntotal = as<int>(x["ntotal"]);
  y  = as<vec>(x["y"]);
  yM = as<vec>(x["yM"]);
  evt  = as<mat>(x["evt"]);
  evtM = as<mat>(x["evtM"]);
  phiM = as<mat>(x["phiM"]);
  indioM = as<uvec>(x["indioM"]);
  int mlen = as<int>(x["mlen"]);
  nM = N*nmc;

  opt = as<List>(x["opt"]);                                      //CHECKME
  optM = as<List>(x["optM"]);                                    //CHECKME

  pc1 = as<uvec>(x["pc1"]);
  covstruct1 = as<mat>(x["covstruct1"]);
  Mcovariables = as<mat>(x["Mcovariables"]);

  nphi1 = as<int>(x["nphi1"]);
  i1 = as<uvec>(x["i1"]);
  Gamma2_phi1 = as<mat>(x["Gamma2_phi1"]);
  mprior_phi1 = as<mat>(x["mprior_phi1"]);
  COV1 = as<mat>(x["COV1"]);
  LCOV1 = as<mat>(x["LCOV1"]);
  COV21 = as<mat>(x["COV21"]);
  MCOV1 = as<mat>(x["MCOV1"]);
  jcov1 = as<uvec>(x["jcov1"]);
  ind_cov1 = as<uvec>(x["ind_cov1"]);
  statphi11 = as<mat>(x["statphi11"]);
  statphi12 = as<mat>(x["statphi12"]);

  nphi0 = as<int>(x["nphi0"]);
  if (nphi0>0) {
  i0 = as<uvec>(x["i0"]);
  Gamma2_phi0 = as<mat>(x["Gamma2_phi0"]);
  mprior_phi0 = as<mat>(x["mprior_phi0"]);
  COV0 = as<mat>(x["COV0"]);
  LCOV0 = as<mat>(x["LCOV0"]);
  COV20 = as<mat>(x["COV20"]);
  MCOV0 = as<mat>(x["MCOV0"]);
  jcov0 = as<uvec>(x["jcov0"]);
  ind_cov0 = as<uvec>(x["ind_cov0"]);
  statphi01 = as<mat>(x["statphi01"]);
  statphi02 = as<mat>(x["statphi02"]);
  }
  fixedIx = as<uvec>(x["fixed.ix"]);

  nlambda1 = as<int>(x["nlambda1"]);
  nlambda0 = as<int>(x["nlambda0"]);
  nlambda = nlambda1 + nlambda0;
  nb_param = nphi1 + nlambda + 1;
  nphi = nphi1+nphi0;
  Plambda.zeros(nlambda);
  ilambda1 = as<uvec>(x["ilambda1"]);
  ilambda0 = as<uvec>(x["ilambda0"]);

  DYF = zeros<mat>(mlen, nM);
  phi.set_size(N, nphi, nmc);

  //FIXME
  nendpnt=as<int>(x["nendpnt"]);
  ix_sorting=as<uvec>(x["ix_sorting"]);
  ys = y(ix_sorting);    //ys: obs sorted by endpnt
  ysM=as<vec>(x["ysM"]);
  y_offset=as<uvec>(x["y_offset"]);
  res_mod = as<vec>(x["res.mod"]);
  ares = as<vec>(x["ares"]);
  bres = as<vec>(x["bres"]);
  ix_endpnt=as<uvec>(x["ix_endpnt"]);
  ix_idM=as<umat>(x["ix_idM"]);
  res_offset=as<uvec>(x["res_offset"]);
  nres = res_offset.max();
  vcsig2.set_size(nres);
  vecares = ares(ix_endpnt);
  vecbres = bres(ix_endpnt);
  for (int b=0; b<nendpnt; ++b) {
    sigma2[b] = 10.0;
    if (res_mod(b)==1)
      sigma2[b] = max(ares(b)*ares(b), 10.0);
    if (res_mod(b)==2)
      sigma2[b] = max(bres(b)*bres(b), 1.0);

    statrese[b] = 0.0;
  }

  print = as<int>(x["print"]);
  par_hist = as<mat>(x["par.hist"]);

  L  = zeros<vec>(nb_param);
  Ha = zeros<mat>(nb_param,nb_param);
  Hb = zeros<mat>(nb_param,nb_param);
  mpost_phi = zeros<mat>(N, nphi);
  cpost_phi = zeros<mat>(N, nphi);

  //handle situation when nphi0=0
  mprior_phi0.set_size(N, nphi0);
  statphi01.set_size(N, nphi0);

  mx.nM = nM;
  mx.yM = yM;
  mx.indioM = indioM;
  mx.evtM = evtM;
  mx.optM = optM;

  distribution=as<int>(x["distribution"]);
  DEBUG=as<int>(x["DEBUG"]);
  phiMFile=as<std::vector< std::string > >(x["phiMFile"]);
  //Rcout << phiMFile[0];

}

void saem_fit() {
  //arma_rng::set_seed(99);
  double double_xmin = 1.0e-200;                               //FIXME hard-coded xmin, also in neldermean.hpp
  ofstream phiFile;
  phiFile.open(phiMFile[0].c_str());

  if (DEBUG>0) Rcout << "initialization successful\n";
  fsave = user_fn(phiM, evtM, optM);
  if (DEBUG>0) Rcout << "initial user_fn successful\n";
  for (unsigned int kiter=0; kiter<(unsigned int)(niter); kiter++) {
    gamma2_phi1=Gamma2_phi1.diag();
    IGamma2_phi1=inv_sympd(Gamma2_phi1);
    D1Gamma21=LCOV1*IGamma2_phi1;
    D2Gamma21=D1Gamma21*LCOV1.t();
    CGamma21=COV21%D2Gamma21;

    gamma2_phi0=Gamma2_phi0.diag();
    IGamma2_phi0=inv_sympd(Gamma2_phi0);
    D1Gamma20=LCOV0*IGamma2_phi0;
    D2Gamma20=D1Gamma20*LCOV0.t();
    CGamma20=COV20%D2Gamma20;

    //    MCMC
    mcmcphi mphi1, mphi0;
    set_mcmcphi(mphi1, i1, nphi1, Gamma2_phi1, IGamma2_phi1, mprior_phi1);
    set_mcmcphi(mphi0, i0, nphi0, Gamma2_phi0, IGamma2_phi0, mprior_phi0);

    // CHG hard coded 20
    int nu1, nu2, nu3;
    if (kiter==0) {
      nu1=20*nu(0); nu2=20*nu(1); nu3=20*nu(2);
    } else {
      nu1=nu(0); nu2=nu(1); nu3=nu(2);
    }

    vec f = fsave;
    vec g = vecares + vecbres % abs(f);                          //make sure g > 0
    g.elem( find( g < double_xmin) ).fill(double_xmin);

    //fsave = f;
    if (distribution == 1) DYF(indioM)=0.5*(((yM-f)/g)%((yM-f)/g))+log(g);
    else
    if (distribution == 2) DYF(indioM)=-yM%log(f)+f;
    else
    if (distribution == 3) DYF(indioM)=-yM%log(f)-(1-yM)%log(1-f);
    else {
        Rcout << "unknown distribution\n";
        return;
    }
    //U_y is a vec of subject llik; summed over obs for each subject
    vec U_y=sum(DYF,0).t();

    if(nphi1>0) {
      vec U_phi;
      do_mcmc(1, nu1, mx, mphi1, DYF, phiM, U_y, U_phi);
      mat dphi=phiM.cols(i1)-mphi1.mprior_phiM;
      U_phi=0.5*sum(dphi%(dphi*IGamma2_phi1),1);
      do_mcmc(2, nu2, mx, mphi1, DYF, phiM, U_y, U_phi);
      do_mcmc(3, nu3, mx, mphi1, DYF, phiM, U_y, U_phi);
    }
    if(nphi0>0) {
      vec U_phi;
      do_mcmc(1, nu1, mx, mphi0, DYF, phiM, U_y, U_phi);
      mat dphi=phiM.cols(i0)-mphi0.mprior_phiM;
      U_phi=0.5*sum(dphi%(dphi*IGamma2_phi0),1);
      do_mcmc(2, nu2, mx, mphi0, DYF, phiM, U_y, U_phi);
      do_mcmc(3, nu3, mx, mphi0, DYF, phiM, U_y, U_phi);
    }
    if (DEBUG>0) Rcout << "mcmc successful\n";
    phiFile << phiM;
    //mat dphi=phiM.cols(i1)-mphi1.mprior_phiM;
    //vec U_phi=0.5*sum(dphi%(dphi*IGamma2_phi1),1);

    //  MCMC stochastic approximation
    mat Statphi11=zeros<mat>(N,nphi1);
    mat Statphi01=zeros<mat>(N,nphi0);
    mat Statphi12=zeros<mat>(nphi1,nphi1);
    mat Statphi02=zeros<mat>(nphi0,nphi0);
    double statr[MAXENDPNT], resk;
    for(int b=0; b<nendpnt; ++b) statr[b]= 0;

    vec D1=zeros<vec>(nb_param);    //CHG!!!
    mat D11=zeros<mat>(nb_param,nb_param);
    mat D2=zeros<mat>(nb_param,nb_param);
    vec resy(nmc);
    mat d2logk=zeros<mat>(nb_param,nb_param);

    d2logk(span(0,nlambda1-1),span(0,nlambda1-1))=-CGamma21;
    if (nphi0>0) {
    d2logk(span(nlambda1,nlambda-1),span(nlambda1,nlambda-1))=-CGamma20;
    }


    vec fsM;
    fsM.set_size(0);
    //integration
    for(int k=0; k<nmc; k++) {
      phi.slice(k)=phiM.rows(span(k*N,(k+1)*N-1));

      Statphi11 += phi.slice(k).cols(i1);
      Statphi01 += phi.slice(k).cols(i0);
      mat phik=phi.slice(k);
      mat phi1k=phik.cols(i1);
      mat phi0k=phik.cols(i0);
      Statphi12=Statphi12+phi1k.t()*phi1k;
      Statphi02=Statphi02+phi0k.t()*phi0k;

      vec fk = fsave(span(k*ntotal, (k+1)*ntotal-1));
      fk = fk(ix_sorting);    //sorted by endpnt
      fsM = join_cols(fsM, fk);
      vec resid_all = ys - fk;
      vec gk, resid;

      //loop thru endpoints here
      for(int b=0; b<nendpnt; ++b) {
        resid = resid_all(span(y_offset(b),y_offset(b+1)-1));
        if (res_mod(b)==2) {
        //double epsilon = std::numeric_limits<double>::epsilon();
          gk = abs(fk(span(y_offset(b),y_offset(b+1)-1)));            //CHK: range & chk resize & .memptr()
          gk.elem( find( gk < double_xmin) ).fill(double_xmin);
          resid = resid/gk;
        }
#if 0
        uvec iix = find(resid>1e9);
        Rcout << b << " " <<iix;
        Rcout << ys(iix) << fk(iix) << gk(iix);
#endif

        if (res_mod(b)<=2)
          resk = dot(resid, resid);
        else
          resk = 1;                                              //FIXME

        statr[b]=statr[b]+resk;
        resy(k) = resk;                                          //FIXME: resy(b,k)?
      }
      if (DEBUG>1) Rcout << "star[] successful\n";

      mat dphi1k=phi1k-mprior_phi1;
      mat dphi0k=phi0k-mprior_phi0;
      vec sdg1=sum(dphi1k%dphi1k,0).t()/gamma2_phi1;
      mat Md1=(IGamma2_phi1*(dphi1k.t()*Mcovariables)).t();
      mat Md0=(IGamma2_phi0*(dphi0k.t()*Mcovariables)).t();
      vec d1_mu_phi1=Md1(ind_cov1);                              //CHK!! vec or mat
      vec d1_mu_phi0=Md0(ind_cov0);                              //CHK!! vec or mat
      vec d1_loggamma2_phi1=0.5*sdg1-0.5*N;
      vec d1_logsigma2;
      d1_logsigma2 << 0.5*resy(k)/sigma2[0]-0.5*ntotal;          //FIXME: sigma2[0], sigma2[b] instead?
      vec d1logk=join_cols(d1_mu_phi1, join_cols(d1_mu_phi0, join_cols(d1_loggamma2_phi1, d1_logsigma2)));
      D1 = D1+d1logk;
      D11= D11+d1logk*d1logk.t();

      vec w2phi=-0.5*sdg1;                                       //CHK!!!
      for(int j=0, l=0; j<nphi1; j++) {
        for(unsigned int jj=0; jj<pc1(j); jj++) {
          double temp=-dot(COV1.col(l),dphi1k.col(j))/gamma2_phi1(j);
          d2logk(l,nlambda+j)=temp;
          d2logk(nlambda+j,l)=temp;
          l=l+1;
        }
        d2logk(nlambda+j,nlambda+j)=w2phi(j);
      }
      d2logk(nb_param-1,nb_param-1)=-0.5*resy(k)/sigma2[0];      //FIXME: sigma2[0], sigma2[b] instead?
      D2=D2+d2logk;
    }//k
    if (DEBUG>0) Rcout << "integration successful\n";

    statphi11=statphi11+pas(kiter)*(Statphi11/nmc-statphi11);
    statphi12=statphi12+pas(kiter)*(Statphi12/nmc-statphi12);
    statphi01=statphi01+pas(kiter)*(Statphi01/nmc-statphi01);
    statphi02=statphi02+pas(kiter)*(Statphi02/nmc-statphi02);
    for(int b=0; b<nendpnt; ++b)
      statrese[b]=statrese[b]+pas(kiter)*(statr[b]/nmc-statrese[b]);

    // update parameters
    vec Plambda1, Plambda0;
    Plambda1=inv_sympd(CGamma21)*sum((D1Gamma21%(COV1.t()*statphi11)),1);
    MCOV1(jcov1)=Plambda1;
    if (nphi0>0) {
    Plambda0=inv_sympd(CGamma20)*sum((D1Gamma20%(COV0.t()*statphi01)),1);
    if (fixedIx.n_elem>0) {
      Plambda0(fixedIx) = MCOV0(jcov0(fixedIx));
    }
    MCOV0(jcov0)=Plambda0;
    }
    mprior_phi1=COV1*MCOV1;
    mprior_phi0=COV0*MCOV0;
    mprior_phi0.set_size(N, nphi0);                              // deal w/ nphi0=0

    mat G1=statphi12/N+mprior_phi1.t()*mprior_phi1/N - statphi11.t()*mprior_phi1/N - mprior_phi1.t()*statphi11/N;
    if (kiter<=(unsigned int)(nb_sa))
      Gamma2_phi1=max(Gamma2_phi1*coef_sa, diagmat(G1));
    else
      Gamma2_phi1=G1;
    Gamma2_phi1=Gamma2_phi1%covstruct1;
    vec Gmin=minv(i1);
    uvec jDmin=find(Gamma2_phi1.diag()<Gmin);
    for(unsigned int jm=0; jm<jDmin.n_elem; jm++)
      Gamma2_phi1(jDmin(jm),jDmin(jm))=Gmin(jDmin(jm));
    if (kiter<=(unsigned int)(nb_correl))
    Gamma2_phi1 = diagmat(Gamma2_phi1);

    if (nphi0>0) {
      if (kiter<=(unsigned int)(niter_phi0)) {
      Gamma2_phi0=statphi02/N+mprior_phi0.t()*mprior_phi0/N - statphi01.t()*mprior_phi0/N - mprior_phi0.t()*statphi01/N;
      Gmin=minv(i0);
      jDmin=find(Gamma2_phi0.diag()<Gmin);
      for(unsigned int jm=0; jm<jDmin.n_elem; jm++)
        Gamma2_phi0(jDmin(jm),jDmin(jm))=Gmin(jDmin(jm));
      dGamma2_phi0=Gamma2_phi0.diag();
      } else
      dGamma2_phi0=dGamma2_phi0*coef_phi0;
      Gamma2_phi0=diagmat(dGamma2_phi0);                         //CHK
    }

    //CHECK the following seg on b & yptr & fptr
    for(int b=0; b<nendpnt; ++b) {
      double sig2=statrese[b]/(y_offset(b+1)-y_offset(b));       //CHK: range
      if (res_mod(b)==1)
        ares(b) = sqrt(sig2);
      else
      if (res_mod(b)==2)
        bres(b) = sqrt(sig2);
      else {
        uvec idx;
        idx = find(ix_endpnt==b);
        vec ysb, fsb;

        ysb = ysM(idx);
        fsb = fsM(idx);

        yptr = ysb.memptr();
        fptr = fsb.memptr();
        len = ysb.n_elem;                                        //CHK: needed by nelder
        vec xmin(2);
        double *pxmin = xmin.memptr();
        int n=2, itmax=50, iconv, it, nfcall, iprint=0;
        double start[2]={sqrt(ares(b)), sqrt(fabs(b))};                  //force are & bres to be positive
        double step[2]={-.2, -.2}, ynewlo;
        nelder_(obj, n, start, step, itmax, 1.0e-4, 1.0, 2.0, .5,        //CHG hard-coded tol
                &iconv, &it, &nfcall, &ynewlo, pxmin, &iprint);
        ares(b) = ares(b) + pas(kiter)*(pxmin[0]*pxmin[0] - ares(b));    //force are & bres to be positive
        bres(b) = bres(b) + pas(kiter)*(pxmin[1]*pxmin[1] - bres(b));    //force are & bres to be positive
      }
      sigma2[b] = sig2;                                          //CHK: sigma2[] use
    }
    vecares = ares(ix_endpnt);
    vecbres = bres(ix_endpnt);
    if (DEBUG>0) Rcout << "par update successful\n";


    //    Fisher information
    DDa=(D1/nmc)*(D1/nmc).t()-D11/nmc-D2/nmc;
    DDb=-D11/nmc-D2/nmc;
    L=L+pash(kiter)*(D1/nmc-L);
    Ha=Ha+pash(kiter)*(DDa- Ha);
    Hb=Hb+pash(kiter)*(DDb- Hb);
    cube phi2 = phi%phi;
    mat sphi1 = sum(phi ,2);
    mat sphi2 = sum(phi2,2);
    mpost_phi=mpost_phi+pash(kiter)*(sphi1/nmc-mpost_phi);
    cpost_phi=cpost_phi+pash(kiter)*(sphi2/nmc-cpost_phi);
    mpost_phi.cols(i0)=mprior_phi0;

    //FIXME: chg according to multiple endpnts; need to chg dim(par_hist)
    for (int b=0; b<nendpnt; ++b) {
		int offset = res_offset[b];
        if (res_mod(b)==1) vcsig2[offset] = sigma2[b];
        if (res_mod(b)==2) vcsig2[offset] = bres(b);
        if (res_mod(b)==3) {
			vcsig2[offset]   = ares(b);
			vcsig2[offset+1] = bres(b);
		}
	}

    Plambda(ilambda1) = Plambda1;
    Plambda(ilambda0) = Plambda0;

    par_hist.row(kiter) = join_cols(join_cols(Plambda, Gamma2_phi1.diag()), vcsig2).t();
    if (print>0 && (kiter==0 || (kiter+1)%print==0))
    Rcout << kiter+1
          << ": "
          << par_hist.row(kiter);
    Rcpp::checkUserInterrupt();
  }//kiter
  phiFile.close();

}


private:

  user_funct user_fn;

  uvec nu;
  int niter;
  int nb_sa;
  int nb_correl;
  int niter_phi0;
  double coef_phi0;
  double rmcmc;
  double coef_sa;
  vec pas, pash;
  vec minv;
  int nmc;
  int nM;

  int ntotal, N;
  vec y, yM, ys;    //ys is y sorted by endpnt
  mat evt, evtM;
  mat phiM;
  uvec indioM;
  mat Mcovariables;
  List opt, optM;

  int nphi0, nphi1, nphi;
  mat covstruct1;
  uvec i1, i0, fixedIx;
  uvec pc1;
  mat COV1, COV0, LCOV1, LCOV0, COV21, COV20, MCOV1, MCOV0;
  mat Gamma2_phi1, Gamma2_phi0, mprior_phi1, mprior_phi0;
  mat IGamma2_phi1, D1Gamma21, D2Gamma21, CGamma21;
  mat IGamma2_phi0, D1Gamma20, D2Gamma20, CGamma20;
  mat Gamma_phi1, Gdiag_phi1, Gamma_phi0, Gdiag_phi0;
  vec gamma2_phi1, gamma2_phi0;
  uvec ind_cov1, ind_cov0, jcov1, jcov0;
  vec dGamma2_phi0;
  vec Plambda;

  int nlambda1, nlambda0, nlambda, nb_param;
  uvec ilambda1, ilambda0;

  mat statphi01, statphi02, statphi11, statphi12;
  double statrese[MAXENDPNT];
  double sigma2[MAXENDPNT];
  vec ares, bres;
  vec vecares, vecbres;
  vec res_mod;

  mat DYF;
  cube phi;

  vec L;
  mat Ha, Hb, DDa, DDb;
  mat mpost_phi, cpost_phi;

  mcmcaux mx;

  int print;
  mat par_hist;
  int distribution;

  int nendpnt;
  uvec ix_endpnt;
  umat ix_idM;
  uvec y_offset;
  uvec res_offset;
  vec vcsig2;
  int nres;
  uvec ix_sorting;
  vec ysM;
  vec fsave;

  int DEBUG;
  std::vector< std::string > phiMFile;

void set_mcmcphi(mcmcphi &mphi1,
                 const uvec i1,
                 const int nphi1,
                 const mat Gamma2_phi1,
                 const mat IGamma2_phi1,
                 const mat mprior_phi1) {
    mphi1.i = i1;
    mphi1.nphi = nphi1;
    mphi1.Gamma_phi=chol(Gamma2_phi1);
    mphi1.IGamma2_phi = IGamma2_phi1;
    mphi1.Gdiag_phi.zeros(nphi1, nphi1);
    mphi1.Gdiag_phi.diag() = sqrt(Gamma2_phi1.diag())*rmcmc;
    mphi1.mprior_phiM = repmat(mprior_phi1,nmc,1);
}

void do_mcmc(const int method,
             const int nu,
             const mcmcaux &mx,
             const mcmcphi &mphi,
             mat &DYF,
             mat &phiM,
             vec &U_y,
             vec &U_phi) {
  vec fc, Uc_y, Uc_phi, deltu;
  uvec ind;
  vec gc;
  uvec i=mphi.i;
  double double_xmin = 1.0e-200;                               //FIXME hard-coded xmin, also in neldermean.hpp

  for (int u=0; u<nu; u++)
  for (int k1=0; k1<mphi.nphi; k1++) {
    mat phiMc=phiM;

    if (method==1)
      phiMc.cols(i)=randn<mat>(mx.nM,mphi.nphi)*mphi.Gamma_phi+mphi.mprior_phiM;
    if (method==2)
      phiMc.cols(i)=phiM.cols(i)+randn<mat>(mx.nM,mphi.nphi)*mphi.Gdiag_phi;
    if (method==3)
      phiMc.col(i(k1))=phiM.col(i(k1))+randn<vec>(mx.nM)*mphi.Gdiag_phi(k1,k1);

    fc = user_fn(phiMc, mx.evtM, mx.optM);
    gc = vecares + vecbres % abs(fc);                            //make sure gc > 0
    gc.elem( find( gc < double_xmin) ).fill(double_xmin);
    if (distribution == 1) DYF(mx.indioM)=0.5*(((mx.yM-fc)/gc)%((mx.yM-fc)/gc))+log(gc);
    else
    if (distribution == 2) DYF(mx.indioM)=-mx.yM%log(fc)+fc;
    else
    if (distribution == 3) DYF(indioM)=-mx.yM%log(fc)-(1-mx.yM)%log(1-fc);

    Uc_y=sum(DYF,0).t();
    if (method==1) deltu=Uc_y-U_y;
    else {
      mat dphic=phiMc.cols(i)-mphi.mprior_phiM;
      Uc_phi=0.5*sum(dphic%(dphic*mphi.IGamma2_phi),1);
      deltu=Uc_y-U_y+Uc_phi-U_phi;
    }

    ind=find( deltu<-log(randu<vec>(mx.nM)) );
    phiM(ind,i)=phiMc(ind,i);
    U_y(ind)=Uc_y(ind);
    if (method>1) U_phi(ind)=Uc_phi(ind);
    ind = getObsIdx(ix_idM.rows(ind));
    fsave(ind)=fc(ind);
    if (method<3) break;
  }
}
};


// closing for #ifndef __SAEM_CLASS_RCPP_HPP__
#endif



/*
- allow phi to psi xform
- chg g<EPS & gc<EPS, chg to EPS
- distribution & res_mod be nested
- clean up sigma2 related code
- chk & fix resy(k) for multiple endpnts
    * resy(k) is only used in Fisher info concerning sigma2; currently NOT used by nlmixr
- two for loops to convert double to integer for EVID & CMT
    * may declare a umat to hold these to avoid loops for faster speed
    * examine all loops
*/


