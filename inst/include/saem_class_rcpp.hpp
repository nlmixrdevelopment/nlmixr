// saem_class_rcpp.hpp: population PK/PD modeling library
//
// Copyright (C) 2014 - 2016  Wenping Wang
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
  double sigma2;
  vec yM;
  mat evtM;
  List optM;
};

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
  vec sig2;
  sig2 << ares << bres;
  return sig2;
}
mat get_par_hist(){
  return par_hist;
}
mat get_eta(){
  mat eta = mpost_phi;
  rowvec th = Plambda(span(0, nphi-1)).t();
  for (int i=0; i<nphi; i++) eta.row(i) -= th;
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

  opt = as<List>(x["opt"]);		//CHECKME
  optM = as<List>(x["optM"]);	//CHECKME

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

  statrese=0;

  nlambda1 = as<int>(x["nlambda1"]);
  nlambda0 = as<int>(x["nlambda0"]);
  nlambda = nlambda1 + nlambda0;
  nb_param = nphi1 + nlambda + 1;
  nphi = nphi1+nphi0;

  DYF = zeros<mat>(mlen, nM);
  phi.set_size(N, nphi, nmc);

  //FIXME
  res_mod = as<int>(x["res.mod"]);
  //Rcout << "res_mod: " << res_mod << "\n";
  ares = as<double>(x["ares"]);
  bres = as<double>(x["bres"]);
  if (res_mod==1) bres=0;
  if (res_mod==2) ares=0;
  sigma2 = max(ares*ares, 10.0); //FIXME

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
}

void saem_fit() {
  //arma_rng::set_seed(99);
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
    mx.sigma2 = sigma2;

    // CHG hard coded 20
    int nu1, nu2, nu3;
    if (kiter==0) {
      nu1=20*nu(0); nu2=20*nu(1); nu3=20*nu(2);
    } else {
      nu1=nu(0); nu2=nu(1); nu3=nu(2);
    }

    vec f=user_fn(phiM, evtM, optM);
    vec g = ares + bres*f;
    DYF(indioM)=0.5*(((yM-f)/g)%((yM-f)/g))+log(g);
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

    //  MCMC stochastic approximation
    mat Statphi11=zeros<mat>(N,nphi1);
    mat Statphi01=zeros<mat>(N,nphi0);
    mat Statphi12=zeros<mat>(nphi1,nphi1);
    mat Statphi02=zeros<mat>(nphi0,nphi0);
    double statr=0, resk=0;

    vec D1=zeros<vec>(nb_param);    //CHG!!!
    mat D11=zeros<mat>(nb_param,nb_param);
    mat D2=zeros<mat>(nb_param,nb_param);
    vec resy(nmc); //FIXME
    mat d2logk=zeros<mat>(nb_param,nb_param);

    d2logk(span(0,nlambda1-1),span(0,nlambda1-1))=-CGamma21;
    if (nphi0>0) {
    d2logk(span(nlambda1,nlambda-1),span(nlambda1,nlambda-1))=-CGamma20;
    }

    for(int k=0; k<nmc; k++) {
      phi.slice(k)=phiM.rows(span(k*N,(k+1)*N-1));

      Statphi11 += phi.slice(k).cols(i1);
      Statphi01 += phi.slice(k).cols(i0);
      mat phik=phi.slice(k);
      mat phi1k=phik.cols(i1);
      mat phi0k=phik.cols(i0);
      Statphi12=Statphi12+phi1k.t()*phi1k;
      Statphi02=Statphi02+phi0k.t()*phi0k;

      vec fk=user_fn(phik, evt, opt);
      vec gk;

      if (res_mod==1) gk.ones();
      if (res_mod==2) {
		gk = fk;
		gk.elem( find(gk < 1e-8) ).fill(1e-8);  //FIXME by eps
      }
      if (res_mod==1) resk=dot(y-fk, y-fk);
      else if (res_mod==2) resk=dot((y-fk)/gk, (y-fk)/gk);
      else resk = 1;	//FIXME
      statr=statr+resk;
      resy(k) = resk;

      mat dphi1k=phi1k-mprior_phi1;
      mat dphi0k=phi0k-mprior_phi0;
      vec sdg1=sum(dphi1k%dphi1k,0).t()/gamma2_phi1;
      mat Md1=(IGamma2_phi1*(dphi1k.t()*Mcovariables)).t();
      mat Md0=(IGamma2_phi0*(dphi0k.t()*Mcovariables)).t();
      vec d1_mu_phi1=Md1(ind_cov1);    //CHK!! vec or mat
      vec d1_mu_phi0=Md0(ind_cov0);    //CHK!! vec or mat
      vec d1_loggamma2_phi1=0.5*sdg1-0.5*N;
      vec d1_logsigma2;
      d1_logsigma2 << 0.5*resy(k)/sigma2-0.5*ntotal;
      vec d1logk=join_cols(d1_mu_phi1, join_cols(d1_mu_phi0, join_cols(d1_loggamma2_phi1, d1_logsigma2)));
      D1 = D1+d1logk;
      D11= D11+d1logk*d1logk.t();

      vec w2phi=-0.5*sdg1;    //CHK!!!
      for(int j=0, l=0; j<nphi1; j++) {
        for(unsigned int jj=0; jj<pc1(j); jj++) {
          double temp=-dot(COV1.col(l),dphi1k.col(j))/gamma2_phi1(j);
          d2logk(l,nlambda+j)=temp;
          d2logk(nlambda+j,l)=temp;
          l=l+1;
        }
        d2logk(nlambda+j,nlambda+j)=w2phi(j);
      }
      d2logk(nb_param-1,nb_param-1)=-0.5*resy(k)/sigma2;
      D2=D2+d2logk;
    }

    statphi11=statphi11+pas(kiter)*(Statphi11/nmc-statphi11);
    statphi12=statphi12+pas(kiter)*(Statphi12/nmc-statphi12);
    statphi01=statphi01+pas(kiter)*(Statphi01/nmc-statphi01);
    statphi02=statphi02+pas(kiter)*(Statphi02/nmc-statphi02);
    statrese=statrese+pas(kiter)*(statr/nmc-statrese);

    // update parameters
    vec Plambda1, Plambda0;
    Plambda1=inv_sympd(CGamma21)*sum((D1Gamma21%(COV1.t()*statphi11)),1);
    MCOV1(jcov1)=Plambda1;
    if (nphi0>0) {
    Plambda0=inv_sympd(CGamma20)*sum((D1Gamma20%(COV0.t()*statphi01)),1);
    MCOV0(jcov0)=Plambda0;
    }
    mprior_phi1=COV1*MCOV1;
    mprior_phi0=COV0*MCOV0;
    mprior_phi0.set_size(N, nphi0);		// deal w/ nphi0=0

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
      Gamma2_phi0=diagmat(dGamma2_phi0);    //CHK
    }

    double sig2=statrese/ntotal;
    if (res_mod==1) ares=sqrt(sig2);
    else if (res_mod==2) bres=sqrt(sig2);
    else {
      yptr = yM.memptr();
      fptr = f.memptr();
      len = ntotal;
      vec xmin(2);
      double *pxmin = xmin.memptr();

      int n=2, itmax=500, iconv, it, nfcall, iprint=0;
      double start[2]={ares,bres}, step[2]={-.2,-.2}, ynewlo;
      nelder_(obj, n, start, step, itmax, 1.0e-6, 1.0, 2.0, .5,
        &iconv, &it, &nfcall, &ynewlo, pxmin, &iprint);
      ares = ares+pas(kiter)*(pxmin[0]-ares);
      bres = bres+pas(kiter)*(pxmin[1]-bres);
    }
    sigma2=sig2;	//FIXME

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

    vec vcsig2;
    if (res_mod==1) vcsig2 << sigma2;
    if (res_mod==2) vcsig2 << bres;
    if (res_mod==3) vcsig2 << ares << bres;

    Plambda = join_cols(Plambda1, Plambda0);

    par_hist.row(kiter) = join_cols(join_cols(Plambda, Gamma2_phi1.diag()), vcsig2).t();
    if (print>0 && (kiter==0 || (kiter+1)%print==0))
    Rcout << kiter+1
          << ": "
          << par_hist.row(kiter);
  }//kiter
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
  vec y, yM;
  mat evt, evtM;
  mat phiM;
  uvec indioM;
  mat Mcovariables;
  List opt, optM;

  int nphi0, nphi1, nphi;
  mat covstruct1;
  uvec i1, i0;
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

  mat statphi01, statphi02, statphi11, statphi12;
  double statrese;
  double sigma2, ares, bres;
  int res_mod;

  mat DYF;
  cube phi;

  vec L;
  mat Ha, Hb, DDa, DDb;
  mat mpost_phi, cpost_phi;

  mcmcaux mx;

  int print;
  mat par_hist;

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
  vec gc;  //FIXME
  uvec i=mphi.i;

  for (int u=0; u<nu; u++)
  for (int k1=0; k1<mphi.nphi; k1++) {
    mat phiMc=phiM;

    if (method==1)
      phiMc.cols(i)=randn<mat>(mx.nM,mphi.nphi)*mphi.Gamma_phi+mphi.mprior_phiM;
    if (method==2)
      phiMc.cols(i)=phiM.cols(i)+randn<mat>(mx.nM,mphi.nphi)*mphi.Gdiag_phi;
    if (method==3)
      phiMc.col(i(k1))=phiM.col(i(k1))+randn<vec>(mx.nM)*mphi.Gdiag_phi(k1,k1);

    fc=user_fn(phiMc, mx.evtM, mx.optM);
    gc = ares + bres*fc;
    DYF(mx.indioM)=0.5*(((mx.yM-fc)/gc)%((mx.yM-fc)/gc))+log(gc);

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
    if (method<3) break;
  }
}
};

#endif

/*
 * prop err
 * combo err
 * clean up sigma2 related code
 */
