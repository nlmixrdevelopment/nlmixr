// [[Rcpp::plugins(openmp)]]
#define STRICT_R_HEADER
#include "armahead.h"
#include "utilc.h"
#include <lbfgsb3c.h>

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("nlmixr", String)
/* replace pkg as appropriate */
#else
#define _(String) (String)
#endif

#define NETAs 20
#define NTHETAs 20
#define NSUBs 100
#define min2( a , b )  ( (a) < (b) ? (a) : (b) )
#define max2( a , b )  ( (a) > (b) ? (a) : (b) )
#define innerOde(id) ind_solve(rx, id, rxInner.dydt_liblsoda, rxInner.dydt_lsoda_dum, rxInner.jdum_lsoda, rxInner.dydt, rxInner.update_inis, rxInner.global_jt)
#define predOde(id) ind_solve(rx, id, rxPred.dydt_liblsoda, rxPred.dydt_lsoda_dum, rxPred.jdum_lsoda, rxPred.dydt, rxPred.update_inis, rxPred.global_jt)
#define getCholOmegaInv() (as<arma::mat>(RxODE::rxSymInvCholEnvCalculate(_rxInv, "chol.omegaInv", R_NilValue)))
#define getOmega() (as<NumericMatrix>(RxODE::rxSymInvCholEnvCalculate(_rxInv, "omega", R_NilValue)))
#define getOmegaMat() (as<arma::mat>(RxODE::rxSymInvCholEnvCalculate(_rxInv, "omega", R_NilValue)))
#define getOmegaInv() (as<arma::mat>(RxODE::rxSymInvCholEnvCalculate(_rxInv, "omegaInv", R_NilValue)))
#define getOmegaDet() (as<double>(RxODE::rxSymInvCholEnvCalculate(_rxInv, "log.det.OMGAinv.5", R_NilValue)))
#define getOmegaN() as<int>(RxODE::rxSymInvCholEnvCalculate(_rxInv, "ntheta", R_NilValue))
#define getOmegaTheta() as<NumericVector>(RxODE::rxSymInvCholEnvCalculate(_rxInv, "theta", R_NilValue));
#define setOmegaTheta(x) RxODE::rxSymInvCholEnvCalculate(_rxInv, "theta", x)
#define tbs(x) _powerD(x,    ind->lambda, (int)(ind->yj), ind->logitLow, ind->logitHi)
#define tbsL(x) _powerL(x,   ind->lambda, (int)(ind->yj), ind->logitLow, ind->logitHi)
#define tbsDL(x) _powerDL(x, ind->lambda, (int)(ind->yj), ind->logitLow, ind->logitHi)
#define tbsD(x) _powerDD(x,  ind->lambda, (int)(ind->yj), ind->logitLow, ind->logitHi)
#define _safe_log(a) (((a) <= 0.0) ? log(DBL_EPSILON) : log(a))
// #define _safe_log(a) log(a)
#define _safe_zero(a) ((a) == 0 ? DBL_EPSILON : (a))
//#define _safe_zero(a) (a)
#define _safe_sqrt(a) ((a) <= 0 ? sqrt(DBL_EPSILON) : sqrt(a))
//#define _safe_sqrt(a) sqrt(a)



#define expit(alpha, low, high) _powerDi(alpha, 1.0, 4, low, high)
#define probitInv(alpha, low, high) _powerDi(alpha, 1.0, 6, low, high)

extern "C" {
  typedef void (*S2_fp) (int *, int *, double *, double *, double *, int *, float *,
			 double *, int *);
  typedef void (*n1qn1_fp)(S2_fp simul, int n[], double x[], double f[], double g[],
			   double var[], double eps[], int mode[], int niter[], int nsim[],
			   int imp[], double zm[], int izs[], float rzs[], double dzs[],
			   int id[]);

  n1qn1_fp n1qn1_;

  typedef double optimfn(int n, double *par, void *ex);

  typedef void optimgr(int n, double *par, double *gr, void *ex);

  void lbfgsbRX(int n, int lmm, double *x, double *lower,
                double *upper, int *nbd, double *Fmin, optimfn fn,
                optimgr gr, int *fail, void *ex, double factr,
                double pgtol, int *fncount, int *grcount,
                int maxit, char *msg, int trace, int nREPORT);


  typedef void (*ind_solve_t)(rx_solve *rx, unsigned int cid, t_dydt_liblsoda dydt_lls,
			      t_dydt_lsoda_dum dydt_lsoda, t_jdum_lsoda jdum,
			      t_dydt c_dydt, t_update_inis u_inis, int jt);
  ind_solve_t ind_solve;
  typedef int (*par_progress_t)(int c, int n, int d, int cores, clock_t t0, int stop);
  par_progress_t par_progress;
  typedef rx_solve* (*getRxSolve_t)();
  typedef int (*isRstudio_t)();
  isRstudio_t isRstudio;
  getRxSolve_t getRx;
  typedef const char *(*rxGetId_t)(int id);
  rxGetId_t rxGetId;
  typedef double (*getTime_t)(int idx, rx_solving_options_ind *ind);
  getTime_t getTimeF;
  typedef void (*sortIds_t)(rx_solve* rx, int ini);
  sortIds_t sortIdsF;
}

typedef int (*iniSubjectI_t)(int solveid, int inLhs, rx_solving_options_ind *ind, rx_solving_options *op, rx_solve *rx,
			     t_update_inis u_inis);

iniSubjectI_t iniSubjectI;


bool assignFn_ = false;

extern void lin_cmt_stanC(double *obs_timeD, const int nobs, double *dose_timeD, const int ndose, double *doseD, double *TinfD,
			  double *paramsD, const int oral, const int infusion, const int ncmt, const int parameterization,
			  const int neta, double *fxD, double *dvdxD, double *fpD);

List _rxInv;

// These are focei inner options
typedef struct {
  //
  // std::string estStr;
  // std::string gradStr;
  // std::string obfStr;
  //
  List mvi;
  double *etaUpper = NULL;
  double *etaLower = NULL;
  int *nbdInner = NULL;
  double *geta = NULL;
  double *getahf =NULL;
  double *getahr = NULL;
  double *goldEta = NULL;
  double *gsaveEta = NULL;
  double *gthetaGrad = NULL;
  bool mGthetaGrad = false;
  // n1qn1 specific vectors
  double *gZm = NULL;
  double *gG = NULL;
  double *gVar = NULL;
  double *gX = NULL;

  double *glp = NULL;
  double *ga = NULL;
  double *gB = NULL;
  double *gc = NULL;
  double *gH = NULL;
  double *gH0 = NULL;
  double *gVid = NULL;

  double *likSav = NULL;

  // Integer of ETAs
  unsigned int gEtaGTransN;
  // Where likelihood is saved.

  int *etaTrans = NULL;
  int *etaFD = NULL;
  double eventFD;
  int predNeq;
  int eventType;

  int neta;
  unsigned int ntheta;
  int npars;
  int thetan;
  int omegan;

  int calcGrad;
  int nF;
  int nF2;
  int nG;
  int derivMethod;
  int covDerivMethod;
  int covMethod;
  int derivMethodSwitch;
  double derivSwitchTol;
  double lastOfv;

  double *fullTheta = NULL;
  double *theta = NULL;
  double *thetaGrad = NULL;
  double *initPar = NULL;
  double *scaleC = NULL;
  double scaleC0;
  int *xPar = NULL;
  NumericVector lowerIn;
  double *lower = NULL;
  NumericVector upperIn;
  double *upper = NULL;
  int *nbd = NULL;

  int *fixedTrans = NULL;
  int *thetaTrans = NULL;

  int scaleType;
  int normType;
  double scaleCmin;
  double scaleCmax;
  double c1;
  double c2;
  double scaleTo;
  double epsilon;

  int maxOuterIterations;
  int maxInnerIterations;

  double odeRecalcFactor;
  int maxOdeRecalc;
  int objfRecalN;
  int stickyRecalcN1;
  int stickyRecalcN2;
  int stickyRecalcN;
  int stickyTol=0;

  int nsim;
  int nzm;

  int imp;
  // int printInner;
  int printOuter;


  mat omega;
  mat omegaInv;
  mat cholOmegaInv;
  mat etaM;
  mat etaS;
  mat eta1SD;
  double n;
  double logDetOmegaInv5;

  int    *gillRetC = NULL;
  int    *gillRet = NULL;
  double *gillDf = NULL;
  double *gillDf2 = NULL;
  double *gillErr = NULL;
  double *rEps = NULL;
  double *aEps = NULL;
  double *rEpsC = NULL;
  double *aEpsC = NULL;

  //
  double factr;
  double pgtol;
  double abstol;
  double reltol;
  int lmm;
  int *skipCov = NULL;
  int skipCovN;

  int outerOpt;
  int eigen;
  int scaleObjective;
  double scaleObjectiveTo;
  int initObj;
  double initObjective;
  // Confidence Interval
  double ci;
  double sigdig;
  //
  clock_t t0 = clock();
  int cur = 0;
  int curTick = 0;
  int totTick = 100;
  int useColor;
  double boundTol;
  int printNcol;
  int noabort;
  int interaction;
  double cholSEtol;
  double hessEps;
  double cholAccept;
  double resetEtaSize;
  int didEtaReset;
  double resetThetaSize = R_PosInf;
  double resetThetaFinalSize = R_PosInf;
  int checkTheta;
  int *muRef = NULL;
  int muRefN;
  int resetHessianAndEta;
  int didHessianReset;
  int cholSEOpt=0;
  int cholSECov;
  int fo;
  int covTryHarder;
  // Gill options
  int gillK;
  double gillStep;
  double gillFtol;
  double gillRtol;
  int gillKcov;
  double gillStepCov;
  double gillFtolCov;
  double covSmall;
  int didGill;
  int smatNorm;
  int rmatNorm;
  int covGillF;
  int optGillF;
  int mixDeriv;
  double gradTrim;
  double gradCalcCentralSmall;
  double gradCalcCentralLarge;
  double etaNudge;
  double etaNudge2;
  int didEtaNudge;
  int reducedTol;
  int reducedTol2;
  int repeatGill;
  int repeatGillN;
  int repeatGillMax;
  int curGill;
  int printTop;
  double resetThetaCheckPer;
  int slow = 0;
  double gradProgressOfvTime;
  bool alloc=false;
  bool zeroGrad = false;
  int nfixed=0;
  NumericVector logitThetaHi;
  NumericVector logitThetaLow;
  double badSolveObjfAdj = 100;
  bool didPredSolve = false;
  bool canDoFD  = false;
} focei_options;

focei_options op_focei;

typedef struct {
  int nInnerF;
  int nInnerG;
  double lik[3]; // lik[0] = liklihood; For central difference: lik[1] = lower lik[2] = upper
  double *eta; // Eta includes the ID number for the patient
  double *etahf;
  double *etahr;
  //
  double *thetaGrad; // Theta gradient; Calculated on the individual level for S matrix calculation
  double thVal[2]; // thVal[0] = lower; thVal[2] = upper
  //
  // F and varaibility
  unsigned int nobs;
  unsigned int setup;

  double *saveEta; // Saved when lik[0] is saved.
  double *oldEta;

  // Likilihood gradient
  double llik;
  double *a;
  double *B;
  double *c;
  double *lp;// = mat(neta,1);

  double *g;
  double *H;
  double *H0;
  double *Vid;

  double tbsLik;

  int mode; // 1 = dont use zm, 2 = use zm.
  double *zm;
  double *var;
  double *x;
  unsigned int uzm;
  int doChol=1;
  int doFD=0;
  int doEtaNudge;
  int badSolve=0;
  double curF;
  double curT;
  double *curS;
} focei_ind;

focei_ind *inds_focei = NULL;

// Parameter table
std::vector<int> niter;
std::vector<int> iterType;
std::vector<double> vPar;
std::vector<double> vGrad;
std::vector<int> niterGrad;
std::vector<int> gradType;

extern "C" void rxOptionsFreeFocei(){

  if (op_focei.etaTrans != NULL) R_Free(op_focei.etaTrans);
  op_focei.etaTrans=NULL;

  if (op_focei.fullTheta != NULL) R_Free(op_focei.fullTheta);
  op_focei.fullTheta = NULL;

  if (op_focei.etaUpper != NULL) R_Free(op_focei.etaUpper);
  op_focei.etaUpper = NULL;

  if (op_focei.gillRet != NULL) R_Free(op_focei.gillRet);
  op_focei.gillRet = NULL;

  if (op_focei.gillDf != NULL) R_Free(op_focei.gillDf);
  op_focei.gillDf = NULL;

  if (op_focei.gthetaGrad != NULL && op_focei.mGthetaGrad) R_Free(op_focei.gthetaGrad);
  op_focei.gthetaGrad = NULL;
  op_focei.mGthetaGrad = false;

  if (inds_focei != NULL) R_Free(inds_focei);
  inds_focei=NULL;

  op_focei.alloc = false;
  op_focei.didPredSolve = false;

  focei_options newf;
  op_focei= newf;

  vGrad.clear();
  vPar.clear();
  iterType.clear();
  gradType.clear();
  niter.clear();
  niterGrad.clear();
}

//[[Rcpp::export]]
void freeFocei(){
  rxOptionsFreeFocei();
}

//' Calculate the inverse preconditioning matrix
//'
//' @param Rin The R matrix input
//'
//' @return The inverse preconditioning matrix
//'
//[[Rcpp::export]]
SEXP preCondInv(SEXP Rin){
  // Assumes Rin is symmetric
  arma::vec eigval;
  arma::mat eigvec;
  arma::mat R = as<arma::mat>(Rin);
  bool success = eig_sym(eigval, eigvec, R);
  if (success){
    // Now calculate the norm
    arma::mat eignorm = normalise(eigvec);
    arma::mat v12 = diagmat(1/(abs(eigval)));
    // This is equation 15 in the Aoki 2016 paper
    R = eignorm*v12;
    SEXP out = wrap(R);
    Rf_setAttrib(out, R_DimNamesSymbol, Rf_getAttrib(Rin, R_DimNamesSymbol));
    return out;
  }
  Rcpp::stop("cannot calculate the eigenvectors/eigenvalues required for preconditioning");
  return R_NilValue;
}

typedef struct {
  //
  // std::string estStr;
  // std::string gradStr;
  // std::string obfStr;
  //
  t_dydt dydt = NULL;
  t_calc_jac calc_jac = NULL;
  t_calc_lhs calc_lhs = NULL;
  t_update_inis update_inis = NULL;
  t_dydt_lsoda_dum dydt_lsoda_dum = NULL;
  t_dydt_liblsoda dydt_liblsoda = NULL;
  t_jdum_lsoda jdum_lsoda = NULL;
  t_set_solve set_solve = NULL;
  t_get_solve get_solve = NULL;
  int global_jt = 2;
  int global_mf = 22;
  int global_debug = 0;
  int neq = NA_INTEGER;
} rxSolveF;

rxSolveF rxInner;
rxSolveF rxPred;

void rxUpdateFuns(SEXP trans, rxSolveF *inner){
  const char *lib, *s_dydt, *s_calc_jac, *s_calc_lhs, *s_inis, *s_dydt_lsoda_dum, *s_dydt_jdum_lsoda,
    *s_ode_solver_solvedata, *s_ode_solver_get_solvedata, *s_dydt_liblsoda;
  lib = CHAR(STRING_ELT(trans, 0));
  s_dydt = CHAR(STRING_ELT(trans, 3));
  s_calc_jac = CHAR(STRING_ELT(trans, 4));
  s_calc_lhs = CHAR(STRING_ELT(trans, 5));
  s_inis = CHAR(STRING_ELT(trans, 8));
  s_dydt_lsoda_dum = CHAR(STRING_ELT(trans, 9));
  s_dydt_jdum_lsoda = CHAR(STRING_ELT(trans, 10));
  s_ode_solver_solvedata = CHAR(STRING_ELT(trans, 11));
  s_ode_solver_get_solvedata = CHAR(STRING_ELT(trans, 12));
  s_dydt_liblsoda = CHAR(STRING_ELT(trans, 13));
  inner->global_jt = 2;
  inner->global_mf = 22;
  inner->global_debug = 0;
  if (strcmp(CHAR(STRING_ELT(trans, 1)),"fulluser") == 0){
    inner->global_jt = 1;
    inner->global_mf = 21;
  } else {
    inner->global_jt = 2;
    inner->global_mf = 22;
  }
  inner->calc_lhs =(t_calc_lhs) R_GetCCallable(lib, s_calc_lhs);
  inner->dydt =(t_dydt) R_GetCCallable(lib, s_dydt);
  inner->calc_jac =(t_calc_jac) R_GetCCallable(lib, s_calc_jac);
  inner->update_inis =(t_update_inis) R_GetCCallable(lib, s_inis);
  inner->dydt_lsoda_dum =(t_dydt_lsoda_dum) R_GetCCallable(lib, s_dydt_lsoda_dum);
  inner->jdum_lsoda =(t_jdum_lsoda) R_GetCCallable(lib, s_dydt_jdum_lsoda);
  inner->set_solve = (t_set_solve)R_GetCCallable(lib, s_ode_solver_solvedata);
  inner->get_solve = (t_get_solve)R_GetCCallable(lib, s_ode_solver_get_solvedata);
  inner->dydt_liblsoda = (t_dydt_liblsoda)R_GetCCallable(lib, s_dydt_liblsoda);
}

void rxClearFuns(rxSolveF *inner){
  inner->calc_lhs              = NULL;
  inner->dydt                  = NULL;
  inner->calc_jac              = NULL;
  inner->update_inis           = NULL;
  inner->dydt_lsoda_dum        = NULL;
  inner->jdum_lsoda            = NULL;
  inner->set_solve             = NULL;
  inner->get_solve             = NULL;
  inner->dydt_liblsoda         = NULL;
}

rx_solve* rx;

////////////////////////////////////////////////////////////////////////////////
// n1qn1 functions
uvec lowerTri(mat H, bool diag = false){
  unsigned int d = H.n_rows;
  mat o(d, d, fill::ones);
  if (!diag){
    return find(trimatl(o,-1));
  } else {
    return find(trimatl(o));
  }
}

void updateZm(focei_ind *indF){
  std::fill(&indF->zm[0], &indF->zm[0]+op_focei.nzm,0.0);
  if (!indF->uzm){
    // Udate the curvature to Hessian to restart n1qn1
    int n = op_focei.neta;
    mat L = eye(n, n);
    mat D = mat(n, n, fill::zeros);
    mat H = mat(n, n);
    unsigned int l_n = n * (n + 1)/2;
    vec zmV(l_n);
    std::copy(&indF->zm[0], &indF->zm[0]+l_n, zmV.begin());
    H.elem(lowerTri(H, true)) = zmV;
    if (n == 1) H = D;
    else{
      L.elem(lowerTri(H,false)) = H.elem(lowerTri(H,0));
      D.diag() = H.diag();
      H = L*D*L.t();
    }
    // Hessian -> c.hess
    vec hessV = H.elem(lowerTri(H, true));
    std::copy(hessV.begin(),hessV.end(),&indF->zm[0]);
    indF->uzm = 1;
    indF->mode=2;
  }
}

static inline double getScaleC(int i){
  if (ISNA(op_focei.scaleC[i])){
    switch (op_focei.xPar[i]){
    case 1: // log
      op_focei.scaleC[i]=1.0;
      break;
    case 2: // diag^2
      op_focei.scaleC[i]=fabs(op_focei.initPar[i]);
      break;
    case 3: // exp(diag)
      op_focei.scaleC[i] = 2.0;
      break;
    case 4: // Identity diagonal chol(Omega ^-1)
    case 5: // off diagonal chol(Omega^-1)
      op_focei.scaleC[i] = 2.0*fabs(op_focei.initPar[i]);
      break;
    default:
      op_focei.scaleC[i]= fabs(op_focei.initPar[i]);
      break;
    }
  }
  return min2(max2(op_focei.scaleC[i], op_focei.scaleCmin),op_focei.scaleCmax);
}

////////////////////////////////////////////////////////////////////////////////
// Likelihood for inner functions
static inline double unscalePar(double *x, int i){
  double scaleTo = op_focei.scaleTo, C=getScaleC(i);
  switch(op_focei.scaleType){
  case 1: // normalized
    return x[i]*op_focei.c2+op_focei.c1;
    break;
  case 2: // log vs linear scales and/or ranges
    if (op_focei.normType <= 5){
      scaleTo = (op_focei.initPar[i]-op_focei.c1)/op_focei.c2;
    } else if (scaleTo == 0){
      scaleTo=op_focei.initPar[i];
    }
    return (x[i]-scaleTo)*C + op_focei.initPar[i];
    break;
  case 3: // simple multiplicative scaling
    if (op_focei.scaleTo != 0){
      return x[i]*op_focei.initPar[i]/scaleTo;
    } else {
      return x[i];
    }
    break;
  case 4: // log non-log multiplicative scaling
    if (op_focei.scaleTo > 0){
      switch (op_focei.xPar[i]){
      case 1:
	return (x[i]-scaleTo) + op_focei.initPar[i];
      default:
	return x[i]*op_focei.initPar[i]/scaleTo;
      }
    } else {
      return x[i];
    }
  default:
    if (op_focei.scaleTo > 0){
      return (x[i]-scaleTo)*1 + op_focei.initPar[i];
    } else {
      return x[i];
    }
  }
  return 0;
}

static inline double scalePar(double *x, int i){
  double scaleTo = op_focei.scaleTo, C=getScaleC(i);
  switch(op_focei.scaleType){
  case 1:
    return (x[i]-op_focei.c1)/op_focei.c2;
  case 2:
    if (op_focei.normType <= 5){
      scaleTo = (op_focei.initPar[i]-op_focei.c1)/op_focei.c2;
    } else if (scaleTo == 0){
      scaleTo=op_focei.initPar[i];
    }
    return (x[i]-op_focei.initPar[i])/C + scaleTo;
    break;
  case 3: // simple multiplicative scaling
    if (op_focei.scaleTo > 0){
      return x[i]/op_focei.initPar[i]*op_focei.scaleTo;
    } else {
      return x[i];
    }
    break;
  case 4: // log non-log multiplicative scaling
    if (op_focei.scaleTo > 0){
      switch (op_focei.xPar[i]){
      case 1:
	return (x[i]-op_focei.initPar[i]) + op_focei.scaleTo;
      default:
	return x[i]/op_focei.initPar[i]*op_focei.scaleTo;
      }
    } else {
      return x[i];
    }
  default:
    if (op_focei.scaleTo > 0){
      return (x[i]-op_focei.initPar[i]) + op_focei.scaleTo;
    } else {
      return x[i];
    }
  }
  return 0;
}


void updateTheta(double *theta){
  // Theta is the acutal theta
  unsigned int j, k;
  for (k = op_focei.npars; k--;){
    j=op_focei.fixedTrans[k];
    op_focei.fullTheta[j] = unscalePar(theta, k);
  }
  // Update theta parameters in each individual
  rx = getRx();
  for (int id = rx->nsub; id--;){
    rx_solving_options_ind *ind = &(rx->subjects[id]);
    for (j = op_focei.ntheta; j--;){
      ind->par_ptr[op_focei.thetaTrans[j]] = op_focei.fullTheta[j];
    }
  }
  // Update setOmegaTheta
  if (op_focei.neta > 0) {
    NumericVector omegaTheta(op_focei.omegan);
    std::copy(&op_focei.fullTheta[0] + op_focei.ntheta,
	      &op_focei.fullTheta[0] + op_focei.ntheta + op_focei.omegan,
	      omegaTheta.begin());
    setOmegaTheta(omegaTheta);
    if (op_focei.fo == 1){
      op_focei.omega = getOmegaMat();
    } else {
      op_focei.omegaInv = getOmegaInv();
      op_focei.cholOmegaInv = getCholOmegaInv();
      op_focei.logDetOmegaInv5 = getOmegaDet();
    }
  }
  //Now Setup Last theta
  if (!op_focei.calcGrad){
    // op_focei.estStr=sc + un + ex;
    std::copy(&theta[0], &theta[0] + op_focei.npars, &op_focei.theta[0]);
  }
}

arma::mat cholSE__(arma::mat A, double tol);

static inline void likM2(focei_ind *fInd, double& limit, double&f, double &r) {
  if (R_FINITE(limit) && !ISNA(limit)) {
    // When limit < f, this is M2
    // When limit >=f, the limit is an upper limit (instead of lower limit)
    fInd->llik += -log(1-0.5*(1+erf(((limit<f)*2-1)*(limit-f)/_safe_sqrt(r)/M_SQRT2)));
  }
}

static inline void likCens(focei_ind *fInd, int &cens, double& limit, double&f, double& dv, double &r) {
  fInd->llik += log(0.5*(1+erf(((double)(cens)*(dv-f))/_safe_sqrt(r)/M_SQRT2)));
  if (R_FINITE(limit)  && !ISNA(limit)){
    fInd->llik += -log(1-0.5*(1+erf((double)(cens)*(limit-f)/_safe_sqrt(r)/M_SQRT2)));
  }
}


typedef void (*gill83fn_type)(double *fp, double *theta, int id);

void gill83fnF(double *fp, double *theta, int);
int gill83(double *hf, double *hphif, double *df, double *df2, double *ef,
	   double *theta, int cpar, double epsR, int K, double gillStep,
	   double fTol, int cid, gill83fn_type gill83fn);

gill83fn_type gill83fnG = &gill83fnF;


void updateEta(double *eta, int cid) {
  rx_solving_options_ind *ind =  &(rx->subjects[cid]);
  for (int i = op_focei.neta; i--;) {
    ind->par_ptr[op_focei.etaTrans[i]]=eta[i];
  }
}

static inline void gill83pred(double *fp, double *eta, int cid, int w) {
  // Calculate derivatives by finite difference
  updateEta(eta, cid);
  predOde(cid);
  focei_ind *fInd = &(inds_focei[cid]);
  rx_solving_options_ind *ind =  &(rx->subjects[cid]);
  rxPred.calc_lhs(cid, fInd->curT, fInd->curS, // Solve space is smaller
				ind->lhs);
  *fp = ind->lhs[w];
}

void gill83etaR(double *fp, double *eta, int cid) {
  gill83pred(fp, eta, cid, 1);
}

void gill83etaF(double *fp, double *eta, int cid) {
  gill83pred(fp, eta, cid, 0);
}

static inline double calcGradForEtaGeneral(double *eta,
					   double *aEps,
					   int cpar, int cid, int w) {
  // First try a gill difference
  if (op_focei.eventType != 1) aEps[cpar] = op_focei.eventFD;  
  if (op_focei.eventType == 1 && aEps[cpar] == 0) {
    double hf, hphif, err, gillDf, gillDf2, gillErr;
    int gillRet;
    if (w == 0) {
      gillRet = gill83(&hf, &hphif, &gillDf, &gillDf2, &gillErr,
		       eta, cpar, op_focei.gillRtol, op_focei.gillK, op_focei.gillStep, op_focei.gillFtol,
		       cid, gill83etaF);
    } else {
      gillRet = gill83(&hf, &hphif, &gillDf, &gillDf2, &gillErr,
		       eta, cpar, op_focei.gillRtol, op_focei.gillK, op_focei.gillStep, op_focei.gillFtol,
		       cid, gill83etaR);
    }
    err = 1/(std::fabs(eta[cpar])+1);
    aEps[cpar]  = hf*err;
    return gillDf;
  } else {
    // now use forward
    double delta = (std::fabs(eta[cpar])*aEps[cpar] + aEps[cpar]);
    focei_ind *fInd = &(inds_focei[cid]);
    rx_solving_options_ind *ind =  &(rx->subjects[cid]);
    double fv;
    ind->par_ptr[op_focei.etaTrans[cpar]]+=delta;
    predOde(cid); // Assumes same order of parameters
    rxPred.calc_lhs(cid, fInd->curT, fInd->curS, // Solve space is smaller
				ind->lhs);
    ind->par_ptr[op_focei.etaTrans[cpar]] -= delta;
    double hv = ind->lhs[w];
    double df = (hv - fInd->curF)/delta;
    if (op_focei.eventType == 2 || fabs(df) < op_focei.gradCalcCentralSmall){
      // Switch to central in close to zero derivatives
      ind->par_ptr[op_focei.etaTrans[cpar]] -= delta;
      predOde(cid);
      rxPred.calc_lhs(cid, fInd->curT, fInd->curS, // Solve space is smaller
		      ind->lhs);
      ind->par_ptr[op_focei.etaTrans[cpar]] += delta;
      df = (hv-ind->lhs[w])/(2*delta);
    }
    return df;
  }
  return 0.0;
}

static inline double calcGradForEtaF(double *eta,
		       double *aEps,
		       int cpar, int cid) {
  return calcGradForEtaGeneral(eta,aEps, cpar, cid, 0);
}

static inline double calcGradForEtaR(double *eta,
		       double *aEps,
		       int cpar, int cid) {
  return calcGradForEtaGeneral(eta,aEps, cpar, cid, 1);
}

double likInner0(double *eta, int id){
  rx = getRx();
  rx_solving_options_ind *ind = &(rx->subjects[id]);
  rx_solving_options *op = rx->op;
  int i, j;
  bool recalc = false;
  focei_ind *fInd= &(inds_focei[id]);
  if (op_focei.neta > 0){
    if (!fInd->setup){
      recalc = true;
      fInd->setup = 1;
    } else {
      // Check to see if old ETA matches.
      for (j = op_focei.neta; j--;){
	if (fInd->oldEta[j] != eta[j]){
	  recalc = true;
	  break;
	}
      }
    }
  } else {
    recalc = true;
  }
  if (recalc){
    for (j = op_focei.neta; j--;){
      ind->par_ptr[op_focei.etaTrans[j]] = eta[j];
    }
    if (op_focei.stickyRecalcN2 <= op_focei.stickyRecalcN){
      op_focei.stickyRecalcN2=0;
    }
    ind->solved = -1;
    // Solve ODE
    bool predSolve = false;
    if (fInd->doFD == 0) {
      innerOde(id);
      j=0;
      while (op_focei.stickyRecalcN2 <= op_focei.stickyRecalcN && op->badSolve && j < op_focei.maxOdeRecalc){
	op_focei.stickyRecalcN2++;
	op_focei.reducedTol=1;
	op_focei.reducedTol2=1;
	// Not thread safe
	RxODE::atolRtolFactor_(op_focei.odeRecalcFactor);
	ind->solved = -1;
	innerOde(id);
	j++;
      }
      if (j != 0) {
	if (op_focei.stickyRecalcN2 <= op_focei.stickyRecalcN){
	  // Not thread safe
	  RxODE::atolRtolFactor_(pow(op_focei.odeRecalcFactor, -j));
	} else {
	  op_focei.stickyTol=1;
	}
      }
    } else {
      predOde(id);
      predSolve=true;
      op_focei.didPredSolve = true;
    }
    if (op->neq > 0 && (ISNA(ind->solve[0]) || std::isnan(ind->solve[0]) ||
			std::isinf(ind->solve[0]))){
      return NA_REAL;
      //throw std::runtime_error("bad solve");
    } else {
      // Update eta.
      arma::mat lp(fInd->lp, op_focei.neta, 1, false, true);
      lp.zeros();
      arma::mat a(fInd->a, ind->n_all_times - ind->ndoses - ind->nevid2, op_focei.neta, false, true);
      arma::mat B(fInd->B, ind->n_all_times - ind->ndoses - ind->nevid2, 1, false, true);
      arma::mat c(fInd->c, ind->n_all_times - ind->ndoses - ind->nevid2,
		  op_focei.neta, false, true);
      arma::mat Vid(fInd->Vid, ind->n_all_times - ind->ndoses - ind->nevid2,
		    ind->n_all_times - ind->ndoses - ind->nevid2, false, true);
      if (op_focei.fo == 1){
	Vid.zeros();
      }
      // RSprintf("ID: %d; Solve #2: %f\n", id, ind->solve[2]);
      // Calculate matricies
      int k = 0, kk=0;//ind->n_all_times - ind->ndoses - ind->nevid2 - 1;
      fInd->llik=0.0;
      fInd->tbsLik=0.0;
      double f, err, r, fpm, fpm2, rp = 0,lnr, limit, dv,dv0, curT;
      int cens = 0;
      int oldNeq = op->neq;
      iniSubjectI(id, 1, ind, op, rx, rxInner.update_inis);
      for (j = 0; j < ind->n_all_times; ++j){
	ind->idx=j;
	kk = ind->ix[j];
	curT = getTimeF(kk, ind);
	dv0 = ind->dv[kk];
	if (isDose(ind->evid[kk])) {
	  // ind->tlast = ind->all_times[ind->ix[ind->idx]];
	  // Need to calculate for advan sensitivities
	  if (predSolve) {
	    rxPred.calc_lhs(id, curT, getSolve(j), ind->lhs);
	    ind->lhs[op_focei.neta + 1] = ind->lhs[1];
	  }
	  else {
	    rxInner.calc_lhs(id, curT, getSolve(j), ind->lhs);
	  }
	} else if (ind->evid[kk] == 0) {
	  if (predSolve) {
	    rxPred.calc_lhs(id, curT, getSolve(j), ind->lhs);
	    ind->lhs[op_focei.neta + 1] = ind->lhs[1];
	  } else {
	    rxInner.calc_lhs(id, curT, getSolve(j), ind->lhs);
	  }
	  f = ind->lhs[0]; // TBS is performed in the RxODE rx_pred_ statement. This allows derivatives of TBS to be propagated
	  dv = tbs(dv0);
	  if (ISNA(f) || std::isnan(f) || std::isinf(f)) {
	    return NA_REAL;
	    //throw std::runtime_error("bad solve");
	  }
	  // fInd->f(k, 0) = ind->lhs[0];
	  err = f - dv;
	  limit = R_NegInf;
	  if (rx->limit) {
	    limit = ind->limit[kk];
	    if (ISNA(limit)) {
	      limit = R_NegInf;
	    } else if (R_FINITE(limit)) {
	      limit = tbs(limit);
	    }
	  }
	  cens = 0;
	  if (rx->cens) cens = ind->cens[kk];
	  fInd->tbsLik+=tbsL(dv0);
	  // fInd->err(k, 0) = ind->lhs[0] - ind->dv[k]; // pred-dv
	  if (ISNA(ind->lhs[op_focei.neta + 1])){
	    return NA_REAL;
	    //throw std::runtime_error("bad solve");
	  }
	  r = ind->lhs[op_focei.neta + 1];
	  if (r == 0.0) {
	    r = 1.0;
	  }
	  if (op_focei.neta == 0) {
	    lnr =_safe_log(r);
	    //llik <- -0.5 * sum(err ^ 2 / R + log(R));
	    if (cens == 0){
	      fInd->llik += err * err/_safe_zero(r) + lnr;
	      likM2(fInd, limit, f, r);
	    } else if (cens != 0) {
	      likCens(fInd, cens, limit, f, dv, r);
	    }
	  } else if (op_focei.fo == 1) {
	    // FO
	    B(k, 0) = err; // res
	    Vid(k, k) = r;
	    for (i = op_focei.neta; i--; ) {
	      if (op_focei.etaFD[i]==0){
		a(k, i) = ind->lhs[i+1];
	      }
	    }
	    op->neq = op_focei.predNeq;
	    fInd->curT = curT;
	    fInd->curF = f;
	    fInd->curS = getSolve(j);
	    for (i = op_focei.neta; i--; ) {
	      if (predSolve || op_focei.etaFD[i]==1){
		// Calculate derivatives by finite difference
		a(k, i) = fpm = calcGradForEtaF(eta, fInd->etahf, i, id);
	      }
	    }
	    op->neq = oldNeq;
	    // Ci = fpm %*% omega %*% t(fpm) + Vi; Vi=diag(r)
	  } else {
	    lnr =_safe_log(ind->lhs[op_focei.neta + 1]);
	    // fInd->r(k, 0) = ind->lhs[op_focei.neta+1];
	    // B(k, 0) = 2.0/ind->lhs[op_focei.neta+1];
	    // lhs 0 = F
	    // lhs 1-eta = df/deta
	    // FIXME faster initialization via copy or elm
	    // RSprintf("id: %d k: %d j: %d\n", id, k, j);
	    B(k, 0) = 2.0/_safe_zero(r);
	    if (op_focei.interaction == 1) {
	      for (i = op_focei.neta; i--; ) {
		if (predSolve || op_focei.etaFD[i]==0) {
		  fpm = a(k, i) = ind->lhs[i + 1]; // Almquist uses different a (see eq #15)
		  rp  = ind->lhs[i + op_focei.neta + 2];
		  c(k, i) = rp/_safe_zero(r);
		}
	      }
	      // Cannot combine for loop with for loop above because
	      // calc_lhs overwrites the lhs memory.
	      op->neq = op_focei.predNeq;
	      op->neq = op_focei.predNeq;
	      fInd->curT = curT;
	      fInd->curS = getSolve(j);
	      for (i = op_focei.neta; i--; ) {
		// Calculate finite difference derivatives if needed
		if (predSolve || op_focei.etaFD[i]==1){
		  // Calculate derivatives by finite difference
		  fInd->curF = f;
		  a(k, i) = fpm = calcGradForEtaF(eta, fInd->etahf, i, id);
		  fInd->curF = r;
		  rp = calcGradForEtaR(eta, fInd->etahr, i, id);
		  if (fpm == 0.0) {
		    a(k, i) = fpm = sqrt(DBL_EPSILON);
		  }
		  if (rp == 0.0) {
		    rp = sqrt(DBL_EPSILON);
		  }
		  c(k, i) = rp/_safe_zero(r);
		} else {
		  fpm = a(k, i);
		  rp  = ind->lhs[i + op_focei.neta + 2];
		  rp = c(k, i)*r;
		}
		// This is calculated at the end; That way it is
		// correct for finite difference and sensitivity.

		//lp is eq 12 in Almquist 2015
		// .5*apply(eps*fp*B + .5*eps^2*B*c - c, 2, sum) - OMGAinv %*% ETA
		if (cens == 0) {
		  // REprintf("t: %f: err: %f; fpm: %f; B(k, 0): %f; c(k, i): %f; rp: %f", ind->all_times[kk], err, fpm, B(k, 0), c(k, i), rp);
		  lp(i, 0)  += 0.25 * err * err * B(k, 0) * c(k, i) -
		    0.5 * c(k, i) - 0.5 * err * fpm * B(k, 0);
		  // REprintf("lp(i,0): %f\n", lp(i,0));
		} else {
		  if (R_FINITE(limit)){
		    // M3 method
		    // logLik = log(phi((QL-f(x)/sqrt(g(x)))))
		    // logLik = log(1/2*(1+erf((cens*(dv-f(x)))/sqrt(g(x))))/sqrt(2))
		    // > D(S("log(1/2*(1+erf((cens*(dv-f(x)))/sqrt(g(x))/M_SQRT2)))"),"x")
		    // (Mul) 2*exp(-(dv - f(x))^2*cens^2/(g(x)*M_SQRT2^2))*(-Derivative(f(x), x)*cens/(sqrt(g(x))*M_SQRT2) + (-1/2)*Derivative(g(x), x)*(dv - f(x))*cens/(g(x)^(3/2)*M_SQRT2))/(sqrt(pi)*(1 + erf((dv - f(x))*cens/(sqrt(g(x))*M_SQRT2))))
		    // 2*exp(-(dv - f)^2*cens^2/(r*M_SQRT2^2))*(-fpm*cens/(sqrt(r)*M_SQRT2) + (-0.5)*rp*(dv - f)*cens/(r^(1.5)*M_SQRT2))/(sqrt(pi)*(1 + erf((dv - f)*cens/(sqrt(r)*M_SQRT2))))
		    double rx_expr_0=dv-f;
		    double rx_expr_1=_safe_sqrt(r);
		    double rx_expr_2=rx_expr_1*M_SQRT2;
		    lp(i, 0) += 2*exp(-(rx_expr_0*rx_expr_0)/_safe_zero(r*2))*
		      (-fpm*(double)(cens)/_safe_zero(rx_expr_2)-0.5*rp*rx_expr_0*(double)(cens)/_safe_zero(R_pow(r,1.5)*M_SQRT2))/(M_SQRT_PI*(1+erf(rx_expr_0*cens/_safe_zero(rx_expr_2))));
		  } else {
		    // M4 method
		    // dv = lloq
		    // logLik = log(phi((dv-f(x))/sqrt(g(x))))+log(phi((limit-f(x))/sqrt(g(x))))-log(1-phi((limit-f(x))/sqrt(g(x))))
		    // > D(S("log(0.5*(1+erf( (dv-f(x))/sqrt(g(x)) /M_SQRT2)))+log(0.5*(1+erf( (limit-f(x))/sqrt(g(x)) /M_SQRT2)))-log(1-0.5*(1+erf( (limit-f(x))/sqrt(g(x)) /M_SQRT2)))"),"x")
		    // (Add) 2.0*exp(-(dv - f(x))^2/(g(x)*M_SQRT2^2))*(-Derivative(f(x), x)/(sqrt(g(x))*M_SQRT2) + (-1/2)*Derivative(g(x), x)*(dv - f(x))/(g(x)^(3/2)*M_SQRT2))/(sqrt(pi)*(1 + erf((dv - f(x))/(sqrt(g(x))*M_SQRT2)))) + 1.0*exp(-(limit - f(x))^2/(g(x)*M_SQRT2^2))*(-Derivative(f(x), x)/(sqrt(g(x))*M_SQRT2) + (-1/2)*Derivative(g(x), x)*(limit - f(x))/(g(x)^(3/2)*M_SQRT2))/(sqrt(pi)*(1 - 0.5*(1 + erf((limit - f(x))/(sqrt(g(x))*M_SQRT2))))) + 2.0*exp(-(limit - f(x))^2/(g(x)*M_SQRT2^2))*(-Derivative(f(x), x)/(sqrt(g(x))*M_SQRT2) + (-1/2)*Derivative(g(x), x)*(limit - f(x))/(g(x)^(3/2)*M_SQRT2))/(sqrt(pi)*(1 + erf((limit - f(x))/(sqrt(g(x))*M_SQRT2))))
		    // 2.0*exp(-(dv - f)^2/(r*M_SQRT2^2))*(-fpm/(sqrt(r)*M_SQRT2) + (-0.5)*rp*(dv - f)/(r^(1.5)*M_SQRT2))/(M_SQRT_PI*(1 + erf((dv - f)/(sqrt(r)*M_SQRT2)))) + 1.0*exp(-(limit - f)^2/(r*M_SQRT2^2))*(-fpm/(sqrt(r)*M_SQRT2) + (-0.5)*rp*(limit - f)/(r^(1.5)*M_SQRT2))/(M_SQRT_PI*(1 - 0.5*(1 + erf((limit - f)/(sqrt(r)*M_SQRT2))))) + 2.0*exp(-(limit - f)^2/(r*M_SQRT2^2))*(-fpm/(sqrt(r)*M_SQRT2) + (-0.5)*rp*(limit - f)/(r^(1.5)*M_SQRT2))/(M_SQRT_PI*(1 + erf((limit - f)/(sqrt(r)*M_SQRT2))))
		    // 2.0*exp(-(dv - f)^2/(r*2))*(-fpm/(sqrt(r)*M_SQRT2) + (-0.5)*rp*(dv - f)/(r^(1.5)*M_SQRT2))/(M_SQRT_PI*(1 + erf((dv - f)/(sqrt(r)*M_SQRT2)))) + 1.0*exp(-(limit - f)^2/(r*2))*(-fpm/(sqrt(r)*M_SQRT2) + (-0.5)*rp*(limit - f)/(r^(1.5)*M_SQRT2))/(M_SQRT_PI*(1 - 0.5*(1 + erf((limit - f)/(sqrt(r)*M_SQRT2))))) + 2.0*exp(-(limit - f)^2/(r*2))*(-fpm/(sqrt(r)*M_SQRT2) + (-0.5)*rp*(limit - f)/(r^(1.5)*M_SQRT2))/(M_SQRT_PI*(1 + erf((limit - f)/(sqrt(r)*M_SQRT2))))
		    double rx_expr_0=r*2;
		    double rx_expr_1=dv-f;
		    double rx_expr_2=_safe_sqrt(r);
		    double rx_expr_3=R_pow(r,1.5);
		    double rx_expr_4=limit-f;
		    double rx_expr_5=-0.5*rp;
		    double rx_expr_6=rx_expr_4*rx_expr_4;
		    double rx_expr_7=rx_expr_2*M_SQRT2;
		    double rx_expr_8=rx_expr_3*M_SQRT2;
		    double rx_expr_9=rx_expr_5*rx_expr_4;
		    double rx_expr_10=exp(-rx_expr_6/_safe_zero(rx_expr_0));
		    double rx_expr_11=(rx_expr_4)/_safe_zero(rx_expr_7);
		    double rx_expr_12=erf(rx_expr_11);
		    double rx_expr_13=1+rx_expr_12;
		    double rx_expr_14=rx_expr_9/_safe_zero(rx_expr_8);
		    lp(i, 0)=2*exp(-(rx_expr_1*rx_expr_1)/_safe_zero(rx_expr_0))*(-fpm/_safe_zero(rx_expr_7)+rx_expr_5*(rx_expr_1)/_safe_zero(rx_expr_8))/(M_SQRT_PI*(1+erf((rx_expr_1)/_safe_zero(rx_expr_7))))+1*rx_expr_10*(-fpm/_safe_zero(rx_expr_7)+rx_expr_14)/_safe_zero(M_SQRT_PI*(1-0.5*(rx_expr_13)))+2*rx_expr_10*(-fpm/_safe_zero(rx_expr_7)+rx_expr_14)/_safe_zero(M_SQRT_PI*(rx_expr_13));
		  }
		}
	      }
	      op->neq = oldNeq;
	      // Eq #10
	      //llik <- -0.5 * sum(err ^ 2 / R + log(R));
	      if (cens == 0){
		fInd->llik += err * err/_safe_zero(r) + lnr;
		likM2(fInd, limit, f, r);
	      } else if (cens != 0) {
		likCens(fInd, cens, limit, f, dv, r);
	      }
	    } else if (op_focei.interaction == 0){
	      for (i = op_focei.neta; i--; ){
		if (op_focei.etaFD[i]==0){
		  a(k, i) = ind->lhs[i + 1];
		}
	      }
	      op->neq = op_focei.predNeq;
	      fInd->curT = curT;
	      fInd->curF = f;
	      fInd->curS = getSolve(j);
	      for (i = op_focei.neta; i--; ) {
		// Calculate derivatives by finite difference
		if (predSolve || op_focei.etaFD[i]==1) {
		  a(k, i) = fpm = calcGradForEtaF(eta, fInd->etahf, i, id);
		} else {
		  fpm = a(k, i);
		}
		if (cens == 0){
		  lp(i, 0) -= 0.5 * err * fpm * B(k, 0);
		} else {
		  if (std::isinf(limit)){
		    // M3 method
		    // logLik = log(phi((QL-f(x)/sqrt(g(x)))))
		    // logLik = log(1/2*(1+erf((cens*(dv-f(x)))/sqrt(g(x))))/sqrt(2))
		    // > D(S("log(1/2*(1+erf((cens*(dv-f(x)))/sqrt(g)/M_SQRT2)))"),"x")
		    // -2*exp(-0.5*(dv - f)^2/(g))*fpm*cens/(sqrt(g)*M_SQRT_PI*M_SQRT2*(1 + erf((dv - f)*cens/(sqrt(g)*M_SQRT2))))
		    double rx_expr_0=dv-f;
		    double rx_expr_1=_safe_sqrt(r);
		    lp(i, 0) -= 2*exp(-0.5*rx_expr_0*rx_expr_0/_safe_zero(r))*fpm*(double)(cens)/
		      _safe_zero(rx_expr_1*M_SQRT_PI*M_SQRT2*(1+erf((rx_expr_0)*cens/(rx_expr_1*M_SQRT2))));
		  } else {
		    // M4
		    // > D(S("log(0.5*(1+erf((dv-f(x))/sqrt(r) /M_SQRT2)))+log(0.5*(1+erf( (limit-f(x))/sqrt(r)/M_SQRT2)))-log(1-0.5*(1+erf( (limit-f(x))/sqrt(r)/M_SQRT2)))"),"x")
		    // -2.0*exp(-0.5*(dv - f)^2/r)*fpm/(sqrt(r)*M_SQRT_PI*M_SQRT2*(1 + erf((dv - f)/(sqrt(r)*M_SQRT2)))) - 1.0*exp(-(limit - f)^2/(r*2))*fpm/(sqrt(r)*M_SQRT_PI*M_SQRT2*(1 - 0.5*(1 + erf((limit - f)/(sqrt(r)*M_SQRT2))))) - 2.0*exp(-(limit - f)^2/(r*2))*fpm/(sqrt(r)*M_SQRT_PI*M_SQRT2*(1 + erf((limit - f)/(sqrt(r)*M_SQRT2))))
		    double rx_expr_0=r*2;
		    double rx_expr_1=dv-f;
		    double rx_expr_2=_safe_sqrt(r);
		    double rx_expr_3=limit-f;
		    double rx_expr_4=rx_expr_3*rx_expr_3;
		    double rx_expr_5=rx_expr_2*M_SQRT2;
		    double rx_expr_6=rx_expr_2*M_SQRT_PI;
		    double rx_expr_7=exp(-rx_expr_4/rx_expr_0);
		    double rx_expr_8=rx_expr_6*M_SQRT2;
		    double rx_expr_9=rx_expr_3/rx_expr_5;
		    double rx_expr_10=erf(rx_expr_9);
		    double rx_expr_11=1+rx_expr_10;
		    lp(i, 0)-=2*exp(-0.5*(rx_expr_1*rx_expr_1)/_safe_zero(r))*fpm/(rx_expr_8*(1+erf(rx_expr_1/(rx_expr_5))))-1*rx_expr_7*fpm/(rx_expr_8*(1-0.5*(rx_expr_11)))-2*rx_expr_7*fpm/(rx_expr_8*(rx_expr_11));
		  }
		}
	      }
	      op->neq = oldNeq;
	      // Eq #10
	      //llik <- -0.5 * sum(err ^ 2 / R + log(R));
	      if (cens == 0){
		fInd->llik += err * err/_safe_zero(r) + lnr;
		likM2(fInd, limit, f, r);
	      } else {
		likCens(fInd, cens, limit, f, dv, r);
	      }
	    }
	  }
	  // k--;
	  k++;
	}
      }
      if (op_focei.neta == 0) {
	fInd->llik = -0.5*fInd->llik;
      } else if (op_focei.fo == 1) {
	if (cens != 0) stop("FO censoring not supported.");
	mat Ci = a * op_focei.omega * trans(a) + Vid;
	mat cholCi = cholSE__(Ci, op_focei.cholSEtol);
	mat CiInv;
	bool success  = inv(CiInv, trimatu(cholCi));
	if (!success){
	  CiInv = pinv(trimatu(cholCi));
	}
	CiInv = CiInv * CiInv.t();
	double lik = 0;
	// 2*sum(log(diag(chol(Ci))))
	for (unsigned int j = cholCi.n_rows; j--;){
	  lik += 2*_safe_log(cholCi(j,j));
	}
	// + t(.$Ri) %*% solve(Ci) %*% .$Ri
	mat rest =trans(B) * CiInv * B;
	lik += rest(0,0);
	// lik = -2*ll
	fInd->llik = -0.5*lik;
	// if (cens == 0){
	// 	fInd->llik += err * err/_safe_zero(r) + lnr;
	// 	likM2(fInd, limit, f, r);
	//       } else {
	// 	likCens(fInd, cens, limit, f, dv, r);
	//       }
      } else {
	fInd->llik = -0.5*fInd->llik;
	// Now finalize lp
	mat etam = arma::mat(op_focei.neta, 1);
	std::copy(&eta[0], &eta[0] + op_focei.neta, etam.begin()); // fill in etam
	// print(wrap(lp));
	// print(wrap(etam));
	// Finalize eq. #12
	lp = -(lp - op_focei.omegaInv * etam);
	// Partially finalize #10
	fInd->llik = -trace(fInd->llik - 0.5*(etam.t() * op_focei.omegaInv * etam));
	// print(wrap(fInd->llik));
	std::copy(&eta[0], &eta[0] + op_focei.neta, &fInd->oldEta[0]);
	// for (int ssi = op_focei.neta; ssi--;){
	//   // RSprintf("ssi: %d :%d;\n",id, ssi);
	//   // RSprintf("eta: %f\n", eta[ssi]);
	//   fInd->oldEta[ssi] = eta[ssi];
	// }
      }
    }
  }
  return fInd->llik;
}

double *lpInner(double *eta, double *g, int id){
  focei_ind *fInd = &(inds_focei[id]);
  likInner0(eta, id);
  std::copy(&fInd->lp[0], &fInd->lp[0] + op_focei.neta,
	    &g[0]);
  return &g[0];
}

//[[Rcpp::export]]
NumericVector foceiInnerLp(NumericVector eta, int id = 1){
  double *etad = new double[eta.size()];
  std::copy(eta.begin(),eta.end(),&etad[0]);
  NumericVector lp(eta.size());
  int curId = id - 1;
  lpInner(etad,&lp[0], curId);
  delete[] etad;
  return lp;
}

//[[Rcpp::export]]
double likInner(NumericVector eta, int id = 1){
  double *etad = new double[eta.size()];
  std::copy(eta.begin(),eta.end(),&etad[0]);
  int curId = id - 1;
  double llik = likInner0(etad, curId);
  delete[] etad;
  return llik;
}


double LikInner2(double *eta, int likId, int id){
  focei_ind *fInd = &(inds_focei[id]);
  // print(wrap(-likInner0(eta)));
  // print(wrap(op_focei.logDetOmegaInv5));
  double lik=0;
  if (op_focei.neta == 0) {
    lik = fInd->llik;
  } else if (op_focei.fo == 1){
    // Already almost completely calculated.
    lik = fInd->llik;
  } else {
    // print(wrap(op_focei.logDetOmegaInv5));
    lik = -likInner0(eta, id) + op_focei.logDetOmegaInv5;
    // print(wrap(lik));
    rx = getRx();
    rx_solving_options_ind *ind = &(rx->subjects[id]);
    rx_solving_options *op = rx->op;
    if (op->neq > 0 && ISNA(ind->solve[0])){
      //return 1e300;
      return NA_REAL;
    }
    // Calculate lik first to calculate components for Hessian
    // Hessian
    mat H(fInd->H, op_focei.neta, op_focei.neta, false, true);
    H.zeros();
    int k, l;
    mat tmp;

    arma::mat a(fInd->a, ind->n_all_times - ind->ndoses - ind->nevid2, op_focei.neta, false, true);
    // std::copy(&fInd->a[0], &fInd->a[0]+a.size(), a.begin());
    arma::mat B(fInd->B, ind->n_all_times - ind->ndoses - ind->nevid2, 1, false, true);
    // std::copy(&fInd->B[0], &fInd->B[0]+B.size(), B.begin());

    // This is actually -H
    if (op_focei.interaction){
      arma::mat c(fInd->c, ind->n_all_times - ind->ndoses - ind->nevid2, op_focei.neta, false, true);
      // std::copy(&fInd->c[0], &fInd->c[0]+c.size(), c.begin());
      for (k = op_focei.neta; k--;){
        for (l = k+1; l--;){
          // tmp = fInd->a.col(l) %  fInd->B % fInd->a.col(k);
          H(k, l) = 0.5*sum(a.col(l) % B % a.col(k) +
                            c.col(l) % c.col(k)) +
            op_focei.omegaInv(k, l);
          H(l, k) = H(k, l);
        }
      }
    } else {
      for (k = op_focei.neta; k--;){
        for (l = k+1; l--;){
          // tmp = a.col(l) %  B % a.col(k);
          H(k, l) = 0.5*sum(a.col(l) % B % a.col(k)) +
            op_focei.omegaInv(k, l);
          H(l, k) = H(k, l);
        }
      }
    }
    arma::mat H0(fInd->H0, op_focei.neta, op_focei.neta, false, true);
    k=0;
    if (fInd->doChol){
      if (!chol(H0, H)) {
	return NA_REAL;
      }
    } else {
      H0=cholSE__(H, op_focei.cholSEtol);
    }
    // - sum(log(H.diag()));
    for (unsigned int j = H0.n_rows; j--;){
      lik -= _safe_log(H0(j,j));
    }
  }
  lik += fInd->tbsLik;
  if (likId == 0){
    fInd->lik[0] = lik;
    if (op_focei.neta != 0) std::copy(&fInd->eta[0], &fInd->eta[0] + op_focei.neta, &fInd->saveEta[0]);
  } else {
    // Use Objective function for Numeric Gradient.
    fInd->lik[likId] = -2*lik;
  }
  return lik;
}

extern "C" double innerOptimF(int n, double *x, void *ex){
  int *id = (int*)ex;
  focei_ind *fInd = &(inds_focei[*id]);
  if (fInd->badSolve == 1) return NA_REAL;
  double f = likInner0(x, *id);
  if (ISNA(f)) {
    fInd->badSolve = 1;
  }
  fInd->nInnerF++;
  return f;
}

extern "C" void innerOptimG(int n, double *x, double *g, void *ex) {
  int *id = (int*)ex;
  focei_ind *fInd = &(inds_focei[*id]);
  if (fInd->badSolve == 1) return;
  lpInner(x, g, *id);
  fInd->nInnerG++;
}

// Scli-lab style cost function for inner
void innerCost(int *ind, int *n, double *x, double *f, double *g, int *ti, float *tr, double *td, int *id){
  rx = getRx();
  // if (*id < 0 || *id >= rx->nsub){
  //   // Stops from accessing bad memory, but it doesn't fix any
  //   // problems here.  Rather, this allows the error without a R
  //   // session crash.
  //   stop("Unexpected id for solving (id=%d and should be between 0 and %d)", *id, rx->nsub);
  // }
  focei_ind *fInd = &(inds_focei[*id]);
  if (fInd->badSolve==1) {
      return;
  }
  if (*ind==2 || *ind==4) {
    // Function
    // Make sure ID remains installed
    *f = likInner0(x, *id);
    if (ISNA(*f))  {
      fInd->badSolve=1;
    }
    fInd->nInnerF++;
    // if (op_focei.printInner != 0 && fInd->nInnerF % op_focei.printInner == 0){
    // RSprintf("%03d: ", *id);
    // for (int i = 0; i < *n; i++) RSprintf(" %#10g", x[i]);
    // RSprintf(" (nG: %d)\n", fInd->nInnerG);
    // }
  }
  if (*ind==3 || *ind==4) {
    // Gradient

    lpInner(x, g, *id);
    fInd->nInnerG++;
  }
}

static inline int innerEval(int id){
  focei_ind *fInd = &(inds_focei[id]);
  // Use eta
  double lik0 = likInner0(fInd->eta, id);
  if (ISNA(lik0)) return 0;
  lik0 = LikInner2(fInd->eta, 0, id);
  if (ISNA(lik0)) return 0;
  return 1;
}

static inline int innerOpt1(int id, int likId) {
  focei_ind *fInd = &(inds_focei[id]);
  focei_options *fop = &op_focei;
  if (op_focei.neta == 0) {
    double lik = likInner0(NULL, id);
    if (ISNA(lik)) return 0;
    lik = LikInner2(NULL, 0, id);
    if (ISNA(lik)) return 0;
    return 1;
  }
  fInd->nInnerF=0;
  fInd->nInnerG=0;
  bool n1qn1Inner = true;
  // Use eta
  // Convert Zm to Hessian, if applicable.
  mat etaMat(fop->neta, 1);
  if (!op_focei.calcGrad){
    if (op_focei.resetEtaSize <= 0){
      if (op_focei.resetHessianAndEta){
	fInd->mode = 1;
	fInd->uzm = 1;
	if (n1qn1Inner) op_focei.didHessianReset=1;
      }
      std::fill(&fInd->eta[0], &fInd->eta[0] + op_focei.neta, 0.0);
      op_focei.didEtaReset=1;
    } else if (R_FINITE(op_focei.resetEtaSize)) {
      std::copy(&fInd->eta[0], &fInd->eta[0] + op_focei.neta, etaMat.begin());
      // Standardized ETAs
      // chol(omega^-1) %*% eta
      mat etaRes = op_focei.cholOmegaInv * etaMat;
      bool doBreak = false;
      for (unsigned int j = etaRes.n_rows; j--;){
	if (std::fabs(etaRes(j, 0)) >= op_focei.resetEtaSize){
	  if (op_focei.resetHessianAndEta){
	    fInd->mode = 1;
	    fInd->uzm = 1;
	    if (n1qn1Inner) op_focei.didHessianReset=1;
	  }
	  std::fill(&fInd->eta[0], &fInd->eta[0] + op_focei.neta, 0.0);
	  op_focei.didEtaReset=1;
	  doBreak=true;
	  break;
	}
      }
      if (!doBreak){
	etaRes = op_focei.eta1SD % etaMat;
	for (unsigned int j = etaRes.n_rows; j--;){
	  if (std::fabs(etaRes(j, 0)) >= op_focei.resetEtaSize){
	    if (op_focei.resetHessianAndEta){
	      fInd->mode = 1;
	      fInd->uzm = 1;
	      if (n1qn1Inner) op_focei.didHessianReset=1;
	    }
	    std::fill(&fInd->eta[0], &fInd->eta[0] + op_focei.neta, 0.0);
	    op_focei.didEtaReset=1;
	    break;
	  }
	}
      }
    }
  }
  if (n1qn1Inner) {
    updateZm(fInd);
    std::fill_n(&fInd->var[0], fop->neta, 0.1);
  }
  int npar = fop->neta;
  // if (std::isnan(fInd->eta[0])) {
  //   std::fill(&fInd->eta[0], &fInd->eta[0]+fop->neta, 0.0);
  //   op_focei.didEtaReset=1;
  // }
  std::copy(&fInd->eta[0], &fInd->eta[0]+fop->neta, fInd->x);
  double f, epsilon = max2(fop->epsilon, sqrt(DBL_EPSILON));

  // Since these are pointers, without reassignment they are modified.
  int mode = fInd->mode, maxInnerIterations=fop->maxInnerIterations,
    nsim=fop->nsim, imp=fop->imp;
  int izs; float rzs; double dzs;

  int nF = fInd->nInnerF;
  if (n1qn1Inner) {
    fInd->badSolve = 0;
    n1qn1_(innerCost, &npar, fInd->x, &f, fInd->g,
	   fInd->var, &epsilon,
	   &mode, &maxInnerIterations, &nsim,
	   &imp, fInd->zm, &izs, &rzs, &dzs, &id);
    if (ISNA(f)) return 0;
    nF = fInd->nInnerF-nF;
    // REprintf("innerCost id: %d, fInd->nInnerF: %d", id, fInd->nInnerF);
    // If stays at zero try another point?
    if (fInd->doEtaNudge == 1 && op_focei.etaNudge != 0.0){
      bool tryAgain=false;
      // if (nF <= 3) tryAgain = true;
      op_focei.didEtaNudge =1;
      if (!tryAgain){
	tryAgain = true;
	for (int i = fop->neta; i--;){
	  if (fInd->x[i] != 0.0){
	    tryAgain=false;
	    break;
	  }
	}
      }
      if (tryAgain){
	fInd->mode = 1;
	fInd->uzm = 1;
	op_focei.didHessianReset=1;
	std::fill_n(fInd->x, fop->neta, op_focei.etaNudge);
	//nF = fInd->nInnerF;
	fInd->badSolve = 0;
	n1qn1_(innerCost, &npar, fInd->x, &f, fInd->g,
	       fInd->var, &epsilon,
	       &mode, &maxInnerIterations, &nsim,
	       &imp, fInd->zm,
	       &izs, &rzs, &dzs, &id);
	if (ISNA(f)) return 0;
	// nF = fInd->nInnerF - nF;
	// if (nF > 3) tryAgain = false;
	if (!tryAgain) {
	  tryAgain = true;
	  for (int i = fop->neta; i--;){
	    if (fInd->x[i] != op_focei.etaNudge){
	      tryAgain=false;
	      break;
	    }
	  }
	}
	if (tryAgain) {
	  fInd->mode = 1;
	  fInd->uzm = 1;
	  op_focei.didHessianReset=1;
	  std::fill_n(fInd->x, fop->neta, -op_focei.etaNudge);
	  nF = fInd->nInnerF;
	  fInd->badSolve = 0;
	  n1qn1_(innerCost, &npar, fInd->x, &f, fInd->g,
		 fInd->var, &epsilon,
		 &mode, &maxInnerIterations, &nsim,
		 &imp, fInd->zm, &izs, &rzs, &dzs, &id);
	  if (ISNA(f)) return 0;
	  // nF = fInd->nInnerF - nF;
	  // if (nF > 3) tryAgain = false;
	  if (!tryAgain){
	    tryAgain = true;
	    for (int i = fop->neta; i--;){
	      if (fInd->x[i] != -op_focei.etaNudge){
		tryAgain=false;
		break;
	      }
	    }
	  }
	  if (tryAgain){
	    fInd->mode = 1;
	    fInd->uzm = 1;
	    op_focei.didHessianReset=1;
	    std::fill_n(fInd->x, fop->neta, -op_focei.etaNudge2);
	    nF = fInd->nInnerF;
	    fInd->badSolve = 0;
	    n1qn1_(innerCost, &npar, fInd->x, &f, fInd->g,
		   fInd->var, &epsilon,
		   &mode, &maxInnerIterations, &nsim,
		   &imp, fInd->zm, &izs, &rzs, &dzs, &id);
	    if (ISNA(f)) return 0;
	    // nF = fInd->nInnerF - nF;
	    // if (nF > 3) tryAgain = false;
	    if (!tryAgain){
	      tryAgain = true;
	      for (int i = fop->neta; i--;){
		if (fInd->x[i] != -op_focei.etaNudge2){
		  tryAgain=false;
		  break;
		}
	      }
	    }
	    if (tryAgain) {
	      fInd->mode = 1;
	      fInd->uzm = 1;
	      op_focei.didHessianReset=1;
	      std::fill_n(fInd->x, fop->neta, +op_focei.etaNudge2);
	      nF = fInd->nInnerF;
	      fInd->badSolve = 0;
	      n1qn1_(innerCost, &npar, fInd->x, &f, fInd->g,
		     fInd->var, &epsilon,
		     &mode, &maxInnerIterations, &nsim,
		     &imp, fInd->zm, &izs, &rzs, &dzs, &id);
	      if (ISNA(f)) return 0;
	      // nF = fInd->nInnerF - nF;
	      // if (nF > 3) tryAgain = false;
	      if (!tryAgain){
		tryAgain = true;
		for (int i = fop->neta; i--;){
		  if (fInd->x[i] != +op_focei.etaNudge2){
		    tryAgain=false;
		    break;
		  }
		}
	      }
	      if (tryAgain) {
		std::fill_n(fInd->x, fop->neta, 0);
		std::fill_n(&fInd->var[0], fop->neta, 0.2);
		nF = fInd->nInnerF;
		fInd->badSolve = 0;
		n1qn1_(innerCost, &npar, fInd->x, &f, fInd->g,
		       fInd->var, &epsilon,
		       &mode, &maxInnerIterations, &nsim,
		       &imp, fInd->zm,
		       &izs, &rzs, &dzs, &id);
		if (ISNA(f)) return 0;
		//nF = fInd->nInnerF-nF;
		// if (nF > 3) tryAgain = false;
		if (!tryAgain){
		  tryAgain = true;
		  for (int i = fop->neta; i--;){
		    if (fInd->x[i] != 0.0){
		      tryAgain=false;
		      break;
		    }
		  }
		}
		if (tryAgain) {
		  std::fill_n(fInd->x, fop->neta, 0);
		}
		std::fill_n(&fInd->var[0], fop->neta, 0.1);
	      }
	    }

	  }
	}
      }
    }
  } else {
    int fail=0, fncount=0, grcount=0;
    char msg[100];
    fInd->badSolve = 0;
    lbfgsb3C(npar, op_focei.lmm, fInd->x, op_focei.etaLower,
	     op_focei.etaUpper, op_focei.nbdInner, &f, innerOptimF, innerOptimG,
	     &fail, (void*)(&id), op_focei.factr,
	     op_focei.pgtol, &fncount, &grcount,
	     op_focei.maxInnerIterations, msg, 0, -1,
	     op_focei.abstol, op_focei.reltol, fInd->g);
    if (ISNA(f)) return 0;
    // if (fail != 6 && fail != 7 && fail != 8 && fail != 27){
    //   // did not converge
    //   if (fInd->doEtaNudge == 1 && op_focei.etaNudge != 0.0){
    // 	std::fill_n(fInd->x, fop->neta, op_focei.etaNudge);
    // 	fail=0;
    // 	lbfgsb3C(npar, op_focei.lmm, fInd->x, op_focei.etaLower,
    // 		 op_focei.etaUpper, op_focei.nbdInner, &f, innerOptimF, innerOptimG,
    // 		 &fail, (void*)(&id), op_focei.factr,
    // 		 op_focei.pgtol, &fncount, &grcount,
    // 		 op_focei.maxInnerIterations, msg, 0, -1,
    // 		 op_focei.abstol, op_focei.reltol, fInd->g);
    // 	if (fail != 6 && fail != 7 && fail != 8 && fail != 27){
    // 	  std::fill_n(fInd->x, fop->neta, -op_focei.etaNudge);
    // 	  lbfgsb3C(npar, op_focei.lmm, fInd->x, op_focei.etaLower,
    // 		 op_focei.etaUpper, op_focei.nbdInner, &f, innerOptimF, innerOptimG,
    // 		 &fail, (void*)(&id), op_focei.factr,
    // 		 op_focei.pgtol, &fncount, &grcount,
    // 		 op_focei.maxInnerIterations, msg, 0, -1,
    // 		 op_focei.abstol, op_focei.reltol, fInd->g);
    // 	  if (fail != 6 && fail != 7 && fail != 8 && fail != 27){
    // 	    std::fill_n(fInd->x, fop->neta, 0);
    // 	  }
    // 	}
    //   }
    // }
  }
  // only nudge once
  fInd->doEtaNudge=0;

  std::copy(&fInd->x[0],&fInd->x[0]+fop->neta,&fInd->eta[0]);
  // Update variances
  std::copy(&fInd->eta[0], &fInd->eta[0] + op_focei.neta, etaMat.begin());
  op_focei.n = op_focei.n + 1.0;
  mat oldM = op_focei.etaM;
  op_focei.etaM = op_focei.etaM + (etaMat - op_focei.etaM)/op_focei.n;
  op_focei.etaS = op_focei.etaS + (etaMat - op_focei.etaM) %  (etaMat - oldM);
  fInd->llik = f;
  // Use saved Hessian on next opimization.
  fInd->mode=2;
  fInd->uzm =0;
  if (ISNA(LikInner2(fInd->eta, likId, id))) return 0;
  return 1;
}

void parHistData(Environment e, bool focei);

static inline void thetaReset00(NumericVector &thetaIni, NumericVector &omegaTheta, arma::mat &etaMat) {
  Function loadNamespace("loadNamespace", R_BaseNamespace);
  Environment nlmixr = loadNamespace("nlmixr");
  Environment thetaReset = nlmixr[".thetaReset"];
  focei_options *fop = &op_focei;
  thetaReset["maxInnerIterations"]=fop->maxInnerIterations;
  thetaReset["etaMat"] = wrap(etaMat);
  thetaReset["thetaIni"]= thetaIni;
  thetaReset["omegaTheta"] = omegaTheta;
  thetaReset["nF"] = op_focei.nF+op_focei.nF2;
  // Save gill info to skip recalc.
  IntegerVector gillRetC(op_focei.npars);
  std::copy(&op_focei.gillRetC[0], &op_focei.gillRetC[0]+op_focei.npars, gillRetC.begin());
  thetaReset["gillRetC"] = gillRetC;
  IntegerVector gillRet(op_focei.npars);
  std::copy(&op_focei.gillRet[0], &op_focei.gillRet[0]+op_focei.npars, gillRet.begin());
  thetaReset["gillRet"] = gillRet;
  NumericVector gillDf(op_focei.npars);
  std::copy(&op_focei.gillDf[0], &op_focei.gillDf[0]+op_focei.npars, gillDf.begin());
  thetaReset["gillDf"] = gillDf;
  NumericVector gillDf2(op_focei.npars);
  std::copy(&op_focei.gillDf2[0], &op_focei.gillDf2[0]+op_focei.npars, gillDf2.begin());
  thetaReset["gillDf2"] = gillDf2;
  NumericVector gillErr(op_focei.npars);
  std::copy(&op_focei.gillErr[0], &op_focei.gillErr[0]+op_focei.npars, gillErr.begin());
  thetaReset["gillErr"] = gillErr;
  NumericVector rEps(op_focei.npars);
  std::copy(&op_focei.rEps[0], &op_focei.rEps[0]+op_focei.npars, rEps.begin());
  thetaReset["rEps"] = rEps;
  NumericVector aEps(op_focei.npars);
  std::copy(&op_focei.aEps[0], &op_focei.aEps[0]+op_focei.npars, aEps.begin());
  thetaReset["aEps"] = aEps;
  NumericVector rEpsC(op_focei.npars);
  std::copy(&op_focei.rEpsC[0], &op_focei.rEpsC[0]+op_focei.npars, rEpsC.begin());
  thetaReset["rEpsC"] = rEpsC;
  NumericVector aEpsC(op_focei.npars);
  std::copy(&op_focei.aEpsC[0], &op_focei.aEpsC[0]+op_focei.npars, aEpsC.begin());
  thetaReset["aEpsC"] = aEpsC;
  thetaReset["c1"] = op_focei.c1;
  thetaReset["c2"] = op_focei.c2;
  parHistData(thetaReset, true);
}

static inline bool isFixedTheta(int m) {
  unsigned int j;
  for (unsigned int k = op_focei.npars; k--;){
    j=op_focei.fixedTrans[k];
    if (m == j) return false; // here the parameter is estimated
  }
  return true; // here the parameter is fixed
}

static inline bool thetaReset0() {
  NumericVector thetaIni(op_focei.ntheta);
  NumericVector thetaUp(op_focei.ntheta);
  NumericVector thetaDown(op_focei.ntheta);
  LogicalVector adjustEta(op_focei.muRefN);
  bool doAdjust = false;
  for (int ii = op_focei.ntheta; ii--;){
    thetaIni[ii] = unscalePar(op_focei.fullTheta, ii);
    if (R_FINITE(op_focei.lower[ii])){
      thetaDown[ii] = unscalePar(op_focei.lower, ii);
    } else {
      thetaDown[ii] = R_NegInf;
    }
    if (R_FINITE(op_focei.upper[ii])) {
      thetaUp[ii]= unscalePar(op_focei.upper, ii);
    } else {
      thetaUp[ii] = R_PosInf;
    }
  }
  double ref=0;
  int ij = 0;
  for (int ii = op_focei.muRefN; ii--;){
    if (op_focei.muRef[ii] != -1 && op_focei.muRef[ii] < (int)op_focei.ntheta) {
      ij = op_focei.muRef[ii];
      if (isFixedTheta(ij)) {
	adjustEta[ii] = false;
      }  else {
	ref = thetaIni[ij] + op_focei.etaM(ii,0);
	if (thetaDown[ij] < ref && thetaUp[ij] > ref) {
	  thetaIni[ij] = ref;
	  adjustEta[ii] = true;
	  doAdjust = true;
	} else {
	  adjustEta[ii] = false;
	}
      }
    } else {
      adjustEta[ii] = false;
    }
  }
  if (!doAdjust) return false;

  arma::mat etaMat(rx->nsub, op_focei.neta);
  for (int ii = rx->nsub; ii--;){
    focei_ind *fInd = &(inds_focei[ii]);
    for (int jj = op_focei.neta; jj--; ){
      if (op_focei.muRef[jj] != -1  && op_focei.muRef[jj] < (int)op_focei.ntheta &&
	  adjustEta[jj]){
	etaMat(ii, jj) = fInd->eta[jj]-op_focei.etaM(jj,0);
      } else {
	etaMat(ii, jj) = fInd->eta[jj];
      }
    }
  }
  // Update omega estimates
  NumericVector omegaTheta(op_focei.omegan);

  std::copy(&op_focei.fullTheta[0] + op_focei.ntheta,
	    &op_focei.fullTheta[0] + op_focei.ntheta + op_focei.omegan,
	    omegaTheta.begin());
  thetaReset00(thetaIni, omegaTheta, etaMat);
  return true;
}

void thetaReset(double size){
  mat etaRes =  op_focei.eta1SD % op_focei.etaM; //op_focei.cholOmegaInv * etaMat;
  for (unsigned int j = etaRes.n_rows; j--;) {
    if (std::fabs(etaRes(j, 0)) >= size) { // Says reset;
      if (thetaReset0()) {
	if (op_focei.didEtaReset==1) {
	  warning(_("mu-referenced Thetas were reset during optimization; (Can control by foceiControl(resetThetaP=.,resetThetaCheckPer=.,resetThetaFinalP=.))"));
	}
	stop("theta reset");
      }
    }
  }
}

void thetaResetZero(){
  thetaReset0();
  warning(_("thetas were reset during optimization because of a zero gradient"));
  stop("theta reset0");
}

void thetaResetObj(Environment e) {
  // Check to see if the last objective function is the lowest, otherwise reset to the lowest thetaResetObj
  if (op_focei.maxOuterIterations > 0) {
    parHistData(e, true);
    if (e.exists("parHistData")) {
      if (TYPEOF(e["parHistData"]) == VECSXP) {
	List parHistData = e["parHistData"];
	IntegerVector iter = parHistData["iter"];
	IntegerVector type = parHistData["type"];
	NumericVector obj = parHistData["objf"];
	int maxiter=-1, minObjId=-1;
	double minObj = R_PosInf;
	for (int i = obj.size(); i--;) {
	  if (type[i] == 5) {
	    if (!ISNA(obj[i])) {
	      if (obj[i] < minObj) {
		minObj = obj[i];
		minObjId = i;
	      }
	    }
	    maxiter = max2(maxiter, iter[i]);
	  }
	}
	if (iter[minObjId] != maxiter) {
	  // Min objective function is not at the last value
	  // REprintf("not at minimum objective function seen\n");
	  // NumericVector thetaIni(op_focei.ntheta);
	  // NumericVector omegaTheta(op_focei.omegan);
	  // for (int j = op_focei.ntheta; j--;){
	  //   NumericVector cur = parHistData[3+j];
	  //   thetaIni[j] = cur[minObjId];
	  // }
	  // for (int j = op_focei.omegan; j--;) {
	  //   NumericVector cur = parHistData[3+op_focei.ntheta+j];
	  //   omegaTheta[j] = cur[minObjId];
	  // }
	  // print(wrap(thetaIni));
	  // print(wrap(omegaTheta));
	  // arma::mat etaMat(rx->nsub, op_focei.neta, arma::fill::zeros);
	  // thetaReset00(thetaIni, omegaTheta, etaMat);
	  warning(_("last objective function was not at minimum, possible problems in optimization"));
	  // stop("theta resetZ");
	}
      }
    }
  }
}

static inline int didInnerResetPointFail(focei_ind *indF, int& id, double point) {
   indF->mode = 1;
   indF->uzm = 1;
   op_focei.didHessianReset=1;
   if (point == 0.0) {
     std::fill(&indF->eta[0], &indF->eta[0] + op_focei.neta, point);
   } else {
     // chol(omega^-1) %*% eta = point for each eta
     //
     //mat etaRes = op_focei.cholOmegaInv * etaMat
     for (int j = op_focei.neta; j--; ) {
       indF->eta[j] =  point/op_focei.cholOmegaInv(j, j);
     }
   }
   return !innerOpt1(id, 0);
}

static inline int didInnerResetFailTryOne(focei_ind *indF, int& id) {
  if (didInnerResetPointFail(indF, id, 0.0)) {
    if (op_focei.etaNudge != 0.0) {
      if (didInnerResetPointFail(indF, id, op_focei.etaNudge) &&
	  didInnerResetPointFail(indF, id, -op_focei.etaNudge)) {
	if (op_focei.etaNudge2 != 0.0) {
	  if (didInnerResetPointFail(indF, id, op_focei.etaNudge2) &&
	      didInnerResetPointFail(indF, id, -op_focei.etaNudge2)) {
	    return 1; // failed
	  } else {
	    return 0; // reset on one of the eta nudge
	  }
	} else {
	  return 1; // failed reset
	}
      } else {
	return 0; // reset
      }
    } else {
      return 1; // failed reset
    }
  }
  return 0;
}

static inline int didInnerResetFail(focei_ind *indF, int& id) {
  if (didInnerResetFailTryOne(indF, id)) {
    if (op_focei.canDoFD) {
      indF->doFD = 1;
      if (didInnerResetFailTryOne(indF, id)) {
	indF->doFD = 0;
	return 1; // failed reset
      } else {
	indF->doFD = 0;
	return 0;
      }
    } else {
      return 1;
    }
  }
  return 0;
}

static inline void resetToZeroWithoutOpt(focei_ind *indF, int& id) {
  // Just use ETA=0
  std::fill(&indF->eta[0], &indF->eta[0] + op_focei.neta, 0.0);
  if (!innerEval(id)) {
    indF->lik[0] = NA_REAL;
    warning(_("bad solve during optimization"));
  }
}

static inline void innerOptId(int id) {
  focei_ind *indF = &(inds_focei[id]);
  if (!innerOpt1(id, 0)) {
    // First try resetting ETA
    if (didInnerResetFail(indF, id)) {
      if(!op_focei.noabort){
	stop("Could not find the best eta even hessian reset and eta reset for ID %d.", id+1);
      } else if (indF->doChol == 1){
	indF->doChol = 0; // Use generalized cholesky decomposition
	if (didInnerResetFail(indF, id)) {
	  resetToZeroWithoutOpt(indF, id);
	}
	indF->doChol = 1; // Use cholesky again.
      } else {
	resetToZeroWithoutOpt(indF, id);
      }
    }
  }
}



void innerOpt(){
// #ifdef _OPENMP
//   int cores = rx->op->cores;
// #endif
  rx = getRx();
  if (op_focei.neta > 0) {
    op_focei.omegaInv=getOmegaInv();
    op_focei.logDetOmegaInv5 = getOmegaDet();
  }
  if (op_focei.maxInnerIterations <= 0){
    std::fill_n(&op_focei.goldEta[0], op_focei.gEtaGTransN, -42.0); // All etas = -42;  Unlikely if normal
    // #ifdef _OPENMP
// #pragma omp parallel for num_threads(cores)
// #endif
    for (int id = 0; id < rx->nsub; id++){
      focei_ind *indF = &(inds_focei[id]);
      indF->doChol = 1;
      if (!innerEval(id)) {
	indF->doChol = 0; // Use generalized cholesky decomposition
        innerEval(id);
	// Not thread safe
	warning(_("non-positive definite individual Hessian at solution(ID=%d); FOCEi objective functions may not be comparable"),id);
        indF->doChol = 1; // Cholesky again.
      }
    }
  } else {
// #ifdef _OPENMP
// #pragma omp parallel for num_threads(cores)
// #endif
    for (int id = 0; id < rx->nsub; id++){
      innerOptId(id);
    }
    // Reset ETA variances for next step
    if (op_focei.neta > 0){
      if (op_focei.zeroGrad) {
	thetaResetZero();
      }
      op_focei.eta1SD = 1/sqrt(op_focei.etaS);
      if (!op_focei.calcGrad && op_focei.maxOuterIterations > 0 &&
	  (!op_focei.initObj || op_focei.checkTheta==1) &&
	  R_FINITE(op_focei.resetThetaSize)){
	// Not thread safe...
	thetaReset(op_focei.resetThetaSize);
      }
      std::fill(op_focei.etaM.begin(),op_focei.etaM.end(), 0.0);
      std::fill(op_focei.etaS.begin(),op_focei.etaS.end(), 0.0);
      op_focei.n = 0.0;
    }
  }
  Rcpp::checkUserInterrupt();
}

static inline double foceiLik0(double *theta){
  updateTheta(theta);
  innerOpt();
  double lik = 0.0;
  double cur;

  for (int id=rx->nsub; id--;){
    focei_ind *fInd = &(inds_focei[id]);
    cur = fInd->lik[0];
    // if (std::isnan(cur)) {
    //   REprintf(_("likelihood of id: %s is NaN\n"), rxGetId(id));
    // } else if (std::isinf(cur)) {
    //   REprintf(_("likelihood of id: %s is infinite\n"), rxGetId(id));
    // }
    if (ISNA(cur) || std::isinf(cur) || std::isnan(cur)) {
      cur = -op_focei.badSolveObjfAdj;
    }
    lik += cur;
  }
  // Now reset the saved ETAs
  if (op_focei.neta !=0) std::fill_n(&op_focei.goldEta[0], op_focei.gEtaGTransN, -42.0); // All etas = -42;  Unlikely if normal
  return lik;
}


static inline double foceiOfv0(double *theta){
  if (op_focei.objfRecalN != 0 && !op_focei.calcGrad) {
    op_focei.stickyRecalcN1++;
    if (op_focei.stickyRecalcN1 <= op_focei.stickyRecalcN){
      RxODE::atolRtolFactor_(pow(op_focei.odeRecalcFactor, -op_focei.objfRecalN));
    } else {
      op_focei.stickyTol=1;
    }
  }
  double ret = -2*foceiLik0(theta);
  while (!op_focei.calcGrad && op_focei.stickyRecalcN1 <= op_focei.stickyRecalcN &&
	 (std::isnan(ret) || std::isinf(ret)) &&
	 op_focei.objfRecalN < op_focei.maxOdeRecalc){
      op_focei.reducedTol=1;
      RxODE::atolRtolFactor_(op_focei.odeRecalcFactor);
      ret = -2*foceiLik0(theta);
      op_focei.objfRecalN++;
  }
  if (!op_focei.initObj){
    op_focei.initObj=1;
    op_focei.initObjective=std::fabs(ret);
    if (std::isnan(ret)){
      stop(_("NaN while evaluating initial objective function"));
    } else if (std::isinf(ret)) {
      stop(_("infinite while evaluating initial objective function"));
    }
    if (op_focei.scaleObjective == 1) op_focei.scaleObjective=2;
  } else {
    if (std::isnan(ret) || std::isinf(ret)){
      ret=5e100;
    }
  }
  if (op_focei.scaleObjective == 2){
    ret = ret / op_focei.initObjective * op_focei.scaleObjectiveTo;
  }
  if (!op_focei.calcGrad){
    if (op_focei.derivMethodSwitch){
      double diff = std::fabs(op_focei.lastOfv-ret);
      if (op_focei.derivMethod==0 && diff <= op_focei.derivSwitchTol){
	op_focei.derivMethod=1;
      } else if (op_focei.derivMethod==1 && diff > op_focei.derivSwitchTol){
	op_focei.derivMethod=0;
      }
    }
    if (fabs((op_focei.lastOfv-ret)/max2(op_focei.lastOfv, ret))*100 < op_focei.resetThetaCheckPer){
      op_focei.checkTheta=1;
    } else {
      op_focei.checkTheta=0;
    }
    op_focei.lastOfv = ret;
  }
  return ret;
}

//[[Rcpp::export]]
double foceiLik(NumericVector theta){
  return foceiLik0(&theta[0]);
}

//[[Rcpp::export]]
double foceiOfv(NumericVector theta){
  return foceiOfv0(&theta[0]);
}

SEXP foceiEtas(Environment e) {
  if (op_focei.neta==0) return R_NilValue;
  List ret(op_focei.neta+2);
  CharacterVector nm(op_focei.neta+2);
  rx = getRx();
  IntegerVector ids(rx->nsub);
  NumericVector ofv(rx->nsub);
  int j,eta;
  for (j = op_focei.neta; j--;){
    ret[j+1]=NumericVector(rx->nsub);
    nm[j+1] = "ETA[" + std::to_string(j+1) + "]";
  }
  NumericVector tmp;
  for (j=rx->nsub; j--;){
    ids[j] = j+1;
    focei_ind *fInd = &(inds_focei[j]);
    ofv[j] = -2*fInd->lik[0];
    for (eta = op_focei.neta; eta--;){
      tmp = ret[eta+1];
      // Save eta is what the ETAs are saved
      tmp[j] = fInd->saveEta[eta];
    }
  }
  if (e.exists("IDlabel")) {
    RObject idl = e["IDlabel"];
    if (idl.sexp_type() == STRSXP) {
      ids.attr("class") = "factor";
      ids.attr("levels") = idl;
    }
  }
  ret[0] = ids;
  nm[0] = "ID";
  ret[op_focei.neta+1]=ofv;
  nm[op_focei.neta+1] = "OBJI";
  ret.attr("names") = nm;
  ret.attr("class") = "data.frame";
  ret.attr("row.names") = IntegerVector::create(NA_INTEGER,-rx->nsub);
  return(wrap(ret));
}


// R style optimfn
extern "C" double outerLikOpim(int n, double *par, void *ex){
  return(foceiOfv0(par));
}

// Gill 1983 Chat
static inline double Chat(double phi, double h, double epsA){
  if (phi == 0) return 2*epsA/(h*h);
   return 2*epsA/(h*fabs(phi));
}

static inline double ChatP(double phi, double h, double epsA){
  if (phi == 0) return 4*epsA/(h*h*h);
  return 4*epsA/(h*h*fabs(phi));
}

static inline double Phi(double fp, double f, double fn, double h){
  return (fp-2*f+fn)/(h*h);
}
static inline double phiC(double fp, double fn, double h){
  return (fp-fn)/(2*h);
}
static inline double phiF(double f, double fp, double h){
  return (fp-f)/h;
}
static inline double phiB(double f, double fn, double h){
  return (f-fn)/h;
}

// Call R function for gill83
int foceiGill=1;
double gillF=NA_REAL;
int gillThetaN=0;
Environment gillRfnE_;
Environment baseEnv = Environment::base_env();
Function doCall = baseEnv["do.call"];
Function gillRfn_ = baseEnv["invisible"];
int gillPar = 0;
double gillLong = false;
double gillRfn(double *theta){
  List par(1);
  NumericVector par0(gillThetaN);
  std::copy(&theta[0], &theta[0]+gillThetaN,par0.begin());
  par[0] = par0;
  NumericVector ret = as<NumericVector>(doCall(_["what"] = gillRfn_, _["args"]=par, _["envir"]=gillRfnE_));
  if (ret.size() == 1){
    return(ret[0]);
  } else {
    return(ret[gillPar]);
  }
}

static inline void gill83tickStep(int &k, int &K) {
  if (foceiGill == 1 && op_focei.slow) {
    if (k < K) {
      op_focei.cur += (K-k);
      op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
    }
  }
}

void gill83fnF(double *fp, double *theta, int) {
  if (foceiGill == 1) {
    updateTheta(theta);
    *fp = foceiOfv0(theta);
    if (op_focei.slow) op_focei.curTick = par_progress(op_focei.cur++, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
  } else {
    *fp = gillRfn(theta);
  }
}



// *hf is the forward difference final estimate
// *hphif is central difference final estimate (when switching from forward to central differences)
// *df is the derivative estimate
// *df2 is the 2nd derivative estimate, useful for pre-conditioning.
// *ef is the err of the final estimate.
// *theta is the theta vector
// cpar is the parameter we are considering
// epsR is the relative error for the problem
// K is the maximum number of iterations before giving up on searching for the best interval.
// Returns 1 -- Success
//         2 -- Large error; Derivative estimate error 50% or more of the derivative
//         3 -- Function constant or nearly constant for this parameter
//         4 -- Function odd or nearly linear, df = K, df2 ~ 0
//         5 -- df2 increases rapidly as h decreases
int gill83(double *hf, double *hphif, double *df, double *df2, double *ef,
	   double *theta, int cpar, double epsR, int K, double gillStep,
	   double fTol, int cid, gill83fn_type gill83fn) {
  if (foceiGill == 1) op_focei.calcGrad=1;

  double f , x, hbar, h0, fp, fn=NA_REAL, phif, phib, phic, phicc = 0, phi, Chf, Chb,
    Ch, hs, hphi, hk, tmp, ehat, lasth, lastht=NA_REAL, lastfpt=NA_REAL, phict=NA_REAL;
  if (foceiGill == 1) {
    f = op_focei.lastOfv;
  } else if (foceiGill == 2) {
    focei_ind *fInd = &(inds_focei[cid]);
    f = fInd->curF;
  } else {
    f = gillF;
  }
  int k = 0;
  // Relative error should be given by the tolerances, I believe.
  double epsA=std::fabs(f)*epsR;
  x = theta[cpar];
  // FD1: // Initialization
  hbar = 2*(1+std::fabs(x))*_safe_sqrt(epsA/(1+std::fabs(f)));
  h0 = gillStep*hbar;
  lasth=h0;
  theta[cpar] = x + h0;
  gill83fn(&fp, theta, cid);

  theta[cpar] = x-h0;
  gill83fn(&fn, theta, cid);

  phif = phiF(f, fp, h0);
  phib = phiB(f, fn, h0);
  phic = phiC(fp, fn, h0);
  phi = Phi(fp, f, fn, h0);
  Chf = Chat(phif, h0, epsA);
  Chb = Chat(phib, h0, epsA);
  Ch = ChatP(phi, h0, epsA);
  hs = -1;
  hphi=hbar; // Not defined in Gill, but used for central difference switch if there are problems
  // FD2:  // Decide if to accept the interval
  hk = h0;
  if (max2(Chf, Chb) <= 0.1){
      hs=h0;
  }
  if (0.001 <= Ch && Ch <= 0.1){
    phicc=phic;
    hphi=h0;
    if (fTol != 0 && fabs(phif) < fTol){
      lastfpt = fp;
      lastht  = lasth;
      phict=phic;
    }
    goto FD5;
  }
  if (fTol != 0 && fabs(phif) < fTol){
    lastfpt = fp;
    lastht  = lasth;
    phict=phic;
  }
  if (Ch < 0.001){
    goto FD4;
  }
 FD3: // Increase h
  k++;
  hk=hk*gillStep;
  lasth=hk;
  // Compute the associated finite difference estimates and their
  // relative condition errors.
  theta[cpar] = x + hk;
  gill83fn(&fp, theta, cid);
  theta[cpar] = x-hk;
  gill83fn(&fn, theta, cid);
  phif = phiF(f, fp, hk);
  phib = phiB(f, fn, hk);
  phic = phiC(fp, fn, hk);
  phi = Phi(fp, f, fn, hk);
  Chf = Chat(phif, hk, epsA);
  Chb = Chat(phib, hk, epsA);
  Ch = ChatP(phi, hk, epsA);
  if (hs < 0 && max2(Chf, Chb) <= 0.1){
    hs = hk;
  }
  if (Ch <= 0.1){
    phicc=phic;
    hphi = hk;
    if (fTol != 0 && fabs(phif) < fTol){
      lastfpt = fp;
      lastht  = lasth;
      phict=phic;
    }
    goto FD5;
  }
  if (fTol != 0 && fabs(phif) < fTol){
    lastfpt = fp;
    lastht  = lasth;
    phict=phic;
  }
  if (k == K) goto FD6;
  goto FD3;
 FD4: // Decrease h
  k++;
  hk=hk/gillStep;
  lasth=hk;
  // Compute the associated finite difference estimates and their
  // relative condition errors.
  theta[cpar] = x + hk;
  gill83fn(&fp, theta, cid);
  theta[cpar] = x-hk;
  gill83fn(&fn, theta, cid);
  phif = phiF(f, fp, hk);
  phib = phiB(f, fn, hk);
  tmp=phic;
  phic = phiC(fp, fn, hk);
  phi = Phi(fp, f, fn, hk);
  Chf = Chat(phif, hk, epsA);
  Chb = Chat(phib, hk, epsA);
  Ch = ChatP(phi, hk, epsA);
  if (Ch > .1){
    phicc=tmp;
    hphi=hk*gillStep; // hphi = h_k-1
    if (fTol != 0 && fabs(phif) < fTol){
      lastfpt = fp;
      lastht  = lasth;
      phict=phic;
    }
    goto FD5;
  }
  if (max2(Chf, Chb) <= 0.1){
    hs = hk;
  }
  if (0.001 <= Ch && Ch <= 1){
    hphi = hk;
    if (fTol != 0 && fabs(phif) < fTol){
      lastfpt = fp;
      lastht  = lasth;
      phict=phic;
    }
    goto FD5;
  }
  if (fTol != 0 && fabs(phif) < fTol){
    lastfpt = fp;
    lastht  = lasth;
    phict=phic;
  }
  if (k == K) goto FD6;
  goto FD4;
 FD5: // Compute the estimate of the optimal interval
  *df2 = phi;
  *hf = 2*_safe_sqrt(epsA/fabs(phi));
  theta[cpar] = x + *hf;
  gill83fn(&fp, theta, cid);
  // Restore theta
  theta[cpar] = x;
  if (foceiGill == 1) updateTheta(theta);
  else if (foceiGill == 2) updateEta(theta, cid);
  *df = phiF(f, fp, *hf);
  *ef = (*hf)*fabs(phi)/2+2*epsA/(*hf);
  *hphif=hphi;
  ehat = fabs(*df-phicc);
  if (max2(*ef, ehat) <= 0.5*(*df)){
    gill83tickStep(k, K);
    return 1;
  } else {
    // warning("The finite difference derivative err more than 50%% of the slope; Consider a different starting point.");
    if (!ISNA(lastht)){
      // Could be used;  Stick with the last below Ftol
      // *hf = lasth;
      // fp = lastfp;
      // *df = phiF(f, fp, *hf);
      // *df2=0;
      // // *df = 0.0; // Doesn't move.
      // *hphif=2*(*hf);
    // } else {
      *hf = lastht;
      fp = lastfpt;
      *df = phiF(f, fp, *hf);
      *df2=phic;
      // *df = 0.0; // Doesn't move.
      *hphif=phict;
    }
    gill83tickStep(k, K);
    return 2;
  }
 FD6: // Check unsatisfactory cases
  if (hs < 0){
    // F nearly constant.
    // Use sqrt(h0) as a last ditch effort.
    *hf = pow(DBL_EPSILON, 0.25);//hbar;
    // *df=phic;
    theta[cpar] = x + *hf;
    gill83fn(&fp, theta, cid);
    *df = phiF(f, fp, *hf);
    *df2=0;
    // *df = 0.0; // Doesn't move.
    *hphif=_safe_sqrt(h0);
    // warning("The surface around the initial estimate is nearly constant in one parameter grad=0.  Consider a different starting point.");
    gill83tickStep(k, K);
    return 3;
  }
  if (Ch > 0.1){ // Odd or nearly linear.
    *hf = h0;
    *df = phic;
    *df2 = 0;
    *ef = 2*epsA/(*hf);
    *hphif=hphi;
    // warning("The surface odd or nearly linear for one parameter; Check your function.");
    gill83tickStep(k, K);
    return 4;
  }
  // f'' is increasing rapidly as h decreases
  *hf = h0;
  *df = phic;
  *df2 = phi;
  *hphif=hphi;
  *ef = (*hf)*fabs(phi)/2+2*epsA/(*hf);
  // warning("The surface around the initial estimate is highly irregular in at least one parameter.  Consider a different starting point.");
  gill83tickStep(k, K);
  return 5;
}


void numericGrad(double *theta, double *g){
  op_focei.mixDeriv=0;
  op_focei.reducedTol2=0;
  op_focei.curGill=0;
  if ((op_focei.repeatGill == 1 || op_focei.nF + op_focei.nF2 == 1) && op_focei.gillK > 0){
    clock_t t = clock() - op_focei.t0;
    op_focei.slow = (op_focei.printOuter == 1) &&
      ((double)t)/CLOCKS_PER_SEC >= op_focei.gradProgressOfvTime;
    op_focei.repeatGill=0;
    op_focei.reducedTol2=0;
    double hf, hphif, err;
    if (op_focei.slow){
      op_focei.cur = 0;
      op_focei.totTick = op_focei.npars * op_focei.gillK;
      op_focei.t0 = clock();
      if (op_focei.repeatGillN != 0){
	RSprintf(_("repeat %d Gill diff/forward difference step size:\n"),
		op_focei.repeatGillN);
      } else {
	RSprintf(_("calculate Gill Difference and optimize forward difference step size:\n"));
      }
    }
    for (int cpar = op_focei.npars; cpar--;){
      op_focei.gillRet[cpar] = gill83(&hf, &hphif, &op_focei.gillDf[cpar], &op_focei.gillDf2[cpar], &op_focei.gillErr[cpar],
      				      theta, cpar, op_focei.gillRtol, op_focei.gillK, op_focei.gillStep, op_focei.gillFtol,
				      -1, gill83fnG);
      err = 1/(std::fabs(theta[cpar])+1);
      if (op_focei.gillDf[cpar] == 0){
	op_focei.scaleC[cpar]=op_focei.scaleC0;
	op_focei.gillRet[cpar] = gill83(&hf, &hphif, &op_focei.gillDf[cpar], &op_focei.gillDf2[cpar], &op_focei.gillErr[cpar],
					theta, cpar, op_focei.gillRtol, op_focei.gillK, op_focei.gillStep, op_focei.gillFtol,
					-1, gill83fnG);
	if (op_focei.gillDf[cpar] == 0){
	  op_focei.scaleC[cpar]=1/op_focei.scaleC0;
	  op_focei.gillRet[cpar] = gill83(&hf, &hphif, &op_focei.gillDf[cpar], &op_focei.gillDf2[cpar], &op_focei.gillErr[cpar],
					  theta, cpar, op_focei.gillRtol, op_focei.gillK, op_focei.gillStep, op_focei.gillFtol,
					  -1, gill83fnG);
	}
      }
      // h=aEps*(|x|+1)/sqrt(1+fabs(f));
      // h*sqrt(1+fabs(f))/(|x|+1) = aEps
      // let err=2*sqrt(epsA/(1+f))
      // err*(aEps+|x|rEps) = h
      // Let aEps = rEps (could be a different ratio)
      // h/err = aEps(1+|x|)
      // aEps=h/err/(1+|x|)
      //
      op_focei.aEps[cpar]  = hf*err;
      op_focei.rEps[cpar]  = hf*err;
      if(op_focei.optGillF){
	op_focei.aEpsC[cpar] = hf*err;
	op_focei.rEpsC[cpar] = hf*err;
      } else {
	op_focei.aEpsC[cpar] = hphif*err;
	op_focei.rEpsC[cpar] = hphif*err;
      }
      g[cpar] = op_focei.gillDf[cpar];
    }
    if(op_focei.slow){
      op_focei.cur=op_focei.totTick;
      op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
      RSprintf("\n");
    }
    op_focei.didGill=1;
    if (op_focei.reducedTol2 && op_focei.repeatGillN < op_focei.repeatGillMax){
      op_focei.repeatGill=1;
      op_focei.repeatGillN++;
      op_focei.reducedTol2=0;
    }
    op_focei.curGill=1;
  } else {
    if(op_focei.slow){
      op_focei.t0 = clock();
      op_focei.cur=0;
      op_focei.curTick=0;
      op_focei.totTick = op_focei.npars * 2;
    }
    op_focei.calcGrad=1;
    rx = getRx();
    int npars = op_focei.npars;
    int cpar;
    double cur, delta, tmp, tmp0=NA_REAL;
    double f=0;
    // Do Forward difference if the OBJF for *theta has already been calculated.
    bool doForward=false;
    if (op_focei.derivMethod == 0){
      doForward=true;
      // If the first derivative wasn't calculated, then calculate it.
      for (cpar = npars; cpar--;){
	if (theta[cpar] != op_focei.theta[cpar]){
	  doForward=false;
	  break;
	}
      }
      if (doForward){
	// Fill in lik0
	f=op_focei.lastOfv;
      } else {
	op_focei.calcGrad=0; // Save OBF and theta
	f = foceiOfv0(theta);
	if(op_focei.slow) op_focei.curTick = par_progress(op_focei.cur++, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
	op_focei.calcGrad=1;
	doForward=true;
      }
    }
    for (cpar = npars; cpar--;) {
      if (doForward){
	delta = (std::fabs(theta[cpar])*op_focei.rEps[cpar] + op_focei.aEps[cpar]);
      } else {
	delta = (std::fabs(theta[cpar])*op_focei.rEpsC[cpar] + op_focei.aEpsC[cpar]);
      }
      cur = theta[cpar];
      theta[cpar] = cur + delta;
      if (doForward){
	tmp = foceiOfv0(theta);
	if(op_focei.slow) op_focei.curTick = par_progress(op_focei.cur++, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
	g[cpar] = (tmp-f)/delta;
      } else {
	tmp0 = foceiOfv0(theta);
	if(op_focei.slow) op_focei.curTick = par_progress(op_focei.cur++, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
	theta[cpar] = cur - delta;
	tmp = foceiOfv0(theta);
	if(op_focei.slow) op_focei.curTick = par_progress(op_focei.cur++, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
	g[cpar] = (tmp0-tmp)/(2*delta);
      }
      if (doForward && fabs(g[cpar]) > op_focei.gradCalcCentralLarge){
	doForward = false;
	theta[cpar] = cur - delta;
	g[cpar] = (tmp-foceiOfv0(theta))/(2*delta);
	if(op_focei.slow) op_focei.curTick = par_progress(op_focei.cur++, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
	op_focei.mixDeriv=1;
      }
      if (std::isnan(g[cpar]) ||  ISNA(g[cpar]) || !R_FINITE(g[cpar])){
	if (doForward){
	  // Switch to Backward difference method
	  op_focei.mixDeriv=1;
	  theta[cpar] = cur - delta;
	  g[cpar] = (f-foceiOfv0(theta))/(delta);
	  if(op_focei.slow) op_focei.curTick = par_progress(op_focei.cur++, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
	  if (R_FINITE(op_focei.gradTrim)){
	    if (g[cpar] > op_focei.gradTrim){
	      g[cpar]=op_focei.gradTrim;
	    } else if (g[cpar] < op_focei.gradTrim){
	      g[cpar]=-op_focei.gradTrim;
	    }
	  }
	} else {
	  // We are using the central difference AND there is an NA in one of the terms
	  // g[cpar] = (tmp0-tmp)/(2*delta);
	  op_focei.mixDeriv=1;
	  if (std::isnan(tmp0) || ISNA(tmp0) || !R_FINITE(tmp0)){
	    // Backward
	    g[cpar] = (f-tmp)/delta;
	  } else {
	    // Forward
	    g[cpar] = (tmp0-f)/delta;
	  }
	  if (g[cpar] == 0.0)  {
	    double hf, hphif, err;
	    // Add a a Gill difference with zero gradient
	    op_focei.gillRet[cpar] = gill83(&hf, &hphif, &op_focei.gillDf[cpar], &op_focei.gillDf2[cpar], &op_focei.gillErr[cpar],
					    theta, cpar, op_focei.gillRtol, op_focei.gillK, op_focei.gillStep, op_focei.gillFtol,
					    -1, gill83fnG);
	    err = 1/(std::fabs(theta[cpar])+1);
	    if (op_focei.gillDf[cpar] == 0){
	      op_focei.scaleC[cpar]=op_focei.scaleC0;
	      op_focei.gillRet[cpar] = gill83(&hf, &hphif, &op_focei.gillDf[cpar], &op_focei.gillDf2[cpar], &op_focei.gillErr[cpar],
					      theta, cpar, op_focei.gillRtol, op_focei.gillK, op_focei.gillStep, op_focei.gillFtol,
					      -1, gill83fnG);
	      if (op_focei.gillDf[cpar] == 0){
		op_focei.scaleC[cpar]=1/op_focei.scaleC0;
		op_focei.gillRet[cpar] = gill83(&hf, &hphif, &op_focei.gillDf[cpar], &op_focei.gillDf2[cpar], &op_focei.gillErr[cpar],
						theta, cpar, op_focei.gillRtol, op_focei.gillK, op_focei.gillStep, op_focei.gillFtol,
						-1, gill83fnG);
	      }
	    }
	    // h=aEps*(|x|+1)/sqrt(1+fabs(f));
	    // h*sqrt(1+fabs(f))/(|x|+1) = aEps
	    // let err=2*sqrt(epsA/(1+f))
	    // err*(aEps+|x|rEps) = h
	    // Let aEps = rEps (could be a different ratio)
	    // h/err = aEps(1+|x|)
	    // aEps=h/err/(1+|x|)
	    //
	    op_focei.aEps[cpar]  = hf*err;
	    op_focei.rEps[cpar]  = hf*err;
	    if(op_focei.optGillF){
	      op_focei.aEpsC[cpar] = hf*err;
	      op_focei.rEpsC[cpar] = hf*err;
	    } else {
	      op_focei.aEpsC[cpar] = hphif*err;
	      op_focei.rEpsC[cpar] = hphif*err;
	    }
	    g[cpar] = op_focei.gillDf[cpar];
	  }
	  if (R_FINITE(op_focei.gradTrim)){
	    if (g[cpar] > op_focei.gradTrim){
	      g[cpar]=op_focei.gradTrim;
	    } else if (g[cpar] < op_focei.gradTrim){
	      g[cpar]=-op_focei.gradTrim;
	    } else if (std::isnan(tmp0) || ISNA(g[cpar])) {
	      g[cpar]=op_focei.gradTrim;
	    }
	  }
	}
      }
      if (R_FINITE(op_focei.gradTrim) && g[cpar] > op_focei.gradTrim){
	if (doForward){
	  op_focei.mixDeriv=1;
	  theta[cpar] = cur - delta;
	  g[cpar] = (tmp-foceiOfv0(theta))/(2*delta);
	  if(op_focei.slow)  op_focei.curTick = par_progress(op_focei.cur++, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
	  if (g[cpar] > op_focei.gradTrim){
	    g[cpar]=op_focei.gradTrim;
	  } else if (g[cpar] < op_focei.gradTrim){
	    g[cpar]=-op_focei.gradTrim;
	  }
	} else {
	  g[cpar]=op_focei.gradTrim;
	}
      } else if (R_FINITE(op_focei.gradTrim) && g[cpar] < -op_focei.gradTrim){
	if (doForward){
	  op_focei.mixDeriv=1;
	  theta[cpar] = cur - delta;
	  g[cpar] = (tmp-foceiOfv0(theta))/(2*delta);
	  if(op_focei.slow) op_focei.curTick = par_progress(op_focei.cur++, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
	  if (g[cpar] > op_focei.gradTrim){
	    g[cpar]=op_focei.gradTrim;
	  } else if (g[cpar] < op_focei.gradTrim){
	    g[cpar]=-op_focei.gradTrim;
	  }
	} else {
	  g[cpar]=-op_focei.gradTrim;
	}
      } else if (doForward && fabs(g[cpar]) < op_focei.gradCalcCentralSmall){
	op_focei.mixDeriv = 1;
	theta[cpar]       = cur - delta;
	tmp = g[cpar];
	g[cpar]           = (tmp-foceiOfv0(theta))/(2*delta);
	if(op_focei.slow) op_focei.curTick = par_progress(op_focei.cur++, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
	if (fabs(tmp) > fabs(g[cpar])) g[cpar] = tmp;
      } else if (doForward) {
	if(op_focei.slow) op_focei.curTick = par_progress(op_focei.cur++, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
      }
      theta[cpar] = cur;
      // Now for zero gradients try Gill differences
      if (g[cpar] == 0.0)  {
	op_focei.mixDeriv=1;
	double hf, hphif, err;
	// Add a a Gill difference with zero gradient
	op_focei.gillRet[cpar] = gill83(&hf, &hphif, &op_focei.gillDf[cpar], &op_focei.gillDf2[cpar], &op_focei.gillErr[cpar],
					theta, cpar, op_focei.gillRtol, op_focei.gillK, op_focei.gillStep, op_focei.gillFtol,
					-1, gill83fnG);
	err = 1/(std::fabs(theta[cpar])+1);
	if (op_focei.gillDf[cpar] == 0){
	  op_focei.scaleC[cpar]=op_focei.scaleC0;
	  op_focei.gillRet[cpar] = gill83(&hf, &hphif, &op_focei.gillDf[cpar], &op_focei.gillDf2[cpar], &op_focei.gillErr[cpar],
					  theta, cpar, op_focei.gillRtol, op_focei.gillK, op_focei.gillStep, op_focei.gillFtol,
					  -1, gill83fnG);
	  if (op_focei.gillDf[cpar] == 0){
	    op_focei.scaleC[cpar]=1/op_focei.scaleC0;
	    op_focei.gillRet[cpar] = gill83(&hf, &hphif, &op_focei.gillDf[cpar], &op_focei.gillDf2[cpar], &op_focei.gillErr[cpar],
					    theta, cpar, op_focei.gillRtol, op_focei.gillK, op_focei.gillStep, op_focei.gillFtol,
					    -1, gill83fnG);
	  }
	}
	// h=aEps*(|x|+1)/sqrt(1+fabs(f));
	// h*sqrt(1+fabs(f))/(|x|+1) = aEps
	// let err=2*sqrt(epsA/(1+f))
	// err*(aEps+|x|rEps) = h
	// Let aEps = rEps (could be a different ratio)
	// h/err = aEps(1+|x|)
	// aEps=h/err/(1+|x|)
	//
	op_focei.aEps[cpar]  = hf*err;
	op_focei.rEps[cpar]  = hf*err;
	if(op_focei.optGillF){
	  op_focei.aEpsC[cpar] = hf*err;
	  op_focei.rEpsC[cpar] = hf*err;
	} else {
	  op_focei.aEpsC[cpar] = hphif*err;
	  op_focei.rEpsC[cpar] = hphif*err;
	}
	g[cpar] = op_focei.gillDf[cpar];
      }
    }
    if(op_focei.slow) {
      op_focei.cur=op_focei.totTick;
      op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
      RSprintf("\n");
    }
    op_focei.calcGrad=0;
  }
}

//[[Rcpp::export]]
NumericVector foceiNumericGrad(NumericVector theta){
  NumericVector ret(theta.size());
  numericGrad(&theta[0], &ret[0]);
  return ret;
}

////////////////////////////////////////////////////////////////////////////////
// Setup FOCEi functions
CharacterVector rxParams_(const RObject &obj);


static inline void foceiSetupTrans_(CharacterVector pars){
  unsigned int k, j,  ps = pars.size();
  k=ps;
  std::string thetaS;
  std::string etaS;
  std::string cur;
  if (op_focei.etaTrans != NULL) R_Free(op_focei.etaTrans);
  op_focei.etaTrans    = R_Calloc(op_focei.neta*3 + 3*(op_focei.ntheta + op_focei.omegan), int); //[neta]
  op_focei.nbdInner    = op_focei.etaTrans + op_focei.neta;
  op_focei.xPar        = op_focei.nbdInner + op_focei.neta; // [ntheta+nomega]
  op_focei.thetaTrans  = op_focei.xPar + op_focei.ntheta + op_focei.omegan; // [ntheta+nomega]
  op_focei.fixedTrans  = op_focei.thetaTrans + op_focei.ntheta + op_focei.omegan; // [ntheta + nomega]
  op_focei.etaFD       = op_focei.fixedTrans + op_focei.ntheta + op_focei.omegan; // [neta]

  if (op_focei.fullTheta != NULL) R_Free(op_focei.fullTheta);
  op_focei.fullTheta   = R_Calloc(4*(op_focei.ntheta+op_focei.omegan), double); // [ntheta+omegan]
  op_focei.theta       = op_focei.fullTheta+op_focei.ntheta+op_focei.omegan; // [ntheta + omegan]
  op_focei.initPar     = op_focei.theta+op_focei.ntheta+op_focei.omegan; // [ntheta + omegan]
  op_focei.scaleC      = op_focei.initPar+op_focei.ntheta+op_focei.omegan; // [ntheta + omegan]

  int neta = 0;
  unsigned int ntheta = 0;
  for (;k--;){
    for (j = ps; j--;){
      // Compare ETAS first since they are smaller strings.
      cur = as<std::string>(pars[k]);
      etaS = "ETA[" + std::to_string(j+1) + "]";
      if (cur == etaS){
        op_focei.etaTrans[j] = k;
        neta++;
        break;
      } else {
        thetaS = "THETA[" + std::to_string(j+1) + "]";
        if (cur == thetaS){
          op_focei.thetaTrans[j] = k;
          ntheta++;
          break;
        }
      }
    }
  }
  if (op_focei.ntheta != ntheta){
    stop("theta mismatch op_focei.ntheta %d, ntheta: %d\n", op_focei.ntheta, ntheta);
  }
  if (op_focei.neta != neta){
    stop("eta mismatch op_focei.neta %d, neta: %d\n", op_focei.neta, neta);
  }
  op_focei.nzm = (op_focei.neta + 1) * (op_focei.neta + 2) / 2 + (op_focei.neta + 1)*6+1;
}

static inline void foceiSetupTheta_(List mvi,
				    NumericVector theta,
				    Nullable<LogicalVector> thetaFixed,
				    double scaleTo,
				    bool alloc){
  // Get the fixed thetas
  // fixedTrans gives the theta->full theta translation
  // initPar is the initial parameters used for parameter scaling.
  op_focei.scaleTo=scaleTo;
  int thetan = theta.size();
  int omegan;
  NumericVector omegaTheta;
  if (op_focei.neta == 0) {
    omegan=0;
  } else {
    omegan = getOmegaN();
    omegaTheta = getOmegaTheta();
  }
  int fixedn = 0;
  int j;
  LogicalVector thetaFixed2;
  if (!thetaFixed.isNull()){
    thetaFixed2 = as<LogicalVector>(thetaFixed);
    for (j = thetaFixed2.size(); j--;){
      if (thetaFixed2[j]) fixedn++;
    }
  } else {
    thetaFixed2 =LogicalVector(0);
  }
  int npars = thetan+omegan-fixedn;
  if (alloc){
    rxUpdateFuns(as<SEXP>(mvi["trans"]), &rxInner);
    foceiSetupTrans_(as<CharacterVector>(mvi["params"]));
  } else if (!op_focei.alloc){
    stop("FOCEi problem not allocated\nThis can happen when sympy<->nlmixr interaction is not working correctly.");
  }
  std::copy(theta.begin(), theta.end(), &op_focei.fullTheta[0]);
  if (op_focei.neta >= 0) {
    std::copy(omegaTheta.begin(), omegaTheta.end(), &op_focei.fullTheta[0]+thetan);
  }
  if ((int)op_focei.ntheta != (int)thetan){
    rxOptionsFreeFocei();
    stop("op_focei.ntheta(%d) != thetan(%d)", op_focei.ntheta, thetan);
  }
  op_focei.ntheta = thetan;
  op_focei.omegan = omegan;
  int k = 0;
  for (j = 0; j < npars+fixedn; j++){
    if (j < thetaFixed2.size() && !thetaFixed2[j]){
      if (j < theta.size()){
        op_focei.initPar[k] = theta[j];
      } else if (j < theta.size() + omegan){
        op_focei.initPar[k] = omegaTheta[j-theta.size()];
      }
      op_focei.fixedTrans[k++] = j;
    } else if (j >= thetaFixed2.size()){
      if (j < theta.size()){
        op_focei.initPar[k] = theta[j];
        op_focei.fixedTrans[k++] = j;
      } else if (j < theta.size() + omegan){
	op_focei.initPar[k] = omegaTheta[j-theta.size()];
	op_focei.fixedTrans[k++] = j;
      }
    }
  }
  op_focei.npars  = npars;
}

static inline void foceiSetupNoEta_(){
  rx = getRx();

  if (inds_focei != NULL) R_Free(inds_focei);
  inds_focei = R_Calloc(rx->nsub, focei_ind);
  op_focei.gEtaGTransN=(op_focei.neta)*rx->nsub;

  if (op_focei.gthetaGrad != NULL && op_focei.mGthetaGrad) R_Free(op_focei.gthetaGrad);
  op_focei.gthetaGrad = R_Calloc(op_focei.gEtaGTransN, double);
  op_focei.mGthetaGrad = true;
  focei_ind *fInd;
  int jj = 0;
  for (int i = rx->nsub; i--;){
    fInd = &(inds_focei[i]);
    fInd->doChol=!(op_focei.cholSEOpt);
    fInd->doFD=0;
    // ETA ini
    fInd->eta = NULL;
    fInd->oldEta = NULL;
    fInd->saveEta = NULL;
    fInd->g = NULL;
    fInd->x = NULL;
    fInd->var = NULL;
    fInd->lp = NULL;
    fInd->H = NULL;
    fInd->H0 = NULL;
    fInd->Vid = NULL;
    fInd->a = NULL;
    fInd->c = NULL;
    fInd->B = NULL;
    fInd->zm = NULL;
    fInd->thetaGrad = &op_focei.gthetaGrad[jj];
    jj+= op_focei.npars;
    fInd->mode = 1;
    fInd->uzm = 1;
    fInd->doEtaNudge=0;
  }
  op_focei.alloc=true;
}

static inline void foceiSetupEta_(NumericMatrix etaMat0){
  rx = getRx();

  if (inds_focei != NULL) R_Free(inds_focei);
  inds_focei = R_Calloc(rx->nsub, focei_ind);
  etaMat0 = transpose(etaMat0);
  op_focei.gEtaGTransN=(op_focei.neta+1)*rx->nsub;
  int nz = ((op_focei.neta+1)*(op_focei.neta+2)/2+6*(op_focei.neta+1)+1)*rx->nsub;

  if (op_focei.etaUpper != NULL) R_Free(op_focei.etaUpper);
  op_focei.etaUpper = R_Calloc(op_focei.gEtaGTransN*9+ op_focei.npars*(rx->nsub + 1)+nz+
			     2*op_focei.neta * rx->nall + rx->nall+ rx->nall*rx->nall +
			     op_focei.neta*5 + 3*op_focei.neta*op_focei.neta*rx->nsub, double);
  op_focei.etaLower =  op_focei.etaUpper + op_focei.neta;
  op_focei.geta = op_focei.etaLower + op_focei.neta;
  op_focei.goldEta = op_focei.geta + op_focei.gEtaGTransN;
  op_focei.getahf = op_focei.goldEta + op_focei.gEtaGTransN;
  op_focei.getahr = op_focei.getahf + op_focei.gEtaGTransN;
  op_focei.gsaveEta = op_focei.getahr + op_focei.gEtaGTransN;
  op_focei.gG = op_focei.gsaveEta + op_focei.gEtaGTransN;
  op_focei.gVar = op_focei.gG + op_focei.gEtaGTransN;
  op_focei.gX = op_focei.gVar + op_focei.gEtaGTransN;
  op_focei.glp = op_focei.gX + op_focei.gEtaGTransN;
  op_focei.gthetaGrad = op_focei.glp + op_focei.gEtaGTransN;  // op_focei.npars*(rx->nsub + 1)
  op_focei.gZm = op_focei.gthetaGrad + op_focei.npars*(rx->nsub + 1); // nz
  op_focei.ga  = op_focei.gZm + nz;//[op_focei.neta * rx->nall]
  op_focei.gc  = op_focei.ga + op_focei.neta * rx->nall;//[op_focei.neta * rx->nall]
  op_focei.gB  = op_focei.gc + op_focei.neta * rx->nall;//[rx->nall]
  op_focei.gH  = op_focei.gB + rx->nall; //[op_focei.neta*op_focei.neta*rx->nsub]
  op_focei.gH0  = op_focei.gB + op_focei.neta*op_focei.neta*rx->nsub; //[op_focei.neta*op_focei.neta*rx->nsub]
  op_focei.gVid = op_focei.gH0 + op_focei.neta*op_focei.neta*rx->nsub;
  double *ptr = op_focei.gB + rx->nall;
  // Could use .zeros() but since I used Calloc, they are already zero.
  op_focei.etaM = mat(ptr, op_focei.neta, 1, false, true);
  ptr += op_focei.neta;
  op_focei.etaS = mat(ptr, op_focei.neta, 1, false, true);
  ptr += op_focei.neta;
  op_focei.eta1SD = mat(ptr, op_focei.neta, 1, false, true);

  // Prefill to 0.1 or 10%
  std::fill_n(&op_focei.gVar[0], op_focei.gEtaGTransN, 0.1);
  std::fill_n(&op_focei.goldEta[0], op_focei.gEtaGTransN, -42.0); // All etas = -42;  Unlikely if normal


  unsigned int i, j = 0, k = 0, ii=0, jj = 0, iA=0, iB=0, iH=0, iVid=0;
  focei_ind *fInd;
  for (i = rx->nsub; i--;){
    fInd = &(inds_focei[i]);
    rx_solving_options_ind *ind = &(rx->subjects[i]);
    fInd->doChol=!(op_focei.cholSEOpt);
    fInd->doFD = 0;
    // ETA ini
    fInd->eta = &op_focei.geta[j];
    fInd->etahf = &op_focei.getahf[j];
    fInd->etahr = &op_focei.getahr[j];
    fInd->oldEta = &op_focei.goldEta[j];
    fInd->saveEta = &op_focei.gsaveEta[j];
    fInd->g = &op_focei.gG[j];
    fInd->x = &op_focei.gX[j];
    fInd->var = &op_focei.gVar[j];
    fInd->lp = &op_focei.glp[j];
    fInd->H = &op_focei.gH[iH];
    fInd->H0 = &op_focei.gH0[iH];
    fInd->Vid = &op_focei.gVid[iVid];
    iH += op_focei.neta*op_focei.neta;
    iVid += (ind->n_all_times - ind->ndoses - ind->nevid2)*(ind->n_all_times - ind->ndoses - ind->nevid2);

    // Copy in etaMat0 to the inital eta stored (0 if unspecified)
    // std::copy(&etaMat0[i*op_focei.neta], &etaMat0[(i+1)*op_focei.neta], &fInd->saveEta[0]);
    std::copy(&etaMat0[i*op_focei.neta], &etaMat0[(i+1)*op_focei.neta], &fInd->eta[0]);

    fInd->eta[op_focei.neta] = i;
    fInd->saveEta[op_focei.neta] = i;
    fInd->oldEta[op_focei.neta] = i;

    j+=op_focei.neta+1;

    k+=op_focei.neta;

    fInd->a = &op_focei.ga[iA];
    fInd->c = &op_focei.gc[iA];
    iA += op_focei.neta * (ind->n_all_times - ind->ndoses - ind->nevid2);

    fInd->B = &op_focei.gB[iB];
    iB += (ind->n_all_times - ind->ndoses - ind->nevid2);

    fInd->zm = &op_focei.gZm[ii];
    ii+=(op_focei.neta+1) * (op_focei.neta + 2) / 2 + 6*(op_focei.neta + 1)+1;

    fInd->thetaGrad = &op_focei.gthetaGrad[jj];
    jj+= op_focei.npars;

    fInd->mode = 1;
    fInd->uzm = 1;
    fInd->doEtaNudge=1;
  }
  op_focei.thetaGrad = &op_focei.gthetaGrad[jj];
  op_focei.alloc=true;
}

extern "C" double foceiOfvOptim(int n, double *x, void *ex);
extern "C" void outerGradNumOptim(int n, double *par, double *gr, void *ex);

// [[Rcpp::export]]
NumericVector foceiSetup_(const RObject &obj,
			  const RObject &data,
			  NumericVector theta,
			  Nullable<LogicalVector> thetaFixed = R_NilValue,
			  Nullable<LogicalVector> skipCov    = R_NilValue,
			  RObject rxInv                      = R_NilValue,
			  Nullable<NumericVector> lower      = R_NilValue,
			  Nullable<NumericVector> upper      = R_NilValue,
			  Nullable<NumericMatrix> etaMat     = R_NilValue,
			  Nullable<List> control             = R_NilValue){
  if (control.isNull()){
    stop("ODE options must be specified.");
  }
  List odeO = as<List>(control);
  // This fills in op_focei.neta
  List mvi;
  if (!Rf_isNull(obj)){
    if (!RxODE::rxDynLoad(obj)){
      stop("Cannot load RxODE dlls for this model.");
    }
    mvi = RxODE::rxModelVars_(obj);
  }
  rxOptionsFreeFocei();
  op_focei.mvi = mvi;

  op_focei.badSolveObjfAdj=fabs(as<double>(odeO["badSolveObjfAdj"]));

  op_focei.zeroGrad = false;
  op_focei.resetThetaCheckPer = as<double>(odeO["resetThetaCheckPer"]);
  op_focei.printTop = as<int>(odeO["printTop"]);
  op_focei.nF2 = as<int>(odeO["nF"]);
  if (op_focei.nF2 == 0){
    vGrad.clear();
    vPar.clear();
    iterType.clear();
    gradType.clear();
    niter.clear();
    niterGrad.clear();
  }
  op_focei.maxOuterIterations = as<int>(odeO["maxOuterIterations"]);
  op_focei.maxInnerIterations = as<int>(odeO["maxInnerIterations"]);
  op_focei.maxOdeRecalc = as<int>(odeO["maxOdeRecalc"]);
  op_focei.objfRecalN=0;
  op_focei.odeRecalcFactor = as<double>(odeO["odeRecalcFactor"]);
  op_focei.reducedTol = 0;
  op_focei.repeatGill=0;
  op_focei.repeatGillN=0;
  op_focei.repeatGillMax=as<int>(odeO["repeatGillMax"]);
  op_focei.stickyRecalcN=as<int>(odeO["stickyRecalcN"]);
  op_focei.neta = as<int>(odeO["neta"]);
  if (op_focei.neta != 0) {
    if (!RxODE::rxIs(rxInv, "rxSymInvCholEnv")){
      stop(_("Omega isn't in the proper format"));
    } else {
      _rxInv = as<List>(rxInv);
    }
  }
  op_focei.omegan = as<int>(odeO["nomega"]);
  op_focei.ntheta = as<int>(odeO["ntheta"]);
  // op_focei.ntheta = op_focei.ntheta;
  op_focei.nfixed = as<int>(odeO["nfixed"]);
  if (op_focei.maxOuterIterations <= 0){
    // No scaling.
    foceiSetupTheta_(mvi, theta, thetaFixed, 0.0, !Rf_isNull(obj));
    op_focei.scaleObjective=0;
  } else {
    foceiSetupTheta_(mvi, theta, thetaFixed, as<double>(odeO["scaleTo"]), !Rf_isNull(obj));
    op_focei.scaleObjectiveTo=as<double>(odeO["scaleObjective"]);
    if (op_focei.scaleObjectiveTo <= 0){
      op_focei.scaleObjective=0;
    } else {
      op_focei.scaleObjective=1;
    }
  }
  // First see if etaMat is null.
  NumericMatrix etaMat0;
  unsigned int nsub=0;
  // Find the number of IDs to create an etaMat
  List df = as<List>(data);
  CharacterVector dfN = df.names();
  int idn = -1;
  std::string cur;
  for (unsigned int j = dfN.size(); j--;){
    cur = as<std::string>(dfN[j]);
    if (cur == "ID" || cur == "id" || cur == "Id" || cur == "iD"){
      idn=j;
      break;
    }
  }
  if  (idn == -1){
    stop("Can't find ID in dataset.");
  }
  IntegerVector ids = as<IntegerVector>(df[idn]);
  int last = ids[ids.size()-1]-1;
  for (unsigned int j = ids.size(); j--;){
    if (last != ids[j]){
      last = ids[j];
      nsub++;
    }
  }
  if (op_focei.neta > 0) {
    etaMat0 = NumericMatrix(nsub, op_focei.neta);
    if (!etaMat.isNull()){
      NumericMatrix etaMat1 = NumericMatrix(etaMat);
      if (etaMat1.nrow() != (int)nsub){
	print(etaMat1);
	stop("The etaMat must have the same number of ETAs (rows) as subjects.");
      }
      if (etaMat1.ncol() != op_focei.neta){
	print(etaMat1);
	stop("The etaMat must have the same number of ETAs (cols) as the model.");
      }
      std::copy(etaMat1.begin(),etaMat1.end(),etaMat0.begin());
    }
  }
  List params(theta.size()+op_focei.neta);
  CharacterVector paramsNames(theta.size()+op_focei.neta);
  unsigned int j;
  for (j = theta.size();j--;){
    params[j] = NumericVector(nsub,theta[j]);
    if (theta.hasAttribute("names")){
      paramsNames[j] = (as<CharacterVector>(theta.names()))[j];
    } else {
      paramsNames[j] = "THETA[" + std::to_string(j + 1) + "]";
    }
  }
  if (op_focei.neta > 0) {
    bool hasDimn = etaMat0.hasAttribute("dimnames");
    CharacterVector dims;
    if (hasDimn){
      List diml = etaMat0.attr("dimnames");
      if (!Rf_isNull(as<RObject>(diml[1]))){
	dims = as<CharacterVector>(diml[1]);
      } else {
	hasDimn=false;
      }
    }
    for (j=op_focei.neta; j--;){
      params[j+theta.size()]= etaMat0(_, j);
      if (hasDimn){
	paramsNames[j+theta.size()] = dims[j];
      } else {
	paramsNames[j+theta.size()] = "ETA[" + std::to_string(j + 1) + "]";
      }
    }
  }
  params.names() = paramsNames;
  params.attr("class") = "data.frame";
  params.attr("row.names") = IntegerVector::create(NA_INTEGER,-nsub);
  // Now pre-fill parameters.
  if (!Rf_isNull(obj)) {
    List rxControl = odeO["rxControl"];
    rxControl[Rxc_cores] = IntegerVector::create(1); // #Force one core; parallelization needs to be taken care of here in inner.cpp
    RxODE::rxSolve_(obj, rxControl,
		    R_NilValue,//const Nullable<CharacterVector> &specParams =
		    R_NilValue,//const Nullable<List> &extraArgs =
		    as<RObject>(params),//const RObject &params =
		    data,//const RObject &events =
		    R_NilValue, // inits
		    1);//const int setupOnly = 0
    rx = getRx();
    if (op_focei.neta == 0) foceiSetupNoEta_();
    else foceiSetupEta_(etaMat0);
  }
  op_focei.epsilon=as<double>(odeO["epsilon"]);
  op_focei.nsim=as<int>(odeO["n1qn1nsim"]);
  op_focei.imp=0;
  op_focei.resetThetaSize = R_PosInf;
  // op_focei.printInner=as<int>(odeO["printInner"]);
  // if (op_focei.printInner < 0) op_focei.printInner = -op_focei.printInner;
  op_focei.printOuter=as<int>(odeO["print"]);
  if (op_focei.printOuter < 0) op_focei.printOuter = -op_focei.printOuter;
  // if (op_focei.printInner > 0){
  //   rx->op->cores=1;
  // }
  int totN=op_focei.ntheta + op_focei.omegan;
  NumericVector cEps=odeO["derivEps"];
  if (cEps.size() != 2){
    stop("derivEps must be 2 elements for determining forward difference step size.");
  }
  NumericVector covDerivEps=odeO["centralDerivEps"];
  if (covDerivEps.size() != 2){
    stop("centralDerivEps must be 2 elements for determining central difference step size.");
  }
  op_focei.derivMethod = as<int>(odeO["derivMethod"]);
  if (op_focei.derivMethod == 3){
    op_focei.derivMethod=0;
    op_focei.derivSwitchTol=as<double>(odeO["derivSwitchTol"]);
    op_focei.derivMethodSwitch=1;
  }

  IntegerVector muRef;
  if (odeO.containsElementNamed("focei.mu.ref")){
    if (RxODE::rxIs(odeO["focei.mu.ref"], "numeric")) {
      muRef = as<IntegerVector>(odeO["focei.mu.ref"]);
    }
  }
  if (muRef.size() == 0){
    op_focei.resetThetaSize = R_PosInf;
    op_focei.muRefN=0;
  } else{
    op_focei.muRefN=muRef.size();
  }

  IntegerVector skipCov1;
  if (skipCov.isNull()){
    op_focei.skipCovN = 0;
  } else {
    skipCov1 = as<IntegerVector>(skipCov);
    op_focei.skipCovN = skipCov1.size();
  }

  if (op_focei.gillRet != NULL) R_Free(op_focei.gillRet);
  op_focei.gillRet = R_Calloc(2*totN+op_focei.npars+
			    op_focei.muRefN + op_focei.skipCovN, int);
  op_focei.gillRetC= op_focei.gillRet + totN;
  op_focei.nbd     = op_focei.gillRetC + totN;//[op_focei.npars]
  op_focei.muRef   = op_focei.nbd + op_focei.npars; //[op_focei.muRefN]
  if (op_focei.muRefN) std::copy(&op_focei.muRef[0], &op_focei.muRef[0]+op_focei.muRefN, muRef.begin());

  op_focei.skipCov   = op_focei.muRef + op_focei.muRefN; //[op_focei.skipCovN]
  if (op_focei.skipCovN) std::copy(skipCov1.begin(),skipCov1.end(),op_focei.skipCov); //

  if (op_focei.gillDf != NULL) R_Free(op_focei.gillDf);
  op_focei.gillDf = R_Calloc(7*totN + 2*op_focei.npars + rx->nsub, double);
  op_focei.gillDf2 = op_focei.gillDf+totN;
  op_focei.gillErr = op_focei.gillDf2+totN;
  op_focei.rEps=op_focei.gillErr + totN;
  op_focei.aEps = op_focei.rEps + totN;
  op_focei.rEpsC = op_focei.aEps + totN;
  op_focei.aEpsC = op_focei.rEpsC + totN;
  op_focei.lower = op_focei.aEpsC + totN;
  op_focei.upper = op_focei.lower +op_focei.npars;
  op_focei.likSav = op_focei.upper + op_focei.npars;//[rx->nsub]

  if (op_focei.nF2 != 0){
    // Restore Gill information
    IntegerVector gillRetC =as<IntegerVector>(odeO["gillRetC"]);
    IntegerVector gillRet =as<IntegerVector>(odeO["gillRet"]);
    NumericVector gillDf = as<NumericVector>(odeO["gillDf"]);
    NumericVector gillDf2 = as<NumericVector>(odeO["gillDf2"]);
    NumericVector gillErr = as<NumericVector>(odeO["gillErr"]);
    NumericVector rEps = as<NumericVector>(odeO["rEps"]);
    NumericVector aEps = as<NumericVector>(odeO["aEps"]);
    NumericVector rEpsC = as<NumericVector>(odeO["rEpsC"]);
    NumericVector aEpsC = as<NumericVector>(odeO["aEpsC"]);
    std::copy(gillRetC.begin(), gillRetC.end(), &op_focei.gillRetC[0]);
    std::copy(gillRet.begin(), gillRet.end(), &op_focei.gillRet[0]);
    std::copy(gillDf.begin(), gillDf.end(), &op_focei.gillDf[0]);
    std::copy(gillDf2.begin(), gillDf2.end(), &op_focei.gillDf2[0]);
    std::copy(gillErr.begin(), gillErr.end(), &op_focei.gillErr[0]);
    std::copy(rEps.begin(), rEps.end(), &op_focei.rEps[0]);
    std::copy(aEps.begin(), aEps.end(), &op_focei.aEps[0]);
    std::copy(rEpsC.begin(), rEpsC.end(), &op_focei.rEpsC[0]);
    std::copy(aEpsC.begin(), aEpsC.end(), &op_focei.aEpsC[0]);
  } else {
    if (op_focei.derivMethod){
      std::fill_n(&op_focei.rEps[0], totN, std::fabs(cEps[0])/2.0);
      std::fill_n(&op_focei.aEps[0], totN, std::fabs(cEps[1])/2.0);
      std::fill_n(&op_focei.rEpsC[0], totN, std::fabs(covDerivEps[0])/2.0);
      std::fill_n(&op_focei.aEpsC[0], totN, std::fabs(covDerivEps[1])/2.0);
    } else {
      std::fill_n(&op_focei.rEps[0], totN, std::fabs(cEps[0]));
      std::fill_n(&op_focei.aEps[0], totN, std::fabs(cEps[1]));
      std::fill_n(&op_focei.rEpsC[0], totN, std::fabs(covDerivEps[0]));
      std::fill_n(&op_focei.aEpsC[0], totN, std::fabs(covDerivEps[1]));
    }
  }
  NumericVector lowerIn(totN);
  NumericVector upperIn(totN);
  if (lower.isNull()){
    std::fill_n(lowerIn.begin(), totN, R_NegInf);
  } else {
    NumericVector lower1=as<NumericVector>(lower);
    if (lower1.size() == 1){
      std::fill_n(lowerIn.begin(), totN, lower1[0]);
    } else if (lower1.size() < totN){
      std::copy(lower1.begin(), lower1.end(), lowerIn.begin());
      std::fill_n(lowerIn.begin()+lower1.size(), totN - lower1.size(), R_NegInf);
    } else if (lower1.size() > totN){
      warning(_("lower bound is larger than the number of parameters being estimated"));
      std::copy(lower1.begin(), lower1.begin()+totN, lowerIn.begin());
    } else {
      lowerIn = lower1;
    }
  }
  if (upper.isNull()){
    std::fill_n(upperIn.begin(), totN, R_PosInf);
  } else {
    NumericVector upper1=as<NumericVector>(upper);
    if (upper1.size() == 1){
      std::fill_n(upperIn.begin(), totN, upper1[0]);
    } else if (upper1.size() < totN){
      std::copy(upper1.begin(), upper1.end(), upperIn.begin());
      std::fill_n(upperIn.begin()+upper1.size(), totN - upper1.size(), R_PosInf);
    } else if (upper1.size() > totN){
      warning(_("upper bound is larger than the number of parameters being estimated"));
      std::copy(upper1.begin(), upper1.begin()+totN, upperIn.begin());
    } else {
      upperIn = upper1;
    }
  }
  op_focei.lowerIn =lowerIn;
  op_focei.upperIn =upperIn;

  std::fill_n(&op_focei.nbd[0], op_focei.npars, 0);

  NumericVector ret(op_focei.npars, op_focei.scaleTo);
  op_focei.calcGrad=0;

  // Outer options
  op_focei.outerOpt	=as<int>(odeO["outerOpt"]);
  // lbfgsb options
  op_focei.factr	= as<double>(odeO["lbfgsFactr"]);
  op_focei.pgtol	= as<double>(odeO["lbfgsPgtol"]);
  op_focei.lmm		= as<int>(odeO["lbfgsLmm"]);
  op_focei.covDerivMethod = as<int>(odeO["covDerivMethod"]);
  op_focei.covMethod = as<int>(odeO["covMethod"]);
  op_focei.eigen = as<int>(odeO["eigen"]);
  op_focei.ci=0.95;
  op_focei.sigdig=as<double>(odeO["sigdig"]);
  op_focei.useColor=as<int>(odeO["useColor"]);
  op_focei.boundTol=as<double>(odeO["boundTol"]);
  op_focei.printNcol=as<int>(odeO["printNcol"]);
  op_focei.noabort=as<int>(odeO["noAbort"]);
  op_focei.interaction=as<int>(odeO["interaction"]);
  op_focei.cholSEtol=as<double>(odeO["cholSEtol"]);
  op_focei.hessEps=as<double>(odeO["hessEps"]);
  op_focei.cholAccept=as<double>(odeO["cholAccept"]);
  op_focei.resetEtaSize=as<double>(odeO["resetEtaSize"]);
  op_focei.resetThetaSize=as<double>(odeO["resetThetaSize"]);
  op_focei.resetThetaFinalSize = as<double>(odeO["resetThetaFinalSize"]);

  op_focei.cholSEOpt=as<double>(odeO["cholSEOpt"]);
  op_focei.cholSECov=as<double>(odeO["cholSECov"]);
  op_focei.fo=as<int>(odeO["fo"]);
  if (op_focei.neta == 0) op_focei.fo = 0;
  if (op_focei.fo) op_focei.maxInnerIterations=0;

  op_focei.covTryHarder=as<int>(odeO["covTryHarder"]);
  op_focei.resetHessianAndEta=as<int>(odeO["resetHessianAndEta"]);
  op_focei.gillK = as<int>(odeO["gillK"]);
  op_focei.gillKcov = as<int>(odeO["gillKcov"]);
  op_focei.gillStep    = fabs(as<double>(odeO["gillStep"]));
  op_focei.gillStepCov = fabs(as<double>(odeO["gillStepCov"]));
  op_focei.gillFtol    = fabs(as<double>(odeO["gillFtol"]));
  op_focei.gillFtolCov = fabs(as<double>(odeO["gillFtolCov"]));
  op_focei.didGill = 0;
  op_focei.gillRtol = as<double>(odeO["gillRtol"]);
  op_focei.scaleType = as<int>(odeO["scaleType"]);
  op_focei.normType = as<int>(odeO["normType"]);
  op_focei.scaleC0=as<double>(odeO["scaleC0"]);
  op_focei.scaleCmin=as<double>(odeO["scaleCmin"]);
  op_focei.scaleCmax=as<double>(odeO["scaleCmax"]);
  op_focei.abstol=as<double>(odeO["abstol"]);
  op_focei.reltol=as<double>(odeO["reltol"]);
  op_focei.smatNorm=as<int>(odeO["smatNorm"]);
  op_focei.rmatNorm=as<int>(odeO["rmatNorm"]);
  op_focei.covGillF=as<int>(odeO["covGillF"]);
  op_focei.optGillF=as<int>(odeO["optGillF"]);
  op_focei.covSmall = as<double>(odeO["covSmall"]);
  op_focei.gradTrim = as<double>(odeO["gradTrim"]);
  op_focei.gradCalcCentralLarge = as<double>(odeO["gradCalcCentralLarge"]);
  op_focei.gradCalcCentralSmall = as<double>(odeO["gradCalcCentralSmall"]);
  op_focei.etaNudge = as<double>(odeO["etaNudge"]);
  op_focei.etaNudge2 = as<double>(odeO["etaNudge2"]);
  op_focei.eventFD = as<double>(odeO["eventFD"]);
  op_focei.eventType = as<int>(odeO["eventType"]);
  op_focei.predNeq = as<int>(odeO["predNeq"]);
  op_focei.gradProgressOfvTime = as<double>(odeO["gradProgressOfvTime"]);
  op_focei.initObj=0;
  op_focei.lastOfv=std::numeric_limits<double>::max();
  for (unsigned int k = op_focei.npars; k--;){
    j=op_focei.fixedTrans[k];
    ret[k] = op_focei.fullTheta[j];
    if (R_FINITE(lowerIn[j])){
      op_focei.lower[k] = lowerIn[j];
      op_focei.lower[k] += 2*(op_focei.lower[k]*op_focei.rEps[k] + op_focei.aEps[k]);
      // lower bound only = 1
      op_focei.nbd[k]=1;
    } else {
      op_focei.lower[k] = R_NegInf;//std::numeric_limits<double>::lowest();
    }
    if (R_FINITE(upperIn[j])){
      op_focei.upper[k] = upperIn[j];
      op_focei.upper[k] -= 2*(op_focei.upper[k]*op_focei.rEps[k] - op_focei.aEps[k]);
      // Upper bound only = 3
      // Upper and lower bound = 2
      op_focei.nbd[k]= 3 - op_focei.nbd[j];
    } else {
      op_focei.upper[k] = R_PosInf;//std::numeric_limits<double>::max();
    }
  }
  double mn = op_focei.initPar[op_focei.npars-1], mx=op_focei.initPar[op_focei.npars-1],mean=0, oN=0, oM=0,s=0;
  double len=0;
  unsigned int k;
  if (op_focei.nF2 > 0 && odeO.containsElementNamed("c1") && odeO.containsElementNamed("c2")){
    op_focei.c1 = odeO["c1"];
    op_focei.c2 = odeO["c2"];
  } else {
    switch (op_focei.normType){
    case 1:
      // OptdesX
      // http://apmonitor.com/me575/uploads/Main/optimization_book.pdf
      for (k = op_focei.npars-1; k--;){
	mn = min2(op_focei.initPar[k],mn);
	mx = max2(op_focei.initPar[k],mx);
      }
      op_focei.c1 = (mx+mn)/2;
      op_focei.c2 = (mx-mn)/2;
      break;
    case 2: // Rescaling (min-max normalization)
      for (k = op_focei.npars-1; k--;){
	mn = min2(op_focei.initPar[k],mn);
	mx = max2(op_focei.initPar[k],mx);
      }
      op_focei.c1 = mn;
      op_focei.c2 = (mx-mn);
      break;
    case 3: // Mean normalization
      for (k = op_focei.npars-1; k--;){
	mn = min2(op_focei.initPar[k],mn);
	mx = max2(op_focei.initPar[k],mx);
	oN++;
	mean += (op_focei.initPar[k]-mean)/oN;
      }
      op_focei.c1 = mean;
      op_focei.c2 = (mx-mn);
      break;
    case 4: // Standardization
      for (k = op_focei.npars-1; k--;){
	mn = min2(op_focei.initPar[k],mn);
	mx = max2(op_focei.initPar[k],mx);
	oM= mean;
	oN++;
	mean += (op_focei.initPar[k]-mean)/oN;
	s += (op_focei.initPar[k]-mean)*(op_focei.initPar[k]-oM);
      }
      op_focei.c1 = mean;
      op_focei.c2 = _safe_sqrt(s/(oN-1));
      break;
    case 5: // Normalize to length.
      for (unsigned int k = op_focei.npars-1; k--;){
	len += op_focei.initPar[k]*op_focei.initPar[k];
      }
      op_focei.c1 = 0;
      op_focei.c2 = _safe_sqrt(len);
      break;
    case 6:
      // No Normalization
      op_focei.c1 = 0;
      op_focei.c2 = 1;
      break;
    default:
      stop("Unrecognized normalization (normType=%d)",op_focei.normType);
    }

  }
    // }
  return ret;
}

LogicalVector nlmixrEnvSetup(Environment e, double fmin){
  if (e.exists("theta") && RxODE::rxIs(e["theta"], "data.frame") &&
      e.exists("omega") && e.exists("etaObf")) {
    bool mixed = RxODE::rxIs(e["omega"], "matrix") && RxODE::rxIs(e["etaObf"], "data.frame");
    int nobs2=0;
    if (e.exists("nobs2")){
      nobs2=as<int>(e["nobs2"]);
    } else {
      nobs2 = rx->nobs2;
    }
    if (mixed){
      arma::mat omega = as<arma::mat>(e["omega"]);
      arma::mat D(omega.n_rows,omega.n_rows,fill::zeros);
      arma::mat cor(omega.n_rows,omega.n_rows);
      D.diag() = (sqrt(omega.diag()));
      arma::vec sd=D.diag();
      D = inv_sympd(D);
      cor = D * omega * D;
      cor.diag()= sd;
      e["omegaR"] = wrap(cor);
    }
    if (op_focei.scaleObjective){
      fmin = fmin * op_focei.initObjective / op_focei.scaleObjectiveTo;
    }
    bool doAdj = false;
    bool doObf = false;
    if (!e.exists("objective")){
      e["objective"] = fmin;
      if (as<bool>(e["adjLik"])){
	doAdj = true;
      }
    } else {
      fmin = as<double>(e["objective"]);
      if (e.exists("adjObf") && as<bool>(e["adjObf"])){
	doObf=true;
      }
    }
    e["OBJF"] = fmin;
    e["objf"] = fmin;
    NumericVector logLik(1);
    double adj= 0;
    if (doAdj){
      adj=nobs2*log(2*M_PI)/2;
    }
    e["adj"]=adj;
    logLik[0]=-fmin/2-adj;
    logLik.attr("df") = op_focei.npars;
    if (e.exists("nobs")){
      logLik.attr("nobs") = e["nobs"];
      e["BIC"] = fmin+2*adj + log(as<double>(e["nobs"]))*op_focei.npars;
    } else {
      logLik.attr("nobs") = nobs2;
      e["BIC"] = fmin + 2*adj + log((double)nobs2)*op_focei.npars;
      e["nobs"] = rx->nobs;
    }
    logLik.attr("class") = "logLik";
    e["logLik"] = logLik;

    e["AIC"] = fmin+2*adj+2*op_focei.npars;
    if (doObf){
      // -2 * object$logLik - object$dim$N * log(2 * pi)
      adj = -2*as<double>(logLik) - (nobs2)*log(2*M_PI);
      e["OBJF"] = adj;
      e["objf"] = adj;
      e["objective"] = adj;
    }
    return true;
  } else {
    stop("Not Setup right.........");
    return false;
  }
}

void foceiOuterFinal(double *x, Environment e){
  double fmin = foceiOfv0(x);

  NumericVector theta(op_focei.ntheta);
  std::copy(&op_focei.fullTheta[0],  &op_focei.fullTheta[0] + op_focei.ntheta,
            theta.begin());

  NumericVector fullTheta(op_focei.ntheta+op_focei.omegan);
  std::copy(&op_focei.fullTheta[0],  &op_focei.fullTheta[0] + op_focei.ntheta + op_focei.omegan,
            fullTheta.begin());
  LogicalVector thetaFixed(op_focei.ntheta);
  std::fill_n(thetaFixed.begin(),op_focei.ntheta, true);
  int j;
  for (int k = op_focei.npars; k--;){
    j=op_focei.fixedTrans[k];
    if (j < thetaFixed.size()) thetaFixed[j]=false;
  }
  // std::copy(&op_focei.thetaFixed[0],  &op_focei.thetaFixed[0] + op_focei.ntheta,
  //           thetaFixed.begin());
  NumericVector lowerIn(op_focei.ntheta);
  NumericVector upperIn(op_focei.ntheta);
  std::copy(&op_focei.lowerIn[0],  &op_focei.lowerIn[0] + op_focei.ntheta,
            lowerIn.begin());
  std::copy(&op_focei.upperIn[0],  &op_focei.upperIn[0] + op_focei.ntheta,
            upperIn.begin());
  e["theta"] = DataFrame::create(_["lower"]=lowerIn, _["theta"]=theta, _["upper"]=upperIn,
				 _["fixed"]=thetaFixed);
  e["fullTheta"] = fullTheta;
  if (op_focei.neta == 0){
    e["omega"] = R_NilValue;
    e["etaObf"] = R_NilValue;
  } else {
    e["omega"] = getOmega();
    e["etaObf"] = foceiEtas(e);
  }
  nlmixrEnvSetup(e, fmin);
}

static inline void foceiPrintLine(int ncol){
  RSprintf("|-----+---------------+");
  for (int i = 0; i < ncol; i++){
    if (i == ncol-1)
      RSprintf("-----------|");
    else
      RSprintf("-----------+");
  }
  RSprintf("\n");
}

////////////////////////////////////////////////////////////////////////////////
// Outer l-BFGS-b from R
extern "C" double foceiOfvOptim(int n, double *x, void *ex){
  double ret = foceiOfv0(x);
  niter.push_back(op_focei.nF2+(++op_focei.nF));
  // Scaled
  vPar.push_back(ret);
  iterType.push_back(5);
  int finalize = 0, i = 0;
  for (i = 0; i < n; i++){
    vPar.push_back(x[i]);
  }
  // Unscaled
  iterType.push_back(6);
  niter.push_back(niter.back());
  if (op_focei.scaleObjective){
    vPar.push_back(op_focei.initObjective * ret / op_focei.scaleObjectiveTo);
  } else {
    vPar.push_back(ret);
  }
  for (i = 0; i < n; i++){
    vPar.push_back(unscalePar(x, i));
  }
  // Back-transformed (7)
  iterType.push_back(7);
  niter.push_back(niter.back());
  if (op_focei.scaleObjective){
    vPar.push_back(op_focei.initObjective * ret / op_focei.scaleObjectiveTo);
  } else {
    vPar.push_back(ret);
  }
  for (i = 0; i < n; i++){
    if (op_focei.xPar[i] == 1){
      vPar.push_back(exp(unscalePar(x, i)));
    } else if (op_focei.xPar[i] < 0){
      int m = -op_focei.xPar[i]-1;
      vPar.push_back(expit(unscalePar(x, i), op_focei.logitThetaLow[m], op_focei.logitThetaHi[m]));
    } else {
      vPar.push_back(unscalePar(x, i));
    }
  }
  if (op_focei.printOuter != 0 && op_focei.nF % op_focei.printOuter == 0){
    if (op_focei.useColor && !isRstudio())
      RSprintf("|\033[1m%5d\033[0m|%#14.8g |", op_focei.nF+op_focei.nF2, ret);
    else
      RSprintf("|%5d|%#14.8g |", op_focei.nF+op_focei.nF2, ret);
    for (i = 0; i < n; i++){
      RSprintf("%#10.4g |", x[i]);
      if ((i + 1) != n && (i + 1) % op_focei.printNcol == 0){
        if (op_focei.useColor && op_focei.printNcol + i  > n){
          RSprintf("\n\033[4m|.....................|");
        } else {
          RSprintf("\n|.....................|");
        }
	finalize=1;
      }
    }
    if (finalize){
      while(true){
        if ((i++) % op_focei.printNcol == 0){
          if (op_focei.useColor) RSprintf("\033[0m");
          RSprintf("\n");
          break;
        } else {
          RSprintf("...........|");
        }
      }
    } else {
      RSprintf("\n");
    }
    if (op_focei.scaleObjective){
      RSprintf("|    U|%14.8g |", op_focei.initObjective * ret / op_focei.scaleObjectiveTo);
    } else {
      RSprintf("|    U|%14.8g |", ret);
    }
    for (i = 0; i < n; i++){
      // new  = (theta[k] - op_focei.scaleTo)*op_focei.scaleC[k] +  op_focei.initPar[k]
      // (new-ini)/c+scaleTo = theta[]
      RSprintf("%#10.4g |", unscalePar(x, i));
      if ((i + 1) != n && (i + 1) % op_focei.printNcol == 0){
	if (op_focei.useColor && op_focei.printNcol + i  > op_focei.npars){
	  RSprintf("\n\033[4m|.....................|");
	} else {
	  RSprintf("\n|.....................|");
	}
      }
    }
    if (finalize){
      while(true){
	if ((i++) % op_focei.printNcol == 0){
	  if (op_focei.useColor) RSprintf("\033[0m");
	  RSprintf("\n");
	  break;
	} else {
	  RSprintf("...........|");
	}
      }
    } else {
      RSprintf("\n");
    }
    if (op_focei.scaleObjective){
      if (op_focei.useColor && !isRstudio()){
	RSprintf("|    X|\033[1m%14.8g\033[0m |", op_focei.initObjective * ret / op_focei.scaleObjectiveTo);
      } else {
	RSprintf("|    X|%14.8g |", op_focei.initObjective * ret / op_focei.scaleObjectiveTo);
      }
    } else {
      if (op_focei.useColor && !isRstudio())
	RSprintf("|    X|\033[1m%14.8g\033[0m |", ret);
      else
	RSprintf("|    X|%14.8g |", ret);
    }
    for (i = 0; i < n; i++){
      if (op_focei.xPar[i] == 1){
	RSprintf("%#10.4g |", exp(unscalePar(x, i)));
      } else if (op_focei.xPar[i] < 0){
	int m = -op_focei.xPar[i]-1;
	RSprintf("%#10.4g |", expit(unscalePar(x, i), op_focei.logitThetaLow[m], op_focei.logitThetaHi[m]));
      } else {
	RSprintf("%#10.4g |", unscalePar(x, i));
      }
      if ((i + 1) != n && (i + 1) % op_focei.printNcol == 0){
	if (op_focei.useColor && op_focei.printNcol + i >= op_focei.npars){
	  RSprintf("\n\033[4m|.....................|");
	} else {
	  RSprintf("\n|.....................|");
	}
      }
    }
    if (finalize){
      while(true){
        if ((i++) % op_focei.printNcol == 0){
          if (op_focei.useColor) RSprintf("\033[0m");
          RSprintf("\n");
          break;
        } else {
          RSprintf("...........|");
        }
      }
    } else {
      RSprintf("\n");
    }
  }
  return ret;
}

//[[Rcpp::export]]
double foceiOuterF(NumericVector &theta){
  int n = theta.size();
  void *ex = NULL;
  return foceiOfvOptim(n, theta.begin(), ex);
}

extern "C" void outerGradNumOptim(int n, double *par, double *gr, void *ex){
  numericGrad(par, gr);
  op_focei.nG++;
  int finalize=0, i = 0;
  niterGrad.push_back(niter.back());
  if (op_focei.derivMethod == 0){
    if (op_focei.curGill){
      gradType.push_back(1);
    } else if (op_focei.mixDeriv){
      gradType.push_back(2);
    } else{
      gradType.push_back(3);
    }
  } else {
    gradType.push_back(4);
  }
  if (op_focei.printOuter != 0 && op_focei.nG % op_focei.printOuter == 0){
    if (op_focei.useColor && op_focei.printNcol >= n){
      switch(gradType.back()){
      case 1:
	RSprintf("|\033[4m    G|    Gill Diff. |");
	break;
      case 2:
	RSprintf("|\033[4m    M|   Mixed Diff. |");
	break;
      case 3:
	RSprintf("|\033[4m    F| Forward Diff. |");
	break;
      case 4:
	RSprintf("|\033[4m    C| Central Diff. |");
	break;
      }
    } else {
      switch(gradType.back()){
      case 1:
	RSprintf("|    G|    Gill Diff. |");
	break;
      case 2:
	RSprintf("|    M|   Mixed Diff. |");
	break;
      case 3:
	RSprintf("|    F| Forward Diff. |");
	break;
      case 4:
	RSprintf("|    C| Central Diff. |");
	break;
      }
    }
    for (i = 0; i < n; i++){
      RSprintf("%#10.4g ", gr[i]);
      if (op_focei.useColor && op_focei.printNcol >= n && i == n-1){
	RSprintf("\033[0m");
      }
      RSprintf("|");
      if ((i + 1) != n && (i + 1) % op_focei.printNcol == 0){
        if (op_focei.useColor && op_focei.printNcol + i  >= op_focei.npars){
          RSprintf("\n\033[4m|.....................|");
        } else {
          RSprintf("\n|.....................|");
        }
        finalize=1;
      }
    }
    if (finalize){
      while(true){
        if ((i++) % op_focei.printNcol == 0){
          if (op_focei.useColor) RSprintf("\033[0m");
          RSprintf("\n");
	  break;
        } else {
          RSprintf("...........|");
	}
      }
    } else {
      RSprintf("\n");
    }
    if (!op_focei.useColor){
      foceiPrintLine(min2(op_focei.npars, op_focei.printNcol));
    }
  }
  vGrad.push_back(NA_REAL); // Gradient doesn't record objf
  for (i = 0; i < n; i++){
    if (gr[i] == 0){
      if (op_focei.nF+op_focei.nF2 == 1){
	stop("On initial gradient evaluation, one or more parameters have a zero gradient\nChange model, try different initial estimates or use outerOpt=\"bobyqa\")");
      } else {
	op_focei.zeroGrad=true;
	gr[i]=sqrt(DBL_EPSILON);
      }
    }
    vGrad.push_back(gr[i]);
  }
}

//[[Rcpp::export]]
NumericVector foceiOuterG(NumericVector &theta){
  int n = theta.size();
  void *ex = NULL;
  NumericVector gr(n);
  outerGradNumOptim(n, theta.begin(), gr.begin(), ex);
  return gr;
}

void foceiLbfgsb3(Environment e){
  void *ex = NULL;
  double Fmin;
  int fail, fncount=0, grcount=0;
  NumericVector x(op_focei.npars);
  NumericVector g(op_focei.npars);
  for (unsigned int k = op_focei.npars; k--;){
    x[k]=scalePar(op_focei.initPar, k);
  }
  char msg[100];
  lbfgsb3C(op_focei.npars, op_focei.lmm, x.begin(), op_focei.lower,
           op_focei.upper, op_focei.nbd, &Fmin, foceiOfvOptim,
           outerGradNumOptim, &fail, ex, op_focei.factr,
           op_focei.pgtol, &fncount, &grcount,
           op_focei.maxOuterIterations, msg, 0, -1,
	   op_focei.abstol, op_focei.reltol, g.begin());
  // Recalculate OFV in case the last calculated OFV isn't at the minimum....
  // Otherwise ETAs may be off
  std::fill_n(&op_focei.goldEta[0], op_focei.gEtaGTransN, -42.0); // All etas = -42;  Unlikely if normal
  // Finalize environment
  foceiOuterFinal(x.begin(), e);
  e["convergence"] = fail;
  e["message"] = msg;
  e["lastGrad"] = g;
}

void foceiLbfgsb(Environment e){
  void *ex = NULL;
  double Fmin;
  int fail, fncount=0, grcount=0;
  NumericVector x(op_focei.npars);
  for (unsigned int k = op_focei.npars; k--;){
    x[k]=scalePar(op_focei.initPar, k);
  }
  char msg[100];
  lbfgsbRX(op_focei.npars, op_focei.lmm, x.begin(), op_focei.lower,
           op_focei.upper, op_focei.nbd, &Fmin, foceiOfvOptim,
           outerGradNumOptim, &fail, ex, op_focei.factr,
           op_focei.pgtol, &fncount, &grcount,
           op_focei.maxOuterIterations, msg, 0, op_focei.maxOuterIterations+1);
  // Recalculate OFV in case the last calculated OFV isn't at the minimum....
  // Otherwise ETAs may be off
  std::fill_n(&op_focei.goldEta[0], op_focei.gEtaGTransN, -42.0); // All etas = -42;  Unlikely if normal
  // Finalize environment
  foceiOuterFinal(x.begin(), e);
  e["convergence"] = fail;
  e["message"] = msg;
}

void foceiCustomFun(Environment e){
  NumericVector x(op_focei.npars);
  NumericVector lower(op_focei.npars);
  NumericVector upper(op_focei.npars);
  for (unsigned int k = op_focei.npars; k--;){
    x[k]=scalePar(op_focei.initPar, k);
  }
  std::copy(&op_focei.upper[0], &op_focei.upper[0]+op_focei.npars, &upper[0]);
  std::copy(&op_focei.lower[0], &op_focei.lower[0]+op_focei.npars, &lower[0]);
  Function loadNamespace("loadNamespace", R_BaseNamespace);
  Environment nlmixr = loadNamespace("nlmixr");
  Function f = as<Function>(nlmixr["foceiOuterF"]);
  Function g = as<Function>(nlmixr["foceiOuterG"]);
  List ctl = e["control"];
  Function opt = as<Function>(ctl["outerOptFun"]);
  //.bobyqa <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...)
  List ret = as<List>(opt(_["par"]=x, _["fn"]=f, _["gr"]=g, _["lower"]=lower,
			  _["upper"]=upper,_["control"]=ctl));
  x = ret["x"];
  // Recalculate OFV in case the last calculated OFV isn't at the minimum....
  // Otherwise ETAs may be off
  if (op_focei.neta != 0) std::fill_n(&op_focei.goldEta[0], op_focei.gEtaGTransN, -42.0); // All etas = -42;  Unlikely if normal
  // Finalize environment
  foceiOuterFinal(x.begin(), e);
  e["convergence"] = ret["convergence"];
  e["message"] = ret["message"];
  e["optReturn"] = ret;
}


////////////////////////////////////////////////////////////////////////////////
// Overall Outer Problem

//[[Rcpp::export]]
Environment foceiOuter(Environment e){
  op_focei.nF=0;
  op_focei.nG=0;
  if (op_focei.maxOuterIterations > 0){
    for (unsigned int k = op_focei.npars; k--;){
      if (R_FINITE(op_focei.lower[k])){
	op_focei.lower[k]=scalePar(op_focei.lower, k);
      }
      if (R_FINITE(op_focei.upper[k])) {
	op_focei.upper[k]=scalePar(op_focei.upper,k);
      }
    }

    switch(op_focei.outerOpt){
    case 0:
      foceiLbfgsb(e);
      break;
    case 1:
      foceiLbfgsb3(e);
      break;
    case -1:
      foceiCustomFun(e);
    }
  } else {
    NumericVector x(op_focei.npars);
    for (unsigned int k = op_focei.npars; k--;){
      x[k]=scalePar(op_focei.initPar, k);
    }
    foceiOuterFinal(x.begin(), e);
    if (op_focei.maxInnerIterations == 0){
      e["fail"] = NA_INTEGER;
      e["message"] = "Likelihood evaluation with provided ETAs";
    } else {
      e["fail"] = 0;
      e["message"] = "Posthoc prediction with provided THETAs";
    }
  }
  return e;
}

//[[Rcpp::export]]
List nlmixrGill83_(Function what, NumericVector args, Environment envir,
		   LogicalVector which,
		   double gillRtol, int gillK=10, double gillStep=2, double gillFtol=0, bool optGillF=true){
  if (args.size()!=which.size()) stop("'args' must have same size as 'which'");
  gillRfn_=what;
  gillThetaN=args.size();
  gillRfnE_=envir;
  double *theta;
  theta = &args[0];
  foceiGill=0;
  NumericVector hfN(args.size());
  NumericVector hphifN(args.size());
  NumericVector gillDfN(args.size());
  NumericVector gillDf2N(args.size());
  NumericVector gillErrN(args.size());
  NumericVector aEps(args.size());
  NumericVector rEps(args.size());
  NumericVector aEpsC(args.size());
  NumericVector rEpsC(args.size());
  IntegerVector retN(args.size());
  gillLong=false;
  NumericVector fN(args.size());
  for (int i = args.size(); i--;){
    if (which[i]){
      gillPar=i;
      if (i == args.size()-1 || gillLong){
	gillF = gillRfn(theta);
      }
      fN[i] = gillF;
      retN[i] = gill83(&hfN[i], &hphifN[i], &gillDfN[i], &gillDf2N[i], &gillErrN[i],
		       theta, i, gillRtol, gillK, gillStep,
		       gillFtol,
		       -1, gill83fnG) + 1;
      double err=1/(std::fabs(theta[i])+1);
      aEps[i]  = hfN[i]*err;
      rEps[i]  = hfN[i]*err;
      if(optGillF){
	aEpsC[i] = hfN[i]*err;
	rEpsC[i] = hfN[i]*err;
      } else {
	aEpsC[i] = hphifN[i]*err;
	rEpsC[i] = hphifN[i]*err;
      }
    } else {
      retN[i] = 1;
      hfN[i] = NA_REAL;
      hphifN[i] = NA_REAL;
      gillDfN[i] = NA_REAL;
      gillDf2N[i] = NA_REAL;
      gillErrN[i] = NA_REAL;
      aEps[i]  = NA_REAL;
      rEps[i]  = NA_REAL;
      aEpsC[i]  = NA_REAL;
      rEpsC[i]  = NA_REAL;
    }
  }
  foceiGill=1;
  List df(11);
  retN.attr("levels") = CharacterVector::create("Not Assessed","Good","High Grad Error",
						"Constant Grad","Odd/Linear Grad",
						"Grad changes quickly");
  retN.attr("class") = "factor";
  df[0] = retN;
  df[1] = hfN;
  df[2] = hphifN;
  df[3] = gillDfN;
  df[4] = gillDf2N;
  df[5] = gillErrN;
  df[6] = aEps;
  df[7] = rEps;
  df[8] = aEpsC;
  df[9] = rEpsC;
  df[10] = fN;
  df.attr("names") = CharacterVector::create("info","hf","hphi","df","df2","err","aEps","rEps","aEpsC","rEpsC","f");
  if (args.hasAttribute("names")){
    df.attr("row.names") = args.attr("names");
  } else {
    df.attr("row.names") = IntegerVector::create(NA_INTEGER, -args.size());
  }
  CharacterVector cls = CharacterVector::create("nlmixrGill83", "data.frame");
  List info = List::create(_["which"]=which,
			   _["gillRtol"]=gillRtol,
			   _["gillK"]=gillK,
			   _["gillStep"]=gillStep,
			   _["gillFtol"]=gillFtol,
			   _["optGillF"]=optGillF);
  info.attr("class") = CharacterVector::create("nlmixrLstSilent");
  cls.attr(".nlmixrGill") = info;
  df.attr("class") = cls;
  return df;
}
//' @rdname nlmixrGradFun
//' @export
//[[Rcpp::export]]
double nlmixrEval_(NumericVector theta, std::string md5){
  Function loadNamespace("loadNamespace", R_BaseNamespace);
  Environment nlmixr = loadNamespace("nlmixr");
  Environment gradInfo = nlmixr[".nlmixrGradInfo"];
  std::string EF = md5 + ".f";
  std::string EE = md5 + ".e";
  std::string EW = md5 + ".w";
  LogicalVector lEW;
  if (!gradInfo.exists(EW)){
    LogicalVector tmp(theta.size());
    for (int i = theta.size(); i--;){
      tmp[i] = true;
    }
    gradInfo[EW] = tmp;
  }
  lEW=gradInfo[EW];
  if (lEW.size() != theta.size()) stop("invalid theta size");
  Function cFun = as<Function>(gradInfo[EF]);
  Environment cEnvir = as<Environment>(gradInfo[EE]);
  double f0;
  List par(1);
  par[0] = theta;
  f0 = as<double>(doCall(_["what"] = cFun, _["args"]=par, _["envir"]=cEnvir));
  std::string f0s = md5 + ".fc";
  std::string f0t = md5 + ".ft";
  std::string cns = md5 + ".n";
  f0 = as<double>(doCall(_["what"] = cFun, _["args"]=par, _["envir"]=cEnvir));
  gradInfo[f0s] = f0;
  gradInfo[f0t] = theta;
  int cn = gradInfo[cns]; cn++;
  gradInfo[cns] = cn;
  bool useColor = as<bool>(gradInfo["useColor"]);
  int printNcol=as<int>(gradInfo["printNcol"]);
  int printN=as<int>(gradInfo["print"]);
  int i, finalize=0, n=theta.size();
  bool isRstudio=as<bool>(gradInfo["isRstudio"]);
  if (cn == 1){
    vGrad.clear();
    vPar.clear();
    iterType.clear();
    gradType.clear();
    niter.clear();
    niterGrad.clear();
    if (printN != 0){
      foceiPrintLine(min2(n, printNcol));
      if (gradInfo.exists("thetaNames")){
	CharacterVector tn;
	tn = gradInfo["thetaNames"];
	if (tn.size()!=lEW.size()){
	  CharacterVector tn2(lEW.size());
	  for (int i = 0; i < lEW.size(); i++){
	    tn2[i] = "t" + std::to_string(i+1);
	  }
	  gradInfo["thetaNames"]=tn2;
	}
      } else {
	CharacterVector tn(lEW.size());
	for (int i = 0; i < lEW.size(); i++){
	  tn[i] = "t" + std::to_string(i+1);
	}
	gradInfo["thetaNames"]=tn;
      }
      CharacterVector thetaNames = gradInfo["thetaNames"];
      RSprintf("|    #| Objective Fun |");
      int i=0, finalize=0;
      std::string tmpS;
      for (i = 0; i < n; i++){
	tmpS = thetaNames[i];
	RSprintf("%#10s |", tmpS.c_str());
	if ((i + 1) != n && (i + 1) % printNcol == 0){
	  if (useColor && printNcol + i  >= n){
	    RSprintf("\n\033[4m|.....................|");
	  } else {
	    RSprintf("\n|.....................|");
	  }
	  finalize=1;
	}
      }
      if (finalize){
	while(true){
	  if ((i++) % printNcol == 0){
	    if (useColor) RSprintf("\033[0m");
	    RSprintf("\n");
	    break;
	  } else {
	    RSprintf("...........|");
	  }
	}
      } else {
	RSprintf("\n");
      }
    }
  }
  bool doUnscaled = false;
  std::string unscaledPar = md5 + ".uPar";
  NumericVector thetaU;
  niter.push_back(cn);
  // Scaled
  vPar.push_back(f0);
  if (gradInfo.exists(unscaledPar)){
    thetaU=as<NumericVector>(gradInfo[unscaledPar]);
    if (thetaU.size() != theta.size()){
      iterType.push_back(6);
    } else {
      doUnscaled=true;
      iterType.push_back(5);
    }
  } else {
    // Actually unscaled
    iterType.push_back(6);
  }
  for (i = 0; i < n; i++){
    vPar.push_back(theta[i]);
  }
  if (printN != 0 && cn % printN == 0){
    if (useColor && isRstudio)
      RSprintf("|\033[1m%5d\033[0m|%#14.8g |", cn, f0);
    else
      RSprintf("|%5d|%#14.8g |", cn, f0);
    for (i = 0; i < n; i++){
      RSprintf("%#10.4g |", theta[i]);
      if ((i + 1) != n && (i + 1) % printNcol == 0){
        if (useColor && printNcol + i  > n){
          RSprintf("\n\033[4m|.....................|");
        } else {
          RSprintf("\n|.....................|");
        }
	finalize=1;
      }
    }
    if (finalize){
      while(true){
        if ((i++) % printNcol == 0){
          if (useColor) RSprintf("\033[0m");
          RSprintf("\n");
          break;
        } else {
          RSprintf("...........|");
        }
      }
    } else {
      RSprintf("\n");
    }
  }
  if (doUnscaled){
    iterType.push_back(6);
    niter.push_back(niter.back());
    finalize=0;
    // No obj scaling currently
    vPar.push_back(f0);
    for (i = 0; i < n; i++){
      vPar.push_back(thetaU[i]);
    }
    if (printN != 0 && cn % printN == 0){
      if (useColor && isRstudio)
	RSprintf("|    U|%#14.8g |", f0);
      else
	RSprintf("|    U|%#14.8g |", f0);
      for (i = 0; i < n; i++){
	RSprintf("%#10.4g |", thetaU[i]);
	if ((i + 1) != n && (i + 1) % printNcol == 0){
	  if (useColor && printNcol + i  > n){
	    RSprintf("\n\033[4m|.....................|");
	  } else {
	    RSprintf("\n|.....................|");
	  }
	  finalize=1;
	}
      }
      if (finalize){
	while(true){
	  if ((i++) % printNcol == 0){
	    if (useColor) RSprintf("\033[0m");
	    RSprintf("\n");
	    break;
	  } else {
	    RSprintf("...........|");
	  }
	}
      } else {
	RSprintf("\n");
      }
    }
  }
  return f0;
}

void nlmixrGradPrint(NumericVector gr, int gradType, int cn, bool useColor,
		     int printNcol, int printN, bool isRstudio){
  int n = gr.size(), finalize=0, i;
  if (printN != 0 && cn % printN == 0){
    if (useColor && printNcol >= n){
      switch(gradType){
      case 1:
	RSprintf("|\033[4m    G|    Gill Diff. |");
	break;
      case 2:
	RSprintf("|\033[4m    M|   Mixed Diff. |");
	break;
      case 3:
	RSprintf("|\033[4m    F| Forward Diff. |");
	break;
      case 4:
	RSprintf("|\033[4m    C| Central Diff. |");
	break;
      }
    } else {
      switch(gradType){
      case 1:
	RSprintf("|    G|    Gill Diff. |");
	break;
      case 2:
	RSprintf("|    M|   Mixed Diff. |");
	break;
      case 3:
	RSprintf("|    F| Forward Diff. |");
	break;
      case 4:
	RSprintf("|    C| Central Diff. |");
	break;
      }
    }
    for (i = 0; i < n; i++){
      RSprintf("%#10.4g ", gr[i]);
      if (useColor && printNcol >= n && i == n-1){
	RSprintf("\033[0m");
      }
      RSprintf("|");
      if ((i + 1) != n && (i + 1) % printNcol == 0){
        if (useColor && printNcol + i  >= n){
          RSprintf("\n\033[4m|.....................|");
        } else {
          RSprintf("\n|.....................|");
        }
        finalize=1;
      }
    }
    if (finalize){
      while(true){
        if ((i++) % printNcol == 0){
          if (useColor) RSprintf("\033[0m");
          RSprintf("\n");
	  break;
        } else {
          RSprintf("...........|");
	}
      }
    } else {
      RSprintf("\n");
    }
    if (!useColor){
      foceiPrintLine(min2(n, printNcol));
    }
  }
}

//' @rdname nlmixrGradFun
//' @export
//[[Rcpp::export]]
RObject nlmixrUnscaled_(NumericVector theta, std::string md5){
  // Unscaled
  Function loadNamespace("loadNamespace", R_BaseNamespace);
  Environment nlmixr = loadNamespace("nlmixr");
  Environment gradInfo = nlmixr[".nlmixrGradInfo"];
  std::string unscaledPar = md5 + ".uPar";
  gradInfo[unscaledPar] = theta;
  return R_NilValue;
}

//' @rdname nlmixrGradFun
//' @export
//[[Rcpp::export]]
NumericVector nlmixrGrad_(NumericVector theta, std::string md5){
  Function loadNamespace("loadNamespace", R_BaseNamespace);
  Environment nlmixr = loadNamespace("nlmixr");
  Environment gradInfo = nlmixr[".nlmixrGradInfo"];

  std::string Egill = md5 + ".g";
  std::string EF = md5 + ".f";
  std::string EE = md5 + ".e";

  Function cFun = as<Function>(gradInfo[EF]);
  Environment cEnvir = as<Environment>(gradInfo[EE]);
  List Lgill;
  std::string EW = md5 + ".w";
  LogicalVector lEW;
  bool useColor = as<bool>(gradInfo["useColor"]);
  int printNcol=as<int>(gradInfo["printNcol"]);
  int printN=as<int>(gradInfo["print"]);
  bool isRstudio=as<bool>(gradInfo["isRstudio"]);
  if (gradInfo.exists(Egill)){
    lEW=gradInfo[EW];
    if (lEW.size() != theta.size()){
      stop("Invalid theta size");
    }
    Lgill = gradInfo[Egill];
  } else {
    if (!gradInfo.exists(EW)){
      LogicalVector tmp(theta.size());
      for (int i = theta.size(); i--;){
	tmp[i] = true;
      }
      gradInfo[EW] = tmp;
    }
    lEW=gradInfo[EW];
    if (lEW.size() != theta.size()){
      stop("Invalid theta size (or which size)");
    }
    std::string Ertol = md5 + ".rtol";
    std::string EK = md5 + ".k";
    std::string Estep = md5 + ".s";
    std::string EFtol = md5 + ".ftol";
    Lgill = nlmixrGill83_(cFun, theta, gradInfo[EE],
			  lEW, gradInfo[Ertol],
			  gradInfo[EK], gradInfo[Estep],
			  gradInfo[EFtol]);
    gradInfo[Egill]=Lgill;
    niterGrad.push_back(niter.back());
    gradType.push_back(1);
    vGrad.push_back(NA_REAL); // Gradient doesn't record objf
    NumericVector gr = as<NumericVector>(Lgill["df"]);
    for (int i = 0; i < gr.size(); i++){
      if (gr[i] == 0){
	stop("On initial gradient evaluation, one or more parameters have a zero gradient\nChange model, try different initial estimates or try derivative free optimization)");
      }
      vGrad.push_back(gr[i]);
    }
    nlmixrGradPrint(gr, gradType.back(), niter.back(), useColor,
		    printNcol, printN, isRstudio);
    return gr;
  }
  NumericVector aEps = as<NumericVector>(Lgill["aEps"]);
  NumericVector rEps = as<NumericVector>(Lgill["rEps"]);
  NumericVector aEpsC = as<NumericVector>(Lgill["aEpsC"]);
  NumericVector rEpsC = as<NumericVector>(Lgill["rEpsC"]);
  NumericVector g(theta.size());
  double f0, delta, cur, tmp=0, tmp0;
  bool doForward=true;
  // FIXME
  List par(1);
  par[0] = theta;
  std::string f0s = md5 + ".fc";
  std::string f0t = md5 + ".ft";
  bool reEval = true;
  if (gradInfo.exists(f0s)){
    NumericVector thetaL = gradInfo[f0t];
    if (thetaL.size() == theta.size()){
      reEval=false;
      for (int i = theta.size(); i--;){
	if (thetaL[i] != theta[i]){
	  reEval=true;
	  break;
	}
      }
      if (!reEval){
	NumericVector tmp  = gradInfo[f0s];
	if (tmp.size() == 1){
	  f0 = tmp[0];
	} else {
	  reEval=true;
	}
      }
    }
  }
  if (reEval){
    f0 = as<double>(doCall(_["what"] = cFun, _["args"]=par, _["envir"]=cEnvir));
  }
  niterGrad.push_back(niter.back());
  vGrad.push_back(NA_REAL); // Gradient doesn't record objf
  bool isMixed=false;
  for (int i = theta.size(); i--;){
    cur = theta[i];
    if (doForward){
      delta = (std::fabs(theta[i])*rEps[i] + aEps[i]);
      theta[i] = cur + delta;
      par[0] = theta;
      tmp = as<double>(doCall(_["what"] = cFun, _["args"]=par, _["envir"]=cEnvir));
      g[i] = (tmp-f0)/delta;
      theta[i] = cur;
    } else {
      delta = (std::fabs(theta[i])*rEpsC[i] + aEpsC[i]);
      theta[i] = cur + delta;
      par[0] = theta;
      tmp0 = as<double>(doCall(_["what"] = cFun, _["args"]=par, _["envir"]=cEnvir));
      theta[i] = cur - delta;
      par[0] = theta;
      tmp = as<double>(doCall(_["what"] = cFun, _["args"]=par, _["envir"]=cEnvir));
      g[i] = (tmp0-tmp)/(2*delta);
      theta[i] = cur;
    }

    // Check for bad grad
    if (std::isnan(g[i]) ||  ISNA(g[i]) || !R_FINITE(g[i])){
      if (doForward){
      	// Switch to Backward difference method
      	// op_focei.mixDeriv=1;
      	theta[i] = cur - delta;
	par[0] = theta;
	tmp0 = as<double>(doCall(_["what"] = cFun, _["args"]=par, _["envir"]=cEnvir));
      	g[i] = (f0-tmp0)/(delta);
	isMixed=true;
      } else {
      	// We are using the central difference AND there is an NA in one of the terms
      	// g[cpar] = (tmp0-tmp)/(2*delta);
      	// op_focei.mixDeriv=1;
	isMixed=true;
      	if (std::isnan(tmp0) || ISNA(tmp0) || !R_FINITE(tmp0)){
      	  // Backward
      	  g[i] = (f0-tmp)/delta;
      	} else {
      	  // Forward
      	  g[i] = (tmp0-f0)/delta;
      	}
      }
    }
  }
  for (int i = 0; i < theta.size(); i++){
    vGrad.push_back(g[i]);
  }
  if (isMixed){
    gradType.push_back(2);
  } else if (doForward) {
    gradType.push_back(3);
  } else {
    gradType.push_back(4);
  }
  nlmixrGradPrint(g, gradType.back(), niter.back(), useColor,
		  printNcol, printN, isRstudio);
  return g;
}
//' @rdname nlmixrGradFun
//' @export
//[[Rcpp::export]]
RObject nlmixrParHist_(std::string md5){
  Function loadNamespace("loadNamespace", R_BaseNamespace);
  Environment nlmixr = loadNamespace("nlmixr");
  Environment gradInfo = nlmixr[".nlmixrGradInfo"];
  std::string EW = md5 + ".w";
  LogicalVector lEW;
  if (!gradInfo.exists(EW)){
    LogicalVector tmp(lEW.size());
    for (int i = lEW.size(); i--;){
      tmp[i] = true;
    }
    gradInfo[EW] = tmp;
  }
  lEW=gradInfo[EW];
  if (gradInfo.exists("thetaNames")){
    CharacterVector tn;
    tn = gradInfo["thetaNames"];
    if (tn.size()!=lEW.size()){
      CharacterVector tn2(lEW.size());
      for (int i = 0; i < lEW.size(); i++){
	tn2[i] = "t" + std::to_string(i+1);
      }
      gradInfo["thetaNames"]=tn2;
    }
  } else {
    CharacterVector tn(lEW.size());
    for (int i = 0; i < lEW.size(); i++){
      tn[i] = "t" + std::to_string(i+1);
    }
    gradInfo["thetaNames"]=tn;
  }
  std::string cns = md5 + ".n";
  gradInfo[cns] = 0;
  parHistData(gradInfo, false);
  return gradInfo["parHistData"];
}

//[[Rcpp::export]]
RObject nlmixrHess_(RObject thetaT, RObject fT, RObject e,
		    RObject gillInfoT){
  par_progress = (par_progress_t) R_GetCCallable("RxODE", "par_progress");
  List par(1);
  NumericVector theta = as<NumericVector>(thetaT);
  Function f = as<Function>(fT);
  List gillInfo = as<List>(gillInfoT);
  arma::mat H(theta.size(), theta.size(), fill::zeros);
  double epsI, epsJ;
  NumericVector rEpsC = as<NumericVector>(gillInfo["rEpsC"]);
  NumericVector aEpsC = as<NumericVector>(gillInfo["aEpsC"]);
  NumericVector nF = as<NumericVector>(gillInfo["f"]);
  double lastOfv=nF[0];
  int n = theta.size();
  double f1,f2,f3,f4;
  double ti, tj;
  int i, j;
  int totTick= 4*n +2*n*(n-1);
  int cur = 0, curTick=0;
  clock_t t0=clock();
  for (i=n; i--;){
    epsI = (std::fabs(theta[i])*rEpsC[i] + aEpsC[i]);
    ti = theta[i];
    theta[i] = ti + 2*epsI;
    par[0]=theta;
    f1 = as<double>(doCall(_["what"] = f, _["args"]=par, _["envir"]=e));
    cur++;
    curTick = par_progress(cur, totTick, curTick, 1, t0, 0);
    theta[i] = ti + epsI;
    par[0]=theta;
    f2 = as<double>(doCall(_["what"] = f, _["args"]=par, _["envir"]=e));
    cur++;
    curTick = par_progress(cur, totTick, curTick, 1, t0, 0);
    theta[i] = ti - epsI;
    par[0]=theta;
    f3 = as<double>(doCall(_["what"] = f, _["args"]=par, _["envir"]=e));
    cur++;
    curTick = par_progress(cur, totTick, curTick, 1, t0, 0);
    theta[i] = ti - 2*epsI;
    par[0]=theta;
    f4 = as<double>(doCall(_["what"] = f, _["args"]=par, _["envir"]=e));
    cur++;
    curTick = par_progress(cur, totTick, curTick, 1, t0, 0);
    theta[i] = ti;
    H(i,i)=(-f1+16*f2-30*lastOfv+16*f3-f4)/(12*epsI*epsI);
    for (j = i; j--;){
      epsJ = (std::fabs(theta[j])*rEpsC[j] + aEpsC[j]);
      // eps = sqrt(epsI*epsJ);// 0.5*epsI+0.5*epsJ;
      // epsI = eps;
      // epsJ = eps;
      tj = theta[j];
      theta[i] = ti + epsI;
      theta[j] = tj + epsJ;
      par[0]=theta;
      f1 = as<double>(doCall(_["what"] = f, _["args"]=par, _["envir"]=e));
      cur++;
      curTick = par_progress(cur, totTick, curTick, 1, t0, 0);
      theta[i] = ti + epsI;
      theta[j] = tj - epsJ;
      par[0]=theta;
      f2 = as<double>(doCall(_["what"] = f, _["args"]=par, _["envir"]=e));
      cur++;
      curTick = par_progress(cur, totTick, curTick, 1, t0, 0);
      theta[i] = ti - epsI;
      theta[j] = tj + epsJ;
      par[0]=theta;
      f3 = as<double>(doCall(_["what"] = f, _["args"]=par, _["envir"]=e));
      cur++;
      curTick = par_progress(cur, totTick, curTick, 1, t0, 0);
      theta[i] = ti - epsI;
      theta[j] = tj - epsJ;
      par[0]=theta;
      f4 = as<double>(doCall(_["what"] = f, _["args"]=par, _["envir"]=e));
      cur++;
      curTick = par_progress(cur, totTick, curTick, 1, t0, 0);
      H(i,j)= (f1-f2-f3+f4)/(4*epsI*epsJ);
      H(j,i) = H(i,j);
      theta[i] = ti;
      theta[j] = tj;
    }
  }
  par_progress(totTick, totTick, cur, 1, t0, 0);
  if (isRstudio){
      RSprintf("\n");
  } else {
      RSprintf("\r                                                                                \r");
  }
  return wrap(H);
}

////////////////////////////////////////////////////////////////////////////////
// Covariance functions
int foceiCalcR(Environment e){
  rx = getRx();
  arma::mat H(op_focei.npars, op_focei.npars);
  arma::vec theta(op_focei.npars);
  unsigned int i, j, k;
  for (k = op_focei.npars; k--;){
    j=op_focei.fixedTrans[k];
    theta[k] = op_focei.fullTheta[j];
  }

  // arma::vec df1(op_focei.npars);
  // arma::vec df2(op_focei.npars);
  double epsI, epsJ;

  bool doForward=false;
  if (op_focei.derivMethod == 0){
    doForward=true;
  }
  double f1,f2,f3,f4;
  double ti, tj;
  if (doForward){
    stop("Not implemented for finite differences.");
  } else {
    //https://pdfs.semanticscholar.org/presentation/e7d5/aff49eb17fd155e75725c295859d983cfda4.pdf
    // https://v8doc.sas.com/sashtml/ormp/chap5/sect28.htm
    double fnscale = 1.0;
    if (op_focei.scaleObjective == 2){
      fnscale = op_focei.initObjective / op_focei.scaleObjectiveTo;
    }
    double parScaleI=1.0, parScaleJ=1.0;
    for (i=op_focei.npars; i--;){
      epsI = (std::fabs(theta[i])*op_focei.rEpsC[i] + op_focei.aEpsC[i]);
      ti = theta[i];
      theta[i] = ti + 2*epsI;
      updateTheta(theta.begin());
      f1 = foceiOfv0(theta.begin());
      if (ISNA(f1)) return 0;
      op_focei.cur++;
      op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
      theta[i] = ti + epsI;
      updateTheta(theta.begin());
      f2 = foceiOfv0(theta.begin());
      if (ISNA(f2)) return 0;
      op_focei.cur++;
      op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
      theta[i] = ti - epsI;
      updateTheta(theta.begin());
      f3 = foceiOfv0(theta.begin());
      if (ISNA(f3)) return 0;
      op_focei.cur++;
      op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
      theta[i] = ti - 2*epsI;
      updateTheta(theta.begin());
      f4 = foceiOfv0(theta.begin());
      if (ISNA(f4)) return 0;
      op_focei.cur++;
      op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
      theta[i] = ti;
      // RSprintf("-- i:%d, i: %d\n", i, i);
      // print(NumericVector::create(f1,f2,f3,f4,op_focei.lastOfv));
      H(i,i)=fnscale*(-f1+16*f2-30*op_focei.lastOfv+16*f3-f4)/(12*epsI*epsI*parScaleI*parScaleI);
      for (j = i; j--;){
	epsJ = (std::fabs(theta[j])*op_focei.rEpsC[j] + op_focei.aEpsC[j]);
	// eps = sqrt(epsI*epsJ);// 0.5*epsI+0.5*epsJ;
	// epsI = eps;
	// epsJ = eps;
	tj = theta[j];
	theta[i] = ti + epsI;
	theta[j] = tj + epsJ;
	updateTheta(theta.begin());
	f1 = foceiOfv0(theta.begin());
	if (ISNA(f1)) return 0;
	op_focei.cur++;
	op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
	theta[i] = ti + epsI;
	theta[j] = tj - epsJ;
	updateTheta(theta.begin());
	f2 = foceiOfv0(theta.begin());
	if (ISNA(f2)) return 0;
	op_focei.cur++;
	op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
	theta[i] = ti - epsI;
	theta[j] = tj + epsJ;
	updateTheta(theta.begin());
	f3 = foceiOfv0(theta.begin());
	if (ISNA(f3)) return 0;
	op_focei.cur++;
	op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
	theta[i] = ti - epsI;
	theta[j] = tj - epsJ;
	updateTheta(theta.begin());
	f4 = foceiOfv0(theta.begin());
	if (ISNA(f4)) return 0;
	op_focei.cur++;
	op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
	// RSprintf("-- i:%d, j: %d\n", i, j);
	// print(NumericVector::create(f1,f2,f3,f4));
	H(i,j)= fnscale*(f1-f2-f3+f4)/(4*epsI*epsJ*parScaleI*parScaleJ);
	H(j,i) = H(i,j);
	theta[i] = ti;
	theta[j] = tj;
      }
    }
  }
  // R matrix = Hessian/2
  H = H*0.5;
  // https://github.com/cran/nmw/blob/59478fcc91f368bb3bbc23e55d8d1d5d53726a4b/R/CovStep.R
  // H = 0.25*H + 0.25*H.t();
  if (e.exists("R.1")){
    // This is the 2nd attempt
    arma::mat H2 = 0.5*H + 0.5*as<arma::mat>(e["R.1"]);
    arma::mat cholR;
    arma::mat RE;
    bool rpd = cholSE0(cholR, RE, H2, op_focei.cholSEtol);
    if (rpd){
      e["R.pd"] =  rpd;
      e["R.E"] =  wrap(RE);
      e["cholR"] = wrap(cholR);
    } else {
      e["R.pd2"] = false;
      e["R.2"] = H2;
      e["R.E2"] = wrap(RE);
      e["cholR2"] = wrap(cholR);
      e["R.pd"] = cholSE0(cholR, RE, H, op_focei.cholSEtol);
      e["R.E"] =  wrap(RE);
      e["cholR"] = wrap(cholR);
    }
  } else {
    e["R.0"] = H;
    arma::mat cholR;
    arma::mat RE;
    bool rpd = cholSE0(cholR, RE, H, op_focei.cholSEtol);
    e["R.pd"] =  wrap(rpd);
    e["R.E"] =  wrap(RE);
    e["cholR"] = wrap(cholR);
  }
  return 1;
}

// Necessary for S-matrix calculation
int foceiS(double *theta, Environment e){
  rx = getRx();
  op_focei.calcGrad=1;
  rx = getRx();
  int npars = op_focei.npars;
  int cpar, gid;
  double cur, delta;
  focei_ind *fInd;
  // Do Forward difference if the OBJF for *theta has already been calculated.
  bool doForward=false;
  if (op_focei.derivMethod == 0){
    doForward=true;
    // If the first derivative wasn't calculated, then calculate it.
    for (cpar = npars; cpar--;){
      if (theta[cpar] != op_focei.theta[cpar]){
        doForward=false;
        break;
      }
    }
    if (doForward){
      // Fill in lik0
      for (gid = rx->nsub; gid--;){
        fInd = &(inds_focei[gid]);
        op_focei.likSav[gid] = -2*fInd->lik[0];
      }
    }
  }
  for (cpar = npars; cpar--;){
    if (op_focei.smatNorm){
      if (doForward){
	delta = (std::fabs(theta[cpar])*op_focei.rEps[cpar] + op_focei.aEps[cpar])/_safe_sqrt(1+std::fabs(min2(op_focei.initObjective, op_focei.lastOfv)));
      } else {
	delta = (std::fabs(theta[cpar])*op_focei.rEpsC[cpar] + op_focei.aEpsC[cpar])/_safe_sqrt(1+std::fabs(min2(op_focei.initObjective, op_focei.lastOfv)));
      }
    } else {
      if (doForward){
	delta = std::fabs(theta[cpar])*op_focei.rEps[cpar] + op_focei.aEps[cpar];
      } else {
	delta = std::fabs(theta[cpar])*op_focei.rEpsC[cpar] + op_focei.aEpsC[cpar];
      }
    }
    if (op_focei.neta != 0) std::fill_n(&op_focei.goldEta[0], op_focei.gEtaGTransN, -42.0); // All etas = -42;  Unlikely if normal
    cur = theta[cpar];
    theta[cpar] = cur + delta;
    updateTheta(theta);
    for (gid = rx->nsub; gid--;){
      if (!innerOpt1(gid,2)) return 0;
      if (doForward){
        fInd = &(inds_focei[gid]);
        fInd->thetaGrad[cpar] = (fInd->lik[2] - op_focei.likSav[gid])/delta;
      }
    }
    if (!doForward){
      if (op_focei.neta != 0) std::fill_n(&op_focei.goldEta[0], op_focei.gEtaGTransN, -42.0);
      theta[cpar] = cur - delta;
      updateTheta(theta);
      for (gid = rx->nsub; gid--;){
        if (!innerOpt1(gid,1)) return 0;
        fInd = &(inds_focei[gid]);
        fInd->thetaGrad[cpar] = (fInd->lik[2] - fInd->lik[1])/(2*delta);
      }
    }
    theta[cpar] = cur;
    op_focei.cur++;
    op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
  }
  op_focei.calcGrad=0;
  // Now calculate S matrix
  arma::mat m1(1, op_focei.npars), S(op_focei.npars, op_focei.npars, fill::zeros), s1(1, op_focei.npars,fill::ones);
  for (gid = rx->nsub; gid--;){
    fInd = &(inds_focei[gid]);
    std::copy(&fInd->thetaGrad[0],&fInd->thetaGrad[0]+op_focei.npars,&m1[0]);
    S = S + m1.t() * m1;
  }
  // S matrix = S/4
  // According to https://github.com/cran/nmw/blob/59478fcc91f368bb3bbc23e55d8d1d5d53726a4b/R/Objs.R
  S=S*0.25;
  e["S0"] = wrap(S);
  arma::mat cholS;
  arma::mat SE;
  e["S.pd"] =  cholSE0(cholS, SE, S, op_focei.cholSEtol);
  e["S.E"] =  wrap(SE);
  e["cholS"] = wrap(cholS);
  return 1;
}
//' Return the square root of general square matrix A
//'
//' @param m Matrix to take the square root of.
//'
//' @return A square root general square matrix of m
//'
//' @export
//[[Rcpp::export]]
NumericMatrix sqrtm(NumericMatrix m){
  arma::cx_mat ret = sqrtmat(as<arma::mat>(m));
  mat im = arma::imag(ret);
  mat re = arma::real(ret);
  if (arma::any(arma::any(im,0))){
    stop("Some components of sqrtm are imaginary.");
  }
  return wrap(re);
}

//[[Rcpp::export]]
NumericMatrix foceiCalcCov(Environment e){
  std::string boundStr = "";
  CharacterVector thetaNames=as<CharacterVector>(e["thetaNames"]);
  try {
    if (op_focei.covMethod) {
      op_focei.derivMethodSwitch=0;
      // Check boundaries
      unsigned int j, k;
      double cur;
      bool boundary=false;
      rx = getRx();
      if (op_focei.neta == 0) op_focei.covMethod = 2; // Always use hessian for NLS
      for (unsigned int k = op_focei.npars; k--;){
	if (R_FINITE(op_focei.lower[k])){
	  op_focei.lower[k]=unscalePar(op_focei.lower,k);
	}
	if (R_FINITE(op_focei.upper[k])) {
	  op_focei.upper[k]=unscalePar(op_focei.upper,k);
	}
      }
      if (op_focei.boundTol > 0){
	// Subtract omegan so that Omega boundaries are not counted.
	for (k = op_focei.npars-op_focei.omegan; k--;){
	  if (op_focei.nbd[k] != 0){
	    // bounds
	    j=op_focei.fixedTrans[k];
	    cur = op_focei.fullTheta[j];
	    if (op_focei.nbd[k] == 1){
	      // Lower only
	      if ((cur-op_focei.lower[k])/cur < op_focei.boundTol) {
		boundary = true;
		boundStr += "\"" + thetaNames[j] + "\" ";
	      }
	    } else if (op_focei.nbd[k] == 2){
	      // Upper and lower
	      if ((cur-op_focei.lower[k])/cur < op_focei.boundTol) {
		boundary = true;
		boundStr += "\"" + thetaNames[j] + "\" ";
	      }
	      if ((op_focei.upper[k]-cur)/cur < op_focei.boundTol){
		boundary = true;
		boundStr += "\"" + thetaNames[j] + "\" ";
	      }
	    } else {
	      // Upper only
	      if ((op_focei.upper[k]-cur)/cur < op_focei.boundTol){
		boundary = true;
		boundStr += "\"" + thetaNames[j] + "\" ";
	      }
	    }
	  }
	}
      }
      for (unsigned int j = rx->nsub; j--;){
	focei_ind *fInd = &(inds_focei[j]);
	fInd->doChol=!(op_focei.cholSECov);
	fInd->doFD = 0;
      }
      op_focei.resetEtaSize = R_PosInf; // Dont reset ETAs
      op_focei.resetEtaSize=0; // Always reset ETAs.
      NumericVector fullT = e["fullTheta"];
      NumericVector fullT2(op_focei.ntheta);
      std::copy(fullT.begin(), fullT.begin()+fullT2.size(), fullT2.begin());
      LogicalVector skipCov(op_focei.ntheta+op_focei.omegan);//skipCovN
      if (op_focei.skipCovN == 0){
	std::fill_n(skipCov.begin(), op_focei.ntheta, false);
	std::fill_n(skipCov.begin()+op_focei.ntheta, skipCov.size() - op_focei.ntheta, true);
      } else {
	std::copy(&op_focei.skipCov[0],&op_focei.skipCov[0]+op_focei.skipCovN,skipCov.begin());
	std::fill_n(skipCov.begin()+op_focei.skipCovN,skipCov.size()-op_focei.skipCovN,true);
      }
      e["skipCov"] = skipCov;
      // Unscaled objective and parameters.
      if (op_focei.scaleObjective){
	op_focei.scaleObjective=0;
	op_focei.lastOfv = op_focei.lastOfv * op_focei.initObjective / op_focei.scaleObjectiveTo;
      }
      // foceiSetupTheta_(op_focei.mvi, fullT2, skipCov, op_focei.scaleTo, false);
      foceiSetupTheta_(op_focei.mvi, fullT2, skipCov, 0, false);
      op_focei.scaleType=10;
      if (op_focei.covMethod && !boundary) {
	rx = getRx();
	op_focei.t0 = clock();
	op_focei.totTick=0;
	op_focei.cur=0;
	op_focei.curTick=0;
	RSprintf(_("calculating covariance matrix\n"));
	// Change options to covariance options
	// op_focei.scaleObjective = 0;
	op_focei.derivMethod = op_focei.covDerivMethod;

	arma::mat Rinv;
	op_focei.totTick=1;
	if (op_focei.covMethod == 1 || op_focei.covMethod == 3){
	  op_focei.totTick+=op_focei.npars;
	}
	if (op_focei.covMethod == 1 || op_focei.covMethod == 2){
	  op_focei.totTick += 2*op_focei.npars +2*(op_focei.npars*op_focei.npars);
	}
	op_focei.totTick += op_focei.npars;
	double hf, hphif, err;
	unsigned int j, k;
	arma::vec theta(op_focei.npars);
	for (k = op_focei.npars; k--;){
	  j=op_focei.fixedTrans[k];
	  theta[k] = op_focei.fullTheta[j];
	}
	std::copy(&theta[0], &theta[0] + op_focei.npars, &op_focei.theta[0]);
	for (int cpar = op_focei.npars; cpar--;){
	  err = op_focei.rmatNorm ? 1/(std::fabs(theta[cpar])+1) : 1;
	  if (op_focei.gillKcov != 0){
	    op_focei.gillRetC[cpar] = gill83(&hf, &hphif, &op_focei.gillDf[cpar], &op_focei.gillDf2[cpar], &op_focei.gillErr[cpar],
					     &theta[0], cpar, op_focei.hessEps, op_focei.gillKcov, op_focei.gillStepCov, op_focei.gillFtolCov,
					     -1, gill83fnG);
	    // h=aEps*(|x|+1)/sqrt(1+fabs(f));
	    // h*sqrt(1+fabs(f))/(|x|+1) = aEps
	    // let err=2*sqrt(epsA/(1+f))
	    // err*(aEps+|x|rEps) = h
	    // Let aEps = rEps (could be a different ratio)
	    // h/err = aEps(1+|x|)
	    // aEps=h/err/(1+|x|)
	    //
	    op_focei.aEps[cpar]  = hf*err;
	    op_focei.rEps[cpar]  = hf*err;
	    if (op_focei.covGillF){
	      op_focei.aEpsC[cpar] = hf*err;
	      op_focei.rEpsC[cpar] = hf*err;
	    } else {
	      op_focei.aEpsC[cpar] = hphif*err;
	      op_focei.rEpsC[cpar] = hphif*err;
	    }
	  } else {
	    hf = op_focei.hessEps;
	    op_focei.aEps[cpar]  = hf*err;
	    op_focei.rEps[cpar]  = hf*err;
	    op_focei.aEpsC[cpar] = hf*err;
	    op_focei.rEpsC[cpar] = hf*err;
	  }
	  op_focei.cur++;
	  op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
	}
	op_focei.didGill+=1;

	bool isPd;
	std::string rstr = "r";
	bool checkSandwich = false;
	if (op_focei.covMethod == 1 || op_focei.covMethod == 2) {
	  // R matrix based covariance
	  arma::mat cholR;
	  if (!e.exists("cholR")){
	    foceiCalcR(e);
	  } else {
	    op_focei.cur += op_focei.npars*2;
	    op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
	  }
	  if (e.exists("cholR")) {
	    isPd = as<bool>(e["R.pd"]);
	    if (!isPd){
	      isPd = true;
	      arma::vec E = as<arma::vec>(e["R.E"]);
	      for (int j = E.size(); j--;){
		if (E[j] > op_focei.cholAccept){
		  isPd=false;
		  break;
		}
	      }
	      if (isPd){
		rstr = "r+";
		checkSandwich = true;
	      }
	    }
	    if (!isPd){
	      // Suggted by https://www.tandfonline.com/doi/pdf/10.1198/106186005X78800
	      mat H0 = as<arma::mat>(e["R.0"]);
	      H0 = H0*H0;
	      cx_mat H1;
	      bool success = sqrtmat(H1,H0);
	      if (success){
		mat im = arma::imag(H1);
		mat re = arma::real(H1);
		if (!arma::any(arma::any(im,0))){
		  success= chol(H0,re);
		  if (success){
		    e["cholR"] = wrap(H0);
		    rstr = "|r|";
		    checkSandwich = true;
		    isPd = true;
		  }
		}
	      }
	    }
	    op_focei.cur += op_focei.npars*2;
	    op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
	    if (!isPd){
	      warning(_("R matrix non-positive definite"));
	      e["R"] = wrap(e["R.0"]);
	      op_focei.covMethod = 3;
	      op_focei.cur += op_focei.npars*2;
	      op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
	    } else {
	      cholR = as<arma::mat>(e["cholR"]);
	      e["R"] = wrap(trans(cholR) * cholR);
	      if (!e.exists("Rinv")){
		bool success  = inv(Rinv, trimatu(cholR));
		if (!success){
		  warning(_("Hessian (R) matrix seems singular; Using pseudo-inverse"));
		  Rinv = pinv(trimatu(cholR));
		  checkSandwich = true;
		}
		Rinv = Rinv * Rinv.t();
		e["Rinv"] = wrap(Rinv);
	      } else {
		Rinv = as<arma::mat>(e["Rinv"]);
	      }
	      op_focei.cur++;
	      op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, rx->op->cores, op_focei.t0, 0);
	      if (!e.exists("covR")){
		e["covR"] = wrap(2*Rinv);
	      }
	      if (op_focei.covMethod == 2){
		e["cov"] = as<NumericMatrix>(e["covR"]);
	      }
	    }
	  } else {
	    RSprintf("\rR matrix calculation failed; Switch to S-matrix covariance.\n");
	    op_focei.covMethod = 3;
	    op_focei.cur += op_focei.npars*2;
	    op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
	  }
	}
	arma::mat cholS;
	int origCov = op_focei.covMethod;
	std::string sstr="s";
	if (op_focei.covMethod == 1 || op_focei.covMethod == 3) {
	  arma::vec theta(op_focei.npars);
	  unsigned int j, k;
	  for (k = op_focei.npars; k--;){
	    j=op_focei.fixedTrans[k];
	    theta[k] = op_focei.fullTheta[j];
	  }
	  if (!e.exists("cholS")){
	    foceiS(&theta[0], e);
	  } else {
	    op_focei.cur += op_focei.npars;
	    op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
	  }
	  if (e.exists("cholS")) {
	    isPd = as<bool>(e["S.pd"]);
	    if (!isPd){
	      isPd=true;
	      arma::vec E = as<arma::vec>(e["S.E"]);
	      for (int j = E.size(); j--;){
		if (E[j] > op_focei.cholAccept){
		  isPd=false;
		  break;
		}
	      }
	      if (isPd){
		sstr="s+";
		checkSandwich = true;
	      }
	    }
	    if (!isPd){
	      // Suggted by https://www.tandfonline.com/doi/pdf/10.1198/106186005X78800
	      mat H0 = as<arma::mat>(e["S0"]);
	      H0 = H0*H0;
	      cx_mat H1;
	      bool success = sqrtmat(H1,H0);
	      if (success){
		mat im = arma::imag(H1);
		mat re = arma::real(H1);
		if (!arma::any(arma::any(im,0))){
		  success= chol(H0,re);
		  if (success){
		    e["cholS"] = wrap(H0);
		    sstr = "|s|";
		    checkSandwich = true;
		    isPd = true;
		  }
		}
	      }
	    }
	    if (!isPd){
	      warning(_("S matrix non-positive definite"));
	      if (op_focei.covMethod == 1){
		e["cov"] = as<NumericMatrix>(e["covR"]);
		op_focei.covMethod = 2;
	      } else {
		warning(_("cannot calculate covariance"));
	      }
	      op_focei.cur += op_focei.npars*2;
	      op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
	    } else {
	      cholS = as<arma::mat>(e["cholS"]);
	      arma::mat S;
	      if (e.exists("S")){
		S = as<arma::mat>(e["S"]);
	      } else {
		S = trans(cholS) * cholS;
		e["S"] = wrap(S);
	      }
	      if (op_focei.covMethod == 1){
		e["covRS"] = Rinv * S *Rinv;
		arma::mat covRS = as<arma::mat>(e["covRS"]);
		if (checkSandwich){
		  mat Sinv;
		  bool success;
		  success = inv(Sinv, trimatu(cholS));
		  if (!success){
		    warning(_("S matrix seems singular; Using pseudo-inverse"));
		    Sinv = pinv(trimatu(cholS));
		  }
		  Sinv = Sinv * Sinv.t();
		  e["covS"]= 4 * Sinv;
		  if (rstr == "r"){
		    // Use covR
		    e["cov"] = as<NumericMatrix>(e["covR"]);
		    op_focei.covMethod=2;
		  } else if (sstr == "s"){
		    // use covS
		    e["cov"] = as<NumericMatrix>(e["covS"]);
		    op_focei.covMethod=3;
		  } else {
		    // Now check sandwich matrix against R and S methods
		    bool covRSsmall = arma::any(abs(covRS.diag()) < op_focei.covSmall);
		    double covRSd= sum(covRS.diag());
		    arma::mat covR = as<arma::mat>(e["covR"]);
		    bool covRsmall = arma::any(abs(covR.diag()) < op_focei.covSmall);
		    double covRd= sum(covR.diag());
		    arma::mat covS = as<arma::mat>(e["covS"]);
		    bool covSsmall = arma::any(abs(covS.diag()) < op_focei.covSmall);
		    double  covSd= sum(covS.diag());
		    if ((covRSsmall && covSsmall && covRsmall)){
		      e["cov"] = covRS;
		    } else if (covRSsmall && !covSsmall && covRsmall) {
		      e["cov"] = covS;
		      op_focei.covMethod=3;
		    } else if (covRSsmall && covSsmall && !covRsmall) {
		      e["cov"] = covR;
		      op_focei.covMethod=2;
		    } else if (covRSd > covRd){
		      // SE(RS) > SE(R)
		      if (covRd > covSd){
			// SE(R) > SE(S)
			e["cov"] = covS;
			op_focei.covMethod=3;
		      } else {
			e["cov"] = covR;
			op_focei.covMethod=2;
		      }
		    } else if (covRSd > covSd){
		      e["cov"] = covS;
		      op_focei.covMethod=3;
		    } else {
		      e["cov"] = covRS;
		    }
		  }
		} else {
		  e["cov"] = covRS;
		}
	      } else {
		mat Sinv;
		bool success;
		success = inv(Sinv, trimatu(cholS));
		if (!success){
		  warning(_("S matrix seems singular; Using pseudo-inverse"));
		  Sinv = pinv(trimatu(cholS));
		}
		Sinv = Sinv * Sinv.t();
		op_focei.cur++;
		op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
		e["cov"]= 4 * Sinv;
	      }
	    }
	  } else {
	    if (op_focei.covMethod == 1){
	      RSprintf("\rS matrix calculation failed; Switch to R-matrix covariance.\n");
	      e["cov"] = wrap(e["covR"]);
	      op_focei.covMethod = 2;
	    } else {
	      op_focei.covMethod=0;
	      RSprintf("\rCould not calculate covariance matrix.\n");
	      warning(_("cannot calculate covariance"));
	      op_focei.cur++;
	      op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
	    }
	  }
	}
	op_focei.cur=op_focei.totTick;
	op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
	if (e.exists("cov")){
	  arma::mat cov = as<arma::mat>(e["cov"]);
	  if (all(cov.diag() < 1e-7)){
	    warning(_("The variance of all elements are unreasonably small, <1e-7"));
	    op_focei.covMethod=0;
	    e.remove("cov");
	  } else {
	    arma::mat Dcov(cov.n_rows,cov.n_rows,fill::zeros);
	    Dcov.diag() = (sqrt(cov.diag()));
	    arma::vec sd2=Dcov.diag();
	    Dcov = inv_sympd(Dcov);
	    arma::mat cor2 = Dcov * cov * Dcov;
	    cor2.diag()= sd2;
	    e["cor"] = cor2;
	  }
	}
	if (op_focei.covMethod==0){
	  warning(_("covariance step failed"));
	  e["covMethod"] = CharacterVector::create("failed");
	  NumericMatrix ret;
	  return ret;
	} else {
	  if (op_focei.covMethod == 1){
	    bool doWarn=false;
	    if (rstr == "|r|"){
	      warning(_("R matrix non-positive definite but corrected by R = sqrtm(R%%*%%R)"));
	      doWarn=true;
	    } else if (rstr == "r+"){
	      warning(_("R matrix non-positive definite but corrected (because of cholAccept)"));
	      doWarn=true;
	    }
	    if (sstr == "|s|"){
	      warning(_("S matrix non-positive definite but corrected by S = sqrtm(S%%*%%S)"));
	      doWarn=true;
	    } else if (sstr == "s+"){
	      warning(_("S matrix non-positive definite but corrected (because of cholAccept)"));
	      doWarn=true;
	    }
	    if (doWarn){
	      warning(_("since sandwich matrix is corrected, you may compare to $covR or $covS if you wish"));
	    }
	    rstr =  rstr + "," + sstr;
	    e["covMethod"] = wrap(rstr);
	  } else if (op_focei.covMethod == 2){
	    if (rstr == "|r|"){
	      warning(_("R matrix non-positive definite but corrected by R = sqrtm(R%%*%%R)"));
	    } else if (rstr == "r+"){
	      warning(_("R matrix non-positive definite but corrected (because of cholAccept)"));
	    }
	    e["covMethod"] = wrap(rstr);
	    if (origCov != 2){
	      if (checkSandwich){
		warning(_("using R matrix to calculate covariance, can check sandwich or S matrix with $covRS and $covS"));
	      } else {
		warning(_("using R matrix to calculate covariance"));
	      }
	    }
	  } else if (op_focei.covMethod == 3){
	    e["covMethod"] = wrap(sstr);
	    if (origCov != 2){
	      if (checkSandwich){
		warning(_("using S matrix to calculate covariance, can check sandwich or R matrix with $covRS and $covR"));
	      } else {
		warning(_("using S matrix to calculate covariance"));
	      }
	    }
	  }
	  return as<NumericMatrix>(e["cov"]);
	}
      } else {
	if (boundary){
	  warning(_("parameter estimate near boundary; covariance not calculated:\n   ") + boundStr + _("\n use 'getVarCov' to calculate anyway"));
	  e["covMethod"] = "Boundary issue; Get SEs with `getVarCov()`: " + boundStr;
	}
	op_focei.cur=op_focei.totTick;
	op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
	NumericMatrix ret;
	return ret;
      }
    }
  } catch (...) {
    warning(_("covariance step failed"));
    e["covMethod"] = CharacterVector::create("failed");
    NumericMatrix ret;
    return ret;
  }
  NumericMatrix ret;
  return ret;
}

void parHistData(Environment e, bool focei){
  if (!e.exists("method") && iterType.size() > 0) {
    CharacterVector thetaNames=as<CharacterVector>(e["thetaNames"]);
    CharacterVector dfNames;
    if (focei){
      CharacterVector dfNames2(3+op_focei.npars);
      dfNames = dfNames2;
    } else {
      CharacterVector dfNames2(3+thetaNames.size());
      dfNames = dfNames2;
    }
    dfNames[0] = "iter";
    dfNames[1] = "type";
    dfNames[2] = "objf";
    int i, j, k=1;
    if (focei){
      for (i = 0; i < op_focei.npars; i++){
	j=op_focei.fixedTrans[i];
	if (j < thetaNames.size()){
	  dfNames[i+3] = thetaNames[j];
	} else {
	  dfNames[i+3] = "o" + std::to_string(k++);
	}
      }
    } else {
      for (i = 0; i < thetaNames.size(); i++){
	dfNames[i+3] = thetaNames[i];
      }
    }
    // iter type parameters
    List ret;
    if (focei){
      ret = List(3+op_focei.npars);
    } else {
      ret = List(3+thetaNames.size());
    }
    int sz = niter.size()+niterGrad.size();
    IntegerVector tmp;
    std::vector<int> iter;
    iter.reserve(sz);
    iter.insert(iter.end(), niter.begin(), niter.end());
    iter.insert(iter.end(), niterGrad.begin(), niterGrad.end());
    ret[0] = iter;
    tmp = IntegerVector(sz);
    std::vector<int> typ;
    typ.reserve(sz);
    typ.insert(typ.end(), iterType.begin(), iterType.end());
    typ.insert(typ.end(), gradType.begin(), gradType.end());
    tmp = typ;
    tmp.attr("levels") = CharacterVector::create("Gill83 Gradient", "Mixed Gradient",
						 "Forward Difference", "Central Difference",
						 "Scaled", "Unscaled", "Back-Transformed");
    tmp.attr("class") = "factor";
    ret[1] = tmp;
    arma::mat cPar(vPar.size()/iterType.size(), iterType.size());
    std::copy(vPar.begin(), vPar.end(), cPar.begin());
    arma::mat vals;
    if (vGrad.size() > 0){
      arma::mat cGrad(vGrad.size()/gradType.size(), gradType.size());
      std::copy(vGrad.begin(), vGrad.end(), cGrad.begin());
      cPar = cPar.t();
      cGrad = cGrad.t();
      vals = arma::join_cols(cPar, cGrad);
    } else {
      cPar = cPar.t();
      vals = cPar;
    }
    if (focei){
      for (i = 0; i < min2(op_focei.npars+1, vals.n_cols); i++){
	ret[i+2]= vals.col(i);
      }
    } else {
      for (i = 0; i < thetaNames.size()+1; i++){
	ret[i+2]= vals.col(i);
      }
    }
    vGrad.clear();
    vPar.clear();
    iterType.clear();
    gradType.clear();
    niter.clear();
    niterGrad.clear();
    ret.attr("names")=dfNames;
    ret.attr("class") = "data.frame";
    ret.attr("row.names")=IntegerVector::create(NA_INTEGER, -sz);
    Function loadNamespace("loadNamespace", R_BaseNamespace);
    Environment nlmixr = loadNamespace("nlmixr");
    Environment thetaReset = nlmixr[".thetaReset"];
    if (thetaReset.exists("parHistData")) {
      // rbind  data.
      if (TYPEOF(thetaReset["parHistData"]) == VECSXP) {
	Function loadNamespace("loadNamespace", R_BaseNamespace);
	Environment nlmixr = loadNamespace("nlmixr");
	Function rbind = nlmixr[".rbindParHistory"];
	ret = rbind(thetaReset["parHistData"], ret);
      }
      thetaReset.remove("parHistData");
    }
    e["parHistData"] = ret;
  }
}

void foceiFinalizeTables(Environment e){
  CharacterVector thetaNames=as<CharacterVector>(e["thetaNames"]);
  arma::mat cov;
  bool covExists = e.exists("cov");
  if (covExists){
    if (RxODE::rxIs(e["cov"], "matrix")){
      cov= as<arma::mat>(e["cov"]);
    } else {
      covExists = false;
    }
  }
  LogicalVector skipCov = e["skipCov"];

  if (covExists && op_focei.eigen){
    arma::vec eigval;
    arma::mat eigvec;

    eig_sym(eigval, eigvec, cov);
    e["eigen"] = eigval;
    e["eigenVec"] = eigvec;
    unsigned int k=0;
    if (eigval.size() > 0){
      double mx=std::fabs(eigval[0]), mn, cur;
      mn=mx;
      for (k = eigval.size(); k--;){
        cur = std::fabs(eigval[k]);
        if (cur > mx){
          mx=cur;
        }
        if (cur < mn){
          mn=cur;
        }
      }
      e["conditionNumber"] = mx/mn;
    } else {
      e["conditionNumber"] = NA_REAL;
    }
  }
  arma::vec se1;
  if (covExists){
    se1 = sqrt(cov.diag());
  }
  DataFrame thetaDf = as<DataFrame>(e["theta"]);
  arma::vec theta = as<arma::vec>(thetaDf["theta"]);
  NumericVector se(theta.size());
  NumericVector cv(theta.size());
  std::fill_n(&se[0], theta.size(), NA_REAL);
  std::fill_n(&cv[0], theta.size(), NA_REAL);
  int j=0, k=0, l=0;
  if (covExists){
    for (int k = 0; k < se.size(); k++){
      if (k >= skipCov.size()) break;
      if (!skipCov[k]){
        se[k] = se1[j++];
        cv[k] = std::fabs(se[k]/theta[k])*100;
      }
    }
  }
  e["se"] = se;
  List popDf = List::create(_["Estimate"]=thetaDf["theta"], _["SE"]=se,
			      _["%RSE"]=cv);
  popDf.attr("class") = "data.frame";
  popDf.attr("row.names") = IntegerVector::create(NA_INTEGER,-theta.size());
  e["popDf"] = popDf;

  e["fixef"]=thetaDf["theta"];
  NumericMatrix tmpNM;
  List tmpL;
  int i;
  if (!Rf_isNull(e["etaObf"])) {
    List etas = e["etaObf"];
    IntegerVector idx = seq_len(etas.length())-1;
    etas = etas[idx != etas.length()-1];
    e["ranef"]=etas;

    // Now put names on the objects
    ////////////////////////////////////////////////////////////////////////////////
    // Eta Names
    CharacterVector etaNames=as<CharacterVector>(e["etaNames"]);
    tmpNM = getOmega();
    tmpNM.attr("dimnames") = List::create(etaNames, etaNames);
    e["omega"] = tmpNM;

    tmpNM = as<NumericMatrix>(e["omegaR"]);
    tmpNM.attr("dimnames") = List::create(etaNames, etaNames);
    e["omegaR"] = tmpNM;

    tmpL  = as<List>(e["ranef"]);
    List tmpL2 = as<List>(e["etaObf"]);
    CharacterVector tmpN  = tmpL.attr("names");
    CharacterVector tmpN2 = tmpL2.attr("names");
    for (i = 0; i < etaNames.size(); i++){
      if (i + 1 <  tmpN.size())  tmpN[i+1] = etaNames[i];
      if (i + 1 < tmpN2.size()) tmpN2[i+1] = etaNames[i];
    }
    ////////////////////////////////////////////////////////////////////////////////
    tmpL.attr("names") = tmpN;
    tmpL2.attr("names") = tmpN2;
    e["ranef"] = tmpL;
    e["etaObf"] = tmpL2;

    // omegaR
    arma::mat omega = as<arma::mat>(e["omega"]);
    arma::mat D(omega.n_rows,omega.n_rows,fill::zeros);
    arma::mat cor(omega.n_rows,omega.n_rows);
    D.diag() = (sqrt(omega.diag()));
    arma::vec sd=D.diag();
    D = inv_sympd(D);
    cor = D * omega * D;
    cor.diag()= sd;
  } else {
    e["ranef"] = R_NilValue;
    e["omegaR"] = R_NilValue;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Theta names
  //
  tmpL = as<List>(e["theta"]);
  tmpL.attr("row.names") = thetaNames;
  e["theta"] = tmpL;

  tmpL=e["popDf"];
  // Add a few columns

  NumericVector Estimate = tmpL["Estimate"];
  NumericVector SE = tmpL["SE"];
  NumericVector RSE = tmpL["%RSE"];
  NumericVector EstBT(Estimate.size());
  NumericVector EstLower(Estimate.size());
  NumericVector EstUpper(Estimate.size());

  CharacterVector EstS(Estimate.size());
  CharacterVector SeS(Estimate.size());
  CharacterVector rseS(Estimate.size());
  CharacterVector btCi(Estimate.size());
  // LogicalVector EstBT(Estimate.size());
  // Rf_pt(stat[7],(double)n1,1,0)
  // FIXME figure out log thetas outside of foceisetup.
  IntegerVector logTheta;
  IntegerVector logitTheta;
  NumericVector logitThetaHi;
  NumericVector logitThetaLow;
  IntegerVector probitTheta;
  NumericVector probitThetaHi;
  NumericVector probitThetaLow;
  if (e.exists("logThetasF")){
    logTheta =  as<IntegerVector>(e["logThetasF"]);
  } else if (e.exists("model")){
    List model = e["model"];
    logTheta =  as<IntegerVector>(model["log.thetas"]);
  }
  if (e.exists("logitThetasF")) {
    logitTheta = as<IntegerVector>(e["logitThetasF"]);
    if (logitTheta.size() != 0) {
      logitThetaHi = as<NumericVector>(e["logitThetasHiF"]);
      logitThetaLow = as<NumericVector>(e["logitThetasLowF"]);
    }
  }
  if (e.exists("probitThetasF")) {
    probitTheta = as<IntegerVector>(e["probitThetasF"]);
    if (probitTheta.size() != 0) {
      probitThetaHi = as<NumericVector>(e["probitThetasHiF"]);
      probitThetaLow = as<NumericVector>(e["probitThetasLowF"]);
    }
  }
  j = logTheta.size()-1;
  k = logitTheta.size()-1;
  l = probitTheta.size()-1;
  double qn= Rf_qnorm5(1.0-(1-op_focei.ci)/2, 0.0, 1.0, 1, 0);
  std::string cur;
  char buff[100];
  LogicalVector thetaFixed =thetaDf["fixed"];
  for (i = Estimate.size(); i--;){
    snprintf(buff, sizeof(buff), "%.*g", (int)op_focei.sigdig, Estimate[i]);
    EstS[i]=buff;
    if (logTheta.size() > 0 && j >= 0 && logTheta[j]-1==i) {
      EstBT[i] = exp(Estimate[i]);
      snprintf(buff, sizeof(buff), "%.*g", (int)op_focei.sigdig, EstBT[i]);
      cur = buff;
      if (ISNA(SE[i])){
        EstLower[i] = NA_REAL;
        EstUpper[i] = NA_REAL;
	if (thetaFixed[i]){
          SeS[i]  = "FIXED";
          rseS[i] = "FIXED";
	} else {
          SeS[i] = "";
          rseS[i]="";
        }
      } else {
        EstLower[i] = exp(Estimate[i]-SE[i]*qn);
        EstUpper[i] = exp(Estimate[i]+SE[i]*qn);
	snprintf(buff, sizeof(buff), "%.*g", (int)op_focei.sigdig, SE[i]);
        SeS[i]=buff;
	snprintf(buff, sizeof(buff), "%.*g", (int)op_focei.sigdig, RSE[i]);
        rseS[i]=buff;
        snprintf(buff, sizeof(buff), "%.*g", (int)op_focei.sigdig, EstLower[i]);
	cur = cur + " (" + buff + ", ";
        snprintf(buff, sizeof(buff), "%.*g", (int)op_focei.sigdig, EstUpper[i]);
	cur = cur + buff + ")";
      }
      btCi[i] = cur;
      j--;
    } else if (logitTheta.size() > 0 && k >= 0 && logitTheta[k]-1==i) {
      EstBT[i] = expit(Estimate[i], logitThetaLow[k], logitThetaHi[k]);
      snprintf(buff, sizeof(buff), "%.*g", (int)op_focei.sigdig, EstBT[i]);
      cur = buff;
      if (ISNA(SE[i])){
        EstLower[i] = NA_REAL;
        EstUpper[i] = NA_REAL;
	if (thetaFixed[i]){
          SeS[i]  = "FIXED";
          rseS[i] = "FIXED";
	} else {
          SeS[i] = "";
          rseS[i]="";
        }
      } else {
        EstLower[i] = expit(Estimate[i]-SE[i]*qn, logitThetaLow[k], logitThetaHi[k]);
        EstUpper[i] = expit(Estimate[i]+SE[i]*qn, logitThetaLow[k], logitThetaHi[k]);
	snprintf(buff, sizeof(buff), "%.*g", (int)op_focei.sigdig, SE[i]);
        SeS[i]=buff;
	snprintf(buff, sizeof(buff), "%.*g", (int)op_focei.sigdig, RSE[i]);
        rseS[i]=buff;
        snprintf(buff, sizeof(buff), "%.*g", (int)op_focei.sigdig, EstLower[i]);
	cur = cur + " (" + buff + ", ";
        snprintf(buff, sizeof(buff), "%.*g", (int)op_focei.sigdig, EstUpper[i]);
	cur = cur + buff + ")";
      }
      btCi[i] = cur;
      k--;
    } else if (probitTheta.size() > 0 && l >= 0 && probitTheta[l]-1==i) {
      EstBT[i] = probitInv(Estimate[i], probitThetaLow[l], probitThetaHi[l]);
      snprintf(buff, sizeof(buff), "%.*g", (int)op_focei.sigdig, EstBT[i]);
      cur = buff;
      if (ISNA(SE[i])){
        EstLower[i] = NA_REAL;
        EstUpper[i] = NA_REAL;
	if (thetaFixed[i]){
          SeS[i]  = "FIXED";
          rseS[i] = "FIXED";
	} else {
          SeS[i] = "";
          rseS[i]="";
        }
      } else {
        EstLower[i] = probitInv(Estimate[i]-SE[i]*qn, probitThetaLow[l], probitThetaHi[l]);
        EstUpper[i] = probitInv(Estimate[i]+SE[i]*qn, probitThetaLow[l], probitThetaHi[l]);
	snprintf(buff, sizeof(buff), "%.*g", (int)op_focei.sigdig, SE[i]);
        SeS[i]=buff;
	snprintf(buff, sizeof(buff), "%.*g", (int)op_focei.sigdig, RSE[i]);
        rseS[i]=buff;
        snprintf(buff, sizeof(buff), "%.*g", (int)op_focei.sigdig, EstLower[i]);
	cur = cur + " (" + buff + ", ";
        snprintf(buff, sizeof(buff), "%.*g", (int)op_focei.sigdig, EstUpper[i]);
	cur = cur + buff + ")";
      }
      btCi[i] = cur;
      l--;
    } else {
      EstBT[i]= Estimate[i];
      snprintf(buff, sizeof(buff), "%.*g", (int)op_focei.sigdig, Estimate[i]);
      EstS[i]=buff;
      snprintf(buff, sizeof(buff), "%.*g", (int)op_focei.sigdig, EstBT[i]);
      cur = buff;
      if (ISNA(SE[i])){
        EstLower[i] = NA_REAL;
        EstUpper[i] = NA_REAL;
        if (thetaFixed[i]){
          SeS[i]  = "FIXED";
          rseS[i] = "FIXED";
        } else {
          SeS[i] = "";
          rseS[i]="";
        }
      } else {
        EstLower[i] = Estimate[i]-SE[i]*qn;
        EstUpper[i] = Estimate[i]+SE[i]*qn;
	snprintf(buff, sizeof(buff), "%.*g", (int)op_focei.sigdig, SE[i]);
        SeS[i]=buff;
	snprintf(buff, sizeof(buff), "%.*g", (int)op_focei.sigdig, RSE[i]);
        rseS[i]=buff;
	snprintf(buff, sizeof(buff), "%.*g", (int)op_focei.sigdig, EstLower[i]);
	cur = cur + " (" + buff + ", ";
        snprintf(buff, sizeof(buff), "%.*g", (int)op_focei.sigdig, EstUpper[i]);
	cur = cur + buff + ")";
      }
      btCi[i] = cur;
    }
  }
  tmpL["Back-transformed"] = EstBT;
  tmpL["CI Lower"] = EstLower;
  tmpL["CI Upper"] = EstUpper;
  tmpL.attr("row.names") = thetaNames;
  tmpL.attr("class") = "data.frame";
  e["popDf"]=tmpL;
  std::string bt = "Back-transformed(" + std::to_string((int)(op_focei.ci*100)) + "%CI)";

  List popDfSig;
  if (e.exists("cov") && RxODE::rxIs(e["cov"], "matrix")){
    popDfSig = List::create(_["Est."]=EstS,
                   _["SE"]=SeS,
                   _["%RSE"]=rseS,
                   _[bt]=btCi);
  } else {
    popDfSig = List::create(_["Est."]=EstS,
			    _["Back-transformed"] = btCi);
  }

  popDfSig.attr("row.names") = thetaNames;
  popDfSig.attr("class") = "data.frame";
  e["popDfSig"]=popDfSig;

  NumericVector tmpNV = e["fixef"];
  tmpNV.names() = thetaNames;
  e["fixef"] = tmpNV;


  tmpNV = e["se"];
  tmpNV.names() = thetaNames;
  e["se"] = tmpNV;

  // Now get covariance names
  if (e.exists("cov")  && RxODE::rxIs(e["cov"], "matrix")){
    tmpNM = as<NumericMatrix>(e["cov"]);
    CharacterVector thetaCovN(tmpNM.nrow());
    LogicalVector skipCov = e["skipCov"];
    int j=0;
    for (unsigned int k = 0; k < thetaNames.size(); k++){
      if (k >= skipCov.size()) break;
      if (j >= thetaCovN.size()) break;
      if (!skipCov[k]){
        thetaCovN[j++] = thetaNames[k];
      }
    }
    List thetaDim = List::create(thetaCovN,thetaCovN);
    tmpNM.attr("dimnames") = thetaDim;
    e["cov"]=tmpNM;
    if (e.exists("cor") && RxODE::rxIs(e["cor"], "matrix")){
      tmpNM = as<NumericMatrix>(e["cor"]);
      tmpNM.attr("dimnames") = thetaDim;
      e["cor"]=tmpNM;
    } else {
      arma::mat cov0 = as<arma::mat>(e["cov"]);
      arma::mat Dcov(cov0.n_rows,cov0.n_rows,fill::zeros);
      Dcov.diag() = (sqrt(cov0.diag()));
      arma::vec sd2=Dcov.diag();
      Dcov = inv_sympd(Dcov);
      arma::mat cor2 = Dcov * cov0 * Dcov;
      cor2.diag()= sd2;
      tmpNM = wrap(cor2);
      tmpNM.attr("dimnames") = thetaDim;
      e["cor"] = tmpNM;
    }
    if (e.exists("Rinv") && RxODE::rxIs(e["Rinv"], "matrix")){
      tmpNM = as<NumericMatrix>(e["Rinv"]);
      tmpNM.attr("dimnames") = thetaDim;
      e["Rinv"]=tmpNM;
    }
    if (e.exists("Sinv") && RxODE::rxIs(e["Sinv"], "matrix")){
      tmpNM = as<NumericMatrix>(e["Sinv"]);
      tmpNM.attr("dimnames") = thetaDim;
      e["Sinv"]=tmpNM;
    }
    if (e.exists("S") && RxODE::rxIs(e["S"], "matrix")){
      tmpNM = as<NumericMatrix>(e["S"]);
      tmpNM.attr("dimnames") = thetaDim;
      e["S"]=tmpNM;
    }
    if (e.exists("R") && RxODE::rxIs(e["R"], "matrix")){
      tmpNM = as<NumericMatrix>(e["R"]);
      tmpNM.attr("dimnames") = thetaDim;
      e["R"]=tmpNM;
    }
    if (e.exists("covR") && RxODE::rxIs(e["covR"], "matrix")){
      tmpNM = as<NumericMatrix>(e["covR"]);
      tmpNM.attr("dimnames") = thetaDim;
      e["covR"]=tmpNM;
    }
    if (e.exists("covRS") && RxODE::rxIs(e["covRS"], "matrix")){
      tmpNM = as<NumericMatrix>(e["covRS"]);
      tmpNM.attr("dimnames") = thetaDim;
      e["covRS"]=tmpNM;
    }
    if (e.exists("covS") && RxODE::rxIs(e["covS"], "matrix")){
      tmpNM = as<NumericMatrix>(e["covS"]);
      tmpNM.attr("dimnames") = thetaDim;
      e["covS"]=tmpNM;
    }
    if (e.exists("R.1") && RxODE::rxIs(e["R.1"], "matrix")){
      tmpNM = as<NumericMatrix>(e["R.1"]);
      tmpNM.attr("dimnames") = thetaDim;
      e["R.1"]=tmpNM;
    }
    if (e.exists("R.2") && RxODE::rxIs(e["R.2"], "matrix")){
      tmpNM = as<NumericMatrix>(e["R.2"]);
      tmpNM.attr("dimnames") = thetaDim;
      e["R.2"]=tmpNM;
    }
    if (e.exists("cholR") && RxODE::rxIs(e["cholR"], "matrix")){
      tmpNM = as<NumericMatrix>(e["cholR"]);
      tmpNM.attr("dimnames") = thetaDim;
      e["cholR"]=tmpNM;
    }
    if (e.exists("cholR2") && RxODE::rxIs(e["cholR2"], "matrix")){
      tmpNM = as<NumericMatrix>(e["cholR2"]);
      tmpNM.attr("dimnames") = thetaDim;
      e["cholR2"]=tmpNM;
    }
    if (e.exists("cholS") && RxODE::rxIs(e["cholS"], "matrix")){
      tmpNM = as<NumericMatrix>(e["cholS"]);
      tmpNM.attr("dimnames") = thetaDim;
      e["cholS"]=tmpNM;
    }
  }
  List objDf;
  if (e.exists("conditionNumber")){
    objDf = List::create(_["OBJF"] = as<double>(e["objective"]), _["AIC"]=as<double>(e["AIC"]),
                         _["BIC"] = as<double>(e["BIC"]), _["Log-likelihood"]=as<double>(e["logLik"]),
                         _["Condition Number"]=as<double>(e["conditionNumber"]));
  } else {
    objDf = List::create(_["OBJF"] = as<double>(e["objective"]), _["AIC"]=as<double>(e["AIC"]),
                         _["BIC"] = as<double>(e["BIC"]), _["Log-likelihood"]=as<double>(e["logLik"]));
  }
  if (op_focei.neta == 0){
    objDf.attr("row.names") = CharacterVector::create("Pop");
    e["ofvType"] = "Pop";
  } else if (op_focei.fo == 1){
    objDf.attr("row.names") = CharacterVector::create("FO");
    e["ofvType"] = "fo";
  } else if (op_focei.interaction){
    objDf.attr("row.names") = CharacterVector::create("FOCEi");
    e["ofvType"] = "focei";
  } else {
    objDf.attr("row.names") = CharacterVector::create("FOCE");
    e["ofvType"] = "foce";
  }
  objDf.attr("class") = "data.frame";
  e["objDf"]=objDf;
  if (!e.exists("method")){
    if (op_focei.neta == 0){
      e["method"] = "Population Only";
    } else if (op_focei.fo == 1){
      e["method"] = "FO";
    } else {
      e["method"] = "FOCE";
    }
  }
  if (!e.exists("extra")){
    if (op_focei.neta == 0){
      e["extra"] = "";
      e["skipTable"] = NA_LOGICAL;
    } else if (op_focei.fo == 1){
      e["extra"] = "";
      e["skipTable"] = LogicalVector::create(true);
    } else if (op_focei.interaction){
      if(op_focei.useColor){
        e["extra"] = "\033[31;1mi\033[0m";
      } else {
        e["extra"] = "i";
      }
    } else {
      e["extra"] = "";
    }
    List ctl = e["control"];
    e["extra"] = as<std::string>(e["extra"]) + " (outer: " + as<std::string>(ctl["outerOptTxt"]) +")";
  }
  // RxODE::rxSolveFree();
  e.attr("class") = "nlmixrFitCore";
}

////////////////////////////////////////////////////////////////////////////////
// FOCEi fit

//' Fit/Evaulate FOCEi
//'
//' This shouldn't be called directly.
//'
//' @param e Enviornment
//'
//' @return A focei fit object
//'
//' @keywords internal
//' @export
//[[Rcpp::export]]
Environment foceiFitCpp_(Environment e){
  if (!assignFn_){
    n1qn1_ = (n1qn1_fp) R_GetCCallable("n1qn1","n1qn1_");
    par_progress = (par_progress_t) R_GetCCallable("RxODE", "par_progress");
    getRx = (getRxSolve_t) R_GetCCallable("RxODE", "getRxSolve_");
    isRstudio = (isRstudio_t) R_GetCCallable("RxODE", "isRstudio");
    ind_solve=(ind_solve_t) R_GetCCallable("RxODE", "ind_solve");
    rxGetId = (rxGetId_t) R_GetCCallable("RxODE", "rxGetId");
    getTimeF = (getTime_t) R_GetCCallable("RxODE", "getTime");
    iniSubjectI = (iniSubjectI_t) R_GetCCallable("RxODE","iniSubjectE");
    sortIdsF = (sortIds_t) R_GetCCallable("RxODE", "sortIds");
    assignFn_=true;
  }
  clock_t t0 = clock();
  List model = e["model"];
  bool doPredOnly = false;
  op_focei.canDoFD = false;
  if (model.containsElementNamed("inner")) {
    RObject inner = model["inner"];
    if (RxODE::rxIs(inner, "RxODE")) {
      foceiSetup_(inner, as<RObject>(e["dataSav"]),
		  as<NumericVector>(e["thetaIni"]), e["thetaFixed"], e["skipCov"],
		  as<RObject>(e["rxInv"]), e["lower"], e["upper"], e["etaMat"],
		  e["control"]);
      if (model.containsElementNamed("pred.nolhs")){
	RObject noLhs = model["pred.nolhs"];
	if (RxODE::rxIs(noLhs, "RxODE")) {
	  List mvp = RxODE::rxModelVars_(noLhs);
	  rxUpdateFuns(as<SEXP>(mvp["trans"]), &rxPred);
	  op_focei.canDoFD = true;
	}
      }
      // Now setup which ETAs need a finite difference
      if (model.containsElementNamed("eventEta")) {
	IntegerVector eventEta = model["eventEta"];
	std::copy(eventEta.begin(), eventEta.end(),&op_focei.etaFD[0]);
      }
    } else if (model.containsElementNamed("pred.only")){
      inner = model["pred.only"];
      if (RxODE::rxIs(inner, "RxODE")){
	foceiSetup_(inner, as<RObject>(e["dataSav"]),
		    as<NumericVector>(e["thetaIni"]), e["thetaFixed"], e["skipCov"],
		    as<RObject>(e["rxInv"]), e["lower"], e["upper"], e["etaMat"],
		    e["control"]);
	doPredOnly = true;
	if (op_focei.neta == 0) doPredOnly = false;
      } else {
	stop(_("cannot run this function"));
      }
    } else {
      doPredOnly=true;
      foceiSetupTrans_(as<CharacterVector>(e[".params"]));
      foceiSetup_(R_NilValue, as<RObject>(e["dataSav"]),
                  as<NumericVector>(e["thetaIni"]), e["thetaFixed"], e["skipCov"],
                  as<RObject>(e["rxInv"]), e["lower"], e["upper"], e["etaMat"],
                  e["control"]);
    }
  } else {
    doPredOnly=true;
    foceiSetupTrans_(as<CharacterVector>(e[".params"]));
    foceiSetup_(R_NilValue, as<RObject>(e["dataSav"]),
                as<NumericVector>(e["thetaIni"]), e["thetaFixed"], e["skipCov"],
                as<RObject>(e["rxInv"]), e["lower"], e["upper"], e["etaMat"],
                e["control"]);
  }
  if (e.exists("setupTime")){
    e["setupTime"] = as<double>(e["setupTime"])+(((double)(clock() - t0))/CLOCKS_PER_SEC);
  } else {
    e["setupTime"] = (((double)(clock() - t0))/CLOCKS_PER_SEC);
  }
  t0 = clock();
  CharacterVector thetaNames=as<CharacterVector>(e["thetaNames"]);
  IntegerVector logTheta;
  IntegerVector logitTheta;
  IntegerVector xType = e["xType"];
  std::fill_n(&op_focei.scaleC[0], op_focei.ntheta+op_focei.omegan, NA_REAL);
  if (e.exists("scaleC")){
    arma::vec scaleC = as<arma::vec>(e["scaleC"]);
    std::copy(scaleC.begin(), scaleC.end(), &op_focei.scaleC[0]);
  }
  if (e.exists("logThetas")){
    logTheta =  as<IntegerVector>(e["logThetas"]);
  } else if (e.exists("model")){
    List model = e["model"];
    if (model.containsElementNamed("log.thetas")){
      if (RxODE::rxIs(model["log.thetas"], "integer") || RxODE::rxIs(model["log.thetas"], "numeric")){
	logTheta =  as<IntegerVector>(model["log.thetas"]);
      }
    }
  }
  if (e.exists("logitThetas") && e.exists("logitThetasLow") &&
      e.exists("logitThetasHi")) {
    logitTheta = as<IntegerVector>(e["logitThetas"]);
    if (logitTheta.size() != 0) {
      op_focei.logitThetaLow = as<NumericVector>(e["logitThetasLow"]);
      op_focei.logitThetaHi = as<NumericVector>(e["logitThetasHi"]);
    } else {
      op_focei.logitThetaLow=NumericVector(0);
      op_focei.logitThetaHi=NumericVector(0);
    }
  } else {
    op_focei.logitThetaLow=NumericVector(0);
    op_focei.logitThetaHi=NumericVector(0);
  }
  int j;
  // Setup which parameters are transformed
  for (unsigned int k = op_focei.npars; k--;){
    j=op_focei.fixedTrans[k];
    op_focei.xPar[k] = 0;
    if ((int)op_focei.ntheta < j){
      op_focei.xPar[k] = xType[j-op_focei.ntheta];
    } else {
      for (unsigned int m=logTheta.size(); m--;){
	if (logTheta[m]-1 == j){
	  op_focei.xPar[k] = 1;
	  break;
	}
      }
      for (int m = logitTheta.size(); m--;) {
	if (logitTheta[m]-1 == j){
	  op_focei.xPar[k] = -m-1;
	  break;
	}
      }
    }
  }
  std::string tmpS;
  if (op_focei.maxOuterIterations > 0 && op_focei.printTop == 1 && op_focei.printOuter != 0){
    if (op_focei.useColor)
      RSprintf("\033[1mKey:\033[0m ");
    else
      RSprintf("Key: ");

    RSprintf("U: Unscaled Parameters; ");
    RSprintf("X: Back-transformed parameters; ");
    RSprintf("G: Gill difference gradient approximation\n");
    RSprintf("F: Forward difference gradient approximation\n");
    RSprintf("C: Central difference gradient approximation\n");
    RSprintf("M: Mixed forward and central difference gradient approximation\n");
    RSprintf("Unscaled parameters for Omegas=chol(solve(omega));\nDiagonals are transformed, as specified by foceiControl(diagXform=)\n");
    op_focei.t0 = clock();
    foceiPrintLine(min2(op_focei.npars, op_focei.printNcol));
    RSprintf("|    #| Objective Fun |");
    int j,  i=0, finalize=0, k=1;

    for (i = 0; i < op_focei.npars; i++){
      j=op_focei.fixedTrans[i];
      if (j < thetaNames.size()){
	tmpS = thetaNames[j];
	RSprintf("%#10s |", tmpS.c_str());
      } else {
	tmpS = "o" +std::to_string(k++);
	RSprintf("%#10s |", tmpS.c_str());
      }
      if ((i + 1) != op_focei.npars && (i + 1) % op_focei.printNcol == 0){
	if (op_focei.useColor && op_focei.printNcol + i  >= op_focei.npars){
	  RSprintf("\n\033[4m|.....................|");
	} else {
	  RSprintf("\n|.....................|");
	}
	finalize=1;
      }
    }
    if (finalize){
      while(true){
	if ((i++) % op_focei.printNcol == 0){
	  if (op_focei.useColor) RSprintf("\033[0m");
	  RSprintf("\n");
	  break;
	} else {
	  RSprintf("...........|");
	}
      }
    } else {
      RSprintf("\n");
    }
    if (!op_focei.useColor){
      foceiPrintLine(min2(op_focei.npars, op_focei.printNcol));
    }
  }
  if (doPredOnly){
    if (e.exists("objective")){
      nlmixrEnvSetup(e, as<double>(e["objective"]));
    } else {
      stop(_("not setup right"));
    }
  } else {
    op_focei.didHessianReset=0;
    op_focei.didEtaNudge =0;
    op_focei.didEtaReset=0;
    op_focei.stickyRecalcN2=0;
    op_focei.stickyRecalcN1=0;
    foceiOuter(e);
    if (op_focei.didHessianReset==1){
      warning(_("Hessian reset during optimization; (can control by foceiControl(resetHessianAndEta=.))"));
    }
    if (op_focei.didEtaNudge==1){
      warning(_("initial ETAs were nudged; (can control by foceiControl(etaNudge=., etaNudge2=))"));
    }
    if (op_focei.didEtaReset==1){
      warning(_("ETAs were reset to zero during optimization; (Can control by foceiControl(resetEtaP=.))"));
    }
    if (op_focei.repeatGillN > 0){
      warning(_("tolerances were reduced during Gill Gradient, so it was repeated %d/%d times\nYou can control this with foceiControl(repeatGillMax=.)"), op_focei.repeatGillN, op_focei.repeatGillMax);
    }
    if (op_focei.neta != 0 && op_focei.maxOuterIterations > 0  && R_FINITE(op_focei.resetThetaFinalSize)) {
      focei_options *fop = &op_focei;
	std::fill(op_focei.etaM.begin(),op_focei.etaM.end(), 0.0);
	std::fill(op_focei.etaS.begin(),op_focei.etaS.end(), 0.0);
      double n = 1.0;
      for (int id=rx->nsub; id--;){
	focei_ind *fInd = &(inds_focei[id]);
	mat etaMat(fop->neta, 1);
	std::copy(&fInd->eta[0], &fInd->eta[0] + op_focei.neta, etaMat.begin());
	mat oldM = op_focei.etaM;
	op_focei.etaM = op_focei.etaM + (etaMat - op_focei.etaM)/n;
	op_focei.etaS = op_focei.etaS + (etaMat - op_focei.etaM) %  (etaMat - oldM);
	n += 1.0;
      }
      op_focei.eta1SD = 1/sqrt(op_focei.etaS);
      thetaReset(op_focei.resetThetaFinalSize);
    }
  }
  NumericVector scaleSave(op_focei.ntheta+op_focei.omegan);
  for (unsigned int i =op_focei.ntheta+op_focei.omegan;i--;){
    scaleSave[i] = getScaleC(i);
  }
  e["scaleC"] = scaleSave;
  parHistData(e, true); // Need to calculate before the parameter translations are mangled
  thetaResetObj(e);
  IntegerVector gillRet(op_focei.ntheta+op_focei.omegan);
  NumericVector gillAEps(op_focei.ntheta+op_focei.omegan,NA_REAL);
  NumericVector gillREps(op_focei.ntheta+op_focei.omegan,NA_REAL);
  NumericVector gillAEpsC(op_focei.ntheta+op_focei.omegan,NA_REAL);
  NumericVector gillREpsC(op_focei.ntheta+op_focei.omegan,NA_REAL);
  NumericVector gillCAEpsC(op_focei.ntheta+op_focei.omegan,NA_REAL);
  NumericVector gillCREpsC(op_focei.ntheta+op_focei.omegan,NA_REAL);
  bool warnGill = false;
  j = op_focei.npars;
  for (int i = op_focei.ntheta+op_focei.omegan; i--;){
    gillRet[i] = op_focei.gillRet[i]+1;
    if (gillRet[i] != 1) {
      gillAEps[i] = op_focei.aEps[--j];
      gillREps[i] = op_focei.rEps[j];
      gillAEpsC[i] = op_focei.aEpsC[j];
      gillREpsC[i] = op_focei.rEpsC[j];
    }
    if (gillRet[i] >= 3) warnGill=true;
  }
  CharacterVector gillLvl = CharacterVector::create("Not Assessed","Good","High Grad Error", "Constant Grad","Odd/Linear Grad",
						    "Grad changes quickly");
  gillRet.attr("levels") = gillLvl;
  gillRet.attr("class") = "factor";
  e["gillRet"] = gillRet;
  t0 = clock();
  foceiCalcCov(e);
  if (op_focei.didPredSolve) {
    warning(_("numerical difficulties solving forward sensitivity inner problem, tried approximating with more inaccurate numeric differences"));
  }
  IntegerVector gillRetC(op_focei.ntheta+op_focei.omegan);
  bool warnGillC = false;
  j = op_focei.npars;
  for (int i = op_focei.ntheta+op_focei.omegan; i--;){
    gillRetC[i] = op_focei.gillRetC[i]+1;
    if (gillRetC[i] >= 3) warnGillC=true;
    if (gillRetC[i] != 1) {
      gillCAEpsC[i] = op_focei.aEpsC[--j];
      gillCREpsC[i] = op_focei.rEpsC[j];
    }
  }
  gillRetC.attr("levels") = gillLvl;
  gillRetC.attr("class") = "factor";
  e["gillRetC"] = gillRetC;
  e["optimTime"] = (((double)(clock() - t0))/CLOCKS_PER_SEC);

  e["covTime"] = (((double)(clock() - t0))/CLOCKS_PER_SEC);
  List timeDf = List::create(_["setup"]=as<double>(e["setupTime"]),
			     _["optimize"]=as<double>(e["optimTime"]),
			     _["covariance"]=as<double>(e["covTime"]));
  timeDf.attr("class") = "data.frame";
  timeDf.attr("row.names") = "";
  e["time"] = timeDf;
  List scaleInfo = List::create(as<NumericVector>(e["fullTheta"]),
				as<NumericVector>(e["scaleC"]), gillRet,
				gillAEps,
				gillREps,
				gillAEpsC,
				gillREpsC,
				gillRetC,
				gillCAEpsC,
				gillCREpsC);
  scaleInfo.attr("names") = CharacterVector::create("est","scaleC","Initial Gradient",
						    "Forward aEps","Forward rEps",
						    "Central aEps","Central rEps",
						    "Covariance Gradient",
						    "Covariance aEps","Covariance rEps");
  scaleInfo.attr("class") = "data.frame";
  scaleInfo.attr("row.names") = IntegerVector::create(NA_INTEGER,-gillRet.size());
  e["scaleInfo"] = scaleInfo;
  if (warnGillC && warnGill){
    warning(_("gradient problems with initial estimate and covariance; see $scaleInfo"));
  } else if (warnGill){
    warning(_("gradient problems with initial estimate; see $scaleInfo"));
  } else if (warnGillC){
    warning(_("gradient problems with covariance; see $scaleInfo"));
  }
  if (op_focei.reducedTol){
    if (op_focei.stickyTol){
      warning(_("tolerances (atol/rtol) were increased (after %d bad solves) for some difficult ODE solving during the optimization.\ncan control with foceiControl(stickyRecalcN=)\nconsider increasing sigdig/atol/rtol changing initial estimates or changing the structural model"), op_focei.stickyRecalcN);
    } else {
      warning(_("tolerances (atol/rtol) were temporarily increased for some difficult ODE solving during the optimization.\nconsider increasing sigdig/atol/rtol changing initial estimates or changing the structural model"));
    }
  }
  if (op_focei.zeroGrad){
    warning(_("zero gradient replaced with small number (%f)"), sqrt(DBL_EPSILON));
  }
  foceiFinalizeTables(e);
  // NumericVector scaleC(op_focei.ntheta+op_focei.omegan);
  // std::copy(&op_focei.scaleC[0], &op_focei.scaleC[0]+op_focei.ntheta+op_focei.omegan, scaleC.begin());
  // e["scaleC"]= scaleC;
  if (op_focei.maxOuterIterations){
    RSprintf(_("done\n"));
  }
  return e;
}

//[[Rcpp::export]]
NumericVector boxCox_(NumericVector x = 1, double lambda=1, int yj = 0){
  NumericVector ret(x.size());
  for (unsigned int i = x.size(); i--;){
    ret[i] = _powerD(x[i], lambda, yj, 0.0, 1.0);
  }
  return ret;
}

//[[Rcpp::export]]
NumericVector iBoxCox_(NumericVector x = 1, double lambda=1, int yj = 0){
  NumericVector ret(x.size());
  for (unsigned int i = x.size(); i--;){
    ret[i] = _powerDi(x[i], lambda, yj, 0.0, 1.0);
  }
  return ret;
}
