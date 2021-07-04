#define STRICT_R_HEADER
#include "armahead.h"


using namespace Rcpp;
using namespace arma;

arma::mat gershNested(arma::mat A, int j, int n) {
  arma::mat g(n, 1, fill::zeros);
  double sumToI, sumAfterI;
  for (int ii = j; ii < n; ++ii){
    if (ii == 0){
      sumToI=0.0;
    } else if (j == ii){
      sumToI=arma::sum(arma::abs(A(ii, span(ii-1, j))));
    } else {
      sumToI=arma::sum(arma::abs(A(ii, span(j, ii-1))));
    }
    if (ii == n-1){
      sumAfterI = 0;
    } else {
      sumAfterI = arma::sum(arma::abs(A(span(ii+1, n-1), ii)));
    }
    g(ii, 0) = sumToI+sumAfterI-A(ii,ii);
  }
  return g;
}
// Suggested from https://gking.harvard.edu/files/help.pdf
// Translated from
//http://www.dynare.org/dynare-matlab-m2html/matlab/chol_SE.html
// Use tau1=sqrt(eps) instead of eps^1/3; In my tests eps^1/3 produces NaNs
bool cholSE0(arma::mat &Ao, arma::mat &E, arma::mat A, double tol) {
  int n = A.n_rows;
  double tau1 = tol;//
  double tau2 = tol;//tau1;
  bool phase1 = true;
  double delta = 0;
  int j;
  arma::mat P(n,1);
  // for (j = n; j--;) p(j,0) = j+1;
  arma::mat g(n,1, fill::zeros);
  // arma::mat E(n,1, fill::zeros);
  E = mat(n, 1, fill::zeros);
  double gamma = A(n-1,n-1);
  if (gamma < 0) phase1 = false;
  for (j = 0; j < n-1; j++){
    if (A(j, j) < 0) phase1 = false;
    if (A(j, j) > gamma) gamma = A(j,j);
  }
  double taugam = tau1*gamma;
  if (!phase1) g = gershNested(A, 0, n);
  // N=1 case
  if (n == 1){
    delta = tau2*std::fabs(A(0,0)) - A(0,0);
    if (delta > 0) E(0,0) = delta;
    if (A(0,0) == 0) E(0,0) = tau2;
    A(0,0)=_safe_sqrt(A(0,0)+E(0,0));
    Ao = A;
    return true;
  }
  int jp1, ii, k;
  double tempjj, temp=1., normj, tmp;
  for (j = 0; j < n-1;  j++){
    // Pivoting not included
    if (phase1){
      jp1 = j+1;
      if (A(j,j)>0){
	arma::mat tmp = (A(span(jp1,n-1),span(jp1,n-1))).diag() - A(span(jp1, n-1),j)%A(span(jp1, n-1),j)/A(j,j);
	double mintmp = tmp[0];
	for (ii = 1; ii < (int)tmp.size(); ii++) mintmp = (mintmp < tmp[ii]) ? mintmp : tmp[ii];
	if (mintmp < taugam) phase1=false;
      } else phase1 = false;

      if (phase1){
	// Do the normal cholesky update if still in phase 1
	A(j,j) = _safe_sqrt(A(j,j));
	tempjj = A(j,j);
	for (ii = jp1; ii < n; ii++){
	  A(ii,j) = A(ii,j)/tempjj;
	}
	for (ii=jp1; ii <n; ii++){
	  temp=A(ii,j);
	  for (k = jp1; k < ii+1; k++){
	    A(ii,k) = A(ii,k)-(temp * A(k,j));
	  }
	}
	if (j == n-2){
	  A(n-1,n-1)=_safe_sqrt(A(n-1,n-1));
	}
      } else {
	// Calculate the negatives of the lower Gershgorin bounds
	g=gershNested(A,j,n);
      }
    }

    if (!phase1){
      if (j != n-2){
        // Calculate delta and add to the diagonal. delta=max{0,-A(j,j) + max{normj,taugam},delta_previous}
	// where normj=sum of |A(i,j)|,for i=1,n, delta_previous is the delta computed at the previous iter and taugam is tau1*gamma.
	normj=arma::sum(arma::abs(A(span(j+1, n-1),j)));
	if (delta < 0) delta = 0;
	tmp  = -A(j,j)+normj;
	if (delta < tmp) delta = tmp;
        tmp  = -A(j,j)+taugam;
        if (delta < tmp) delta = tmp;
	// get adjustment based on formula on bottom of p. 309 of Eskow/Schnabel (1991)
	E(j,0) =  delta;
	A(j,j) = A(j,j) + E(j,0);
	// Update the Gershgorin bound estimates (note: g(i) is the negative of the Gershgorin lower bound.)
	if (A(j,j) != normj){
	  temp = (normj/A(j,j)) - 1;
	  for (ii = j+1; ii < n; ii++){
	    g(ii) = g(ii) + std::fabs(A(ii,j)) * temp;
	  }
	}
	for (int ii = j+1; ii < n; ii++){
	  g(ii,0) = g(ii,0) + std::fabs(A(ii,j)) * temp;
	}
	// Do the cholesky update
	A(j,j) = _safe_sqrt(A(j,j));
	tempjj = A(j,j);
	for (ii = j+1; ii < n; ii++){
	  A(ii,j) = A(ii,j) / tempjj;
	}
	for (ii = j+1; ii < n; ii++){
	  temp = A(ii,j);
	  for (k = j+1; k < ii+1; k++){
	    A(ii,k) = A(ii,k) - (temp * A(k,j));
	  }
	}
      } else {
	// Find eigenvalues of final 2 by 2 submatrix
        // Find delta such that:
	// 1.  the l2 condition number of the final 2X2 submatrix + delta*I <= tau2
	// 2. delta >= previous delta,
	// 3. min(eigvals) + delta >= tau2 * gamma, where min(eigvals) is the smallest eigenvalue of the final 2X2 submatrix
	// A(n-2,n-1)=A(n-1,n-2);
	//set value above diagonal for computation of eigenvalues
        A(n-2,n-1)=A(n-1,n-2); //set value above diagonal for computation of eigenvalues
	vec eigvals  = eig_sym(A(span(n-2, n-1),span(n-2, n-1)));
        // Formula 5.3.2 of Schnabel/Eskow (1990)
	if (delta < 0) delta = 0;
	tmp= (max(eigvals)-min(eigvals))/(1-tau1);
	if (tmp < gamma) tmp = gamma;
	tmp=tau2*tmp;
	tmp =tmp - min(eigvals);
	if (delta < tmp) delta = tmp;
	if (delta > 0){
	  A(n-2, n-2) = A(n-2,n-2) + delta;
	  A(n-1, n-1) = A(n-1,n-1) + delta;
          E(n-2, 0) = delta;
          E(n-1, 0) = delta;
        }
	// Final update
	A(n-2,n-2) = _safe_sqrt(A(n-2,n-2));
        A(n-1,n-2) = A(n-1,n-2)/A(n-2,n-2);
        A(n-1,n-1) = A(n-1,n-1) - A(n-1,n-2)*A(n-1,n-2);
        A(n-1,n-1) = _safe_sqrt(A(n-1,n-1));
      }
    }
  }
  Ao = (trimatl(A)).t();
  return phase1;
}

arma::mat cholSE__(arma::mat A, double tol) {
  arma::mat Ao, E;
  cholSE0(Ao, E, A, tol);
  return Ao;
}
//[[Rcpp::export]]
NumericMatrix cholSE_(NumericMatrix A, double tol){
  arma::mat Ao, E;
  cholSE0(Ao, E, as<arma::mat>(A), tol);
  return wrap(Ao);
}
