#include <R_ext/Applic.h>
// This is to avoid BLAS conflicts with arma.
void lbfgsbRX(int n, int lmm, double *x, double *lower,
	      double *upper, int *nbd, double *Fmin, optimfn fn,
	      optimgr gr, int *fail, void *ex, double factr,
	      double pgtol, int *fncount, int *grcount,
	      int maxit, char *msg, int trace, int nREPORT){
  /*
  From:  src/appl/lbfgsb.c

      n is an integer variable.
	 On entry n is the number of variables.
	 On exit n is unchanged.

      lmm is an integer giving the number of BFGS 
         updates retained in the "L-BFGS-B" method, It defaults to 5.


      x is a double precision array of dimension n.
	 On entry x is an approximation to the solution.
	 On exit x is the current approximation.

      lower, upper 
         Bounds on the variables for the "L-BFGS-B" method
 

      nbd is an integer array of dimension n.
	 On entry nbd represents the type of bounds imposed on the
	   variables, and must be specified as follows:
	   nbd(i)=0 if x(i) is unbounded,
		  1 if x(i) has only a lower bound,
		  2 if x(i) has both lower and upper bounds,
		  3 if x(i) has only an upper bound.
	 On exit nbd is unchanged.

     Fmin pointer to the minimum value.

     fn function

     gr gradient

     fail pointer to integer to give fail/success code.
    
     pgtol is a double precision variable.
	 On entry pgtol >= 0 is specified by the user.	The iteration
	   will stop when
		   max{|proj g_i | i = 1, ..., n} <= pgtol
	   where pg_i is the ith component of the projected gradient.
	 On exit pgtol is unchanged.
         This defaults to zero, when the check is suppressed.

     *ex -- Extra information for minimization
     
     factr -- controls the convergence of the "L-BFGS-B" method. 
       Convergence occurs when the reduction in the objective is 
       within this factor of the machine tolerance. Default is 1e7, 
       that is a tolerance of about 1e-8.

     fncount -- pointer to function count

     grcount -- pointer to gradient count

     maxit -- Maximum number of iterations

     msg -- Minimization message

     trace -- Non-negative integer. If positive, tracing information 
         on the progress of the optimization is produced. Higher values 
         may produce more tracing information: for method "L-BFGS-B" there 
         are six levels of tracing. 

     nREPORT -- Number of iterations before printout.


  */
  lbfgsb(n, lmm, x, lower, upper, nbd, Fmin, fn,
	 gr, fail, ex, factr, pgtol, fncount, grcount, 
	 maxit, msg, trace, nREPORT);
}
