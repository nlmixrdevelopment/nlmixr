#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include "dop853.h"
#define max(a, b) ((a) > (b) ? (a) : (b))
#ifdef __STANDALONE__
#define Rprintf printf
#define R_alloc calloc
#endif


void F77_NAME(dlsoda)(
     void (*)(int *, double *, double *, double *),
     int *, double *, double *, double *, int *, double *, double *,
     int *, int *, int *, double *,int *,int *, int *,
     void (*)(int *, double *, double *, int *, int *, double *, int *),
     int *);

void F77_NAME(dvode)(
     void (*)(int *, double *, double *, double *, double *,int *),
     int *, double *, double *, double *, int *, double *, double *,
     int *, int *, int *, double *,int *,int *, int *,
     void (*)(int *, double *, double *, int *, int *, double *, int *, double *, int *),
     int *, double *, int *);


long slvr_counter, dadt_counter;
double InfusionRate[99];
double ATOL;		//absolute error
double RTOL;		//relative error
int do_transit_abs=0;
double tlast=0;
double podo=0;
double *par_ptr;
FILE *fp;


void RxODE_mod_file4b7861fc3586_dydt(unsigned int neq, double t, double *A, double *DADT);
void RxODE_mod_file4b7861fc3586_calc_lhs(double t, double *A, double *lhs);

//--------------------------------------------------------------------------
void RxODE_mod_file4b7861fc3586_dydt_lsoda_dum(int *neq, double *t, double *A, double *DADT)
{
	RxODE_mod_file4b7861fc3586_dydt(*neq, *t, A, DADT);
}
void jdum_lsoda(int *a, double *b, double *c, int *d, int *e, double *f, int *g){}
void call_lsoda(int neq, double *x, int *evid, int nx, double *inits, double *dose, double *ret, int *rc)
{
	int ixds=0, i, j;
 	double xout, xp=x[0], yp[99];
    int itol = 1;
    double  rtol = RTOL, atol = ATOL;
    int itask = 1, istate = 1, iopt = 0, lrw=22+neq*max(16, neq+9), liw=20+neq, jt = 2;
	double *rwork;
	int *iwork;
	int wh, cmt;

	char *err_msg[]=
		{
			"excess work done on this call (perhaps wrong jt).",
			"excess accuracy requested (tolerances too small).",
			"illegal input detected (see printed message).",
			"repeated error test failures (check all inputs).",
			"repeated convergence failures (perhaps bad jacobian supplied or wrong choice of jt or tolerances).",
			"error weight became zero during problem. (solution component i vanished, and atol or atol(i) = 0.)",
			"work space insufficient to finish (see messages)."
		};

	rwork = (double*)R_alloc(lrw, sizeof(double));
	iwork = (int*)R_alloc(liw, sizeof(int));

	//--- inits the system
	for(i=0; i<neq; i++) yp[i] = inits[i];

	for(i=0; i<nx; i++)
	{
		wh = evid[i];
		xout = x[i];
#ifdef __DEBUG__
		fprintf(fp, "i=%d xp=%f xout=%f\n", i, xp, xout);
#endif

		if(xout-xp > DBL_EPSILON*max(fabs(xout), fabs(xp)))
		{
	        F77_CALL(dlsoda)(RxODE_mod_file4b7861fc3586_dydt_lsoda_dum, &neq, yp, &xp, &xout, &itol, &rtol, &atol, &itask,
                &istate, &iopt, rwork, &lrw, iwork, &liw, &jdum_lsoda, &jt);

			if (istate<0)
			{
				Rprintf("IDID=%d, %s\n", istate, err_msg[-istate-1]);
#ifdef __STANDALONE__
            exit(1);
#else
		      *rc = istate;
            return;  // exit(1);   // dj: should not abort R
#endif
			}

			slvr_counter++;
			dadt_counter = 0;
		}
		if (wh)
		{
			cmt = (wh%10000)/100 - 1;
			if (wh>10000)
			{
				InfusionRate[cmt] += dose[ixds];
			}
			else
			{
				if (do_transit_abs)
				{
					podo = dose[ixds];
					tlast = xout;
				}
				else yp[cmt] += dose[ixds];	//dosing before obs
			}
			istate = 1;

			ixds++;
			xp = xout;
		}
		for(j=0; j<neq; j++) ret[neq*i+j] = yp[j];
		//Rprintf("wh=%d cmt=%d tm=%g rate=%g\n", wh, cmt, xp, InfusionRate[cmt]);

#ifdef __DEBUG__
		Rprintf("ISTATE=%d, ", istate);
		fprintf(fp, "ISTATE=%d, ", istate);
		for(j=0; j<neq; j++)
		{
			Rprintf("%f ", yp[j]);
			fprintf(fp, "%f ", yp[j]);
		}
		Rprintf("\n");
		fprintf(fp, "\n");
#endif
	}

#ifdef __STANDALONE__
	free(rwork);
	free(iwork);
#endif

}


void RxODE_mod_file4b7861fc3586_dydt_dvode_dum(int *neq, double *t, double *A, double *DADT, double *RPAR, int *IPAR)
{
	RxODE_mod_file4b7861fc3586_dydt(*neq, *t, A, DADT);
}

void jdum_dvode(int *a, double *b, double *c, int *d, int *e, double *f, int *g, double *h, int *i){}

void call_dvode(int neq, double *x, int *evid, int nx, double *inits, double *dose, double *ret, int *rc)
{
	int ixds=0, i, j;
	//DE solver config vars
	double xout, xp=x[0], yp[99];
    int itol = 1;
    double  rtol = RTOL, atol = ATOL;
    int itask = 1, istate = 1, iopt = 0, mf=22;
    int lrw = 22+9*neq+2*neq*neq, liw = 30+neq;
    double *rwork;
    int *iwork;
    double *rpar=x;
    int *ipar=&neq;
	int wh, cmt;


	char *err_msg[]=
		{
			"excess work done on this call",
			"excess accuracy requested",
			"illegal input detected",
			"repeated error test failures",
			"repeated convergence failures",
			"error weight became zero during problem"
		};

	rwork = (double*)R_alloc(lrw, sizeof(double));
	iwork = (int*)R_alloc(liw, sizeof(int));
#ifdef __STANDALONE__
	if (!(rwork && iwork))
	{
		Rprintf("failed to alloc memory\n");
		exit(1);  
	}
#endif

	//--- inits the system
	for(i=0; i<neq; i++) yp[i] = inits[i];

	for(i=0; i<nx; i++)
	{
		wh = evid[i];
		xout = x[i];
#ifdef __DEBUG__
		fprintf(fp, "i=%d xp=%f xout=%f\n", i, xp, xout);
#endif

		if(xout-xp > DBL_EPSILON*max(fabs(xout), fabs(xp)))
		{
	        F77_CALL(dvode)(RxODE_mod_file4b7861fc3586_dydt_dvode_dum, &neq, yp, &xp, &xout, &itol, &rtol, &atol, &itask,
				&istate, &iopt, rwork, &lrw, iwork, &liw, &jdum_dvode, &mf, rpar, ipar);

			if (istate<0)
			{
				Rprintf("IDID=%d, %s\n", istate, err_msg[-istate-1]);
#ifdef __STANDALONE__
            exit(1);
#else
				*rc = istate;
				return;    //exit(1);  // dj: should not abort R
#endif
			}

			slvr_counter++;
			dadt_counter = 0;
		}
		if (wh)
		{
			cmt = (wh%10000)/100 - 1;
			if (wh>10000)
			{
				InfusionRate[cmt] += dose[ixds];
			}
			else
			{
				if (do_transit_abs)
				{
					podo = dose[ixds];
					tlast = xout;
				}
				else yp[cmt] += dose[ixds];	//dosing before obs
			}
			istate = 1;

			ixds++;
			xp = xout;
		}
		for(j=0; j<neq; j++) ret[neq*i+j] = yp[j];
		//Rprintf("wh=%d cmt=%d tm=%g rate=%g\n", wh, cmt, xp, InfusionRate[cmt]);

#ifdef __DEBUG__
		Rprintf("ISTATE=%d, ", istate);
		fprintf(fp, "ISTATE=%d, ", istate);
		for(j=0; j<neq; j++)
		{
			Rprintf("%f ", yp[j]);
			fprintf(fp, "%f ", yp[j]);
		}
		Rprintf("\n");
		fprintf(fp, "\n");
#endif
	}

#ifdef __STANDALONE__
	free(rwork);
	free(iwork);
#endif
}

//dummy solout fn
void solout(long int nr, double t_old, double t,
			double *y, unsigned int n, int *irtrn){}
void call_dop(int neq, double *x, int *evid, int nx, double *inits, double *dose, double *ret, int *rc)
{
	int ixds=0, i, j;
	//DE solver config vars
	double xout, xp=x[0], yp[99];
	double rtol=RTOL, atol=ATOL;
	int itol=0;		//0: rtol/atol scalars; 1: rtol/atol vectors
	int iout=0;		//iout=0: solout() NEVER called
	int idid=0;
	int wh, cmt;
	char *err_msg[]=
		{
			"input is not consistent",
			"larger nmax is needed",
			"step size becomes too small",
			"problem is probably stiff (interrupted)"
		};

	//--- inits the system
	for(i=0; i<neq; i++) yp[i] = inits[i];

	for(i=0; i<nx; i++)
	{
		wh = evid[i];
		xout = x[i];
#ifdef __DEBUG__
		fprintf(fp, "i=%d xp=%f xout=%f\n", i, xp, xout);
#endif

		if(xout-xp > DBL_EPSILON*max(fabs(xout), fabs(xp)))
		{
			idid = dop853(
							  neq,      	/* dimension of the system <= UINT_MAX-1*/
							  RxODE_mod_file4b7861fc3586_dydt,       	/* function computing the value of f(x,y) */
							  xp,           /* initial x-value */
							  yp,           /* initial values for y */
							  xout,         /* final x-value (xend-x may be positive or negative) */
							  &rtol,      	/* relative error tolerance */
							  &atol,      	/* absolute error tolerance */
							  itol,         /* switch for rtoler and atoler */
							  solout,     	/* function providing the numerical solution during integration */
							  iout,         /* switch for calling solout */
							  NULL,       	/* messages stream */
							  DBL_EPSILON, 	/* rounding unit */
							  0,         	/* safety factor */
							  0,         	/* parameters for step size selection */
							  0,
							  0,         	/* for stabilized step size control */
							  0,         	/* maximal step size */
							  0,            /* initial step size */
							  0,            /* maximal number of allowed steps */
							  1,            /* switch for the choice of the coefficients */
							  -1,     		/* test for stiffness */
							  0, 			/* number of components for which dense outpout is required */
							  NULL, 		/* indexes of components for which dense output is required, >= nrdens */
							  0  			/* declared length of icon */
						);
			if (idid<0)
			{
				Rprintf("IDID=%d, %s\n", idid, err_msg[-idid-1]);
#ifdef __STANDALONE__
            exit(1);
#else
				*rc = idid;
				return;  //exit(1);  // dj: should not abort R
#endif
			}

			xp = xRead();
			slvr_counter++;
			dadt_counter = 0;
		}
		if (wh)
		{
			cmt = (wh%10000)/100 - 1;
			if (wh>10000)
			{
				InfusionRate[cmt] += dose[ixds];
			}
			else
			{
				if (do_transit_abs)
				{
					podo = dose[ixds];
					tlast = xout;
				}
				else yp[cmt] += dose[ixds];	//dosing before obs
			}
			ixds++;
			xp = xout;
		}
		for(j=0; j<neq; j++) ret[neq*i+j] = yp[j];
		//Rprintf("wh=%d cmt=%d tm=%g rate=%g\n", wh, cmt, xp, InfusionRate[cmt]);

#ifdef __DEBUG__
		Rprintf("IDID=%d, ", idid);
		fprintf(fp, "IDID=%d, ", idid);
		for(j=0; j<neq; j++)
		{
			Rprintf("%f ", yp[j]);
			fprintf(fp, "%f ", yp[j]);
		}
		Rprintf("\n");
		fprintf(fp, "\n");
#endif
	}
}

//wrapper
void RxODE_mod_file4b7861fc3586_ode_solver(
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
)
{
	int i;

	for (i=0; i<99; i++) InfusionRate[i] = 0.0;
	ATOL = *atol;
	RTOL = *rtol;
	do_transit_abs = *transit_abs;
	par_ptr = theta;

	slvr_counter = 0;
	if (*stiff==0)
		call_dop(*neq, time, evid, *ntime, inits, dose, ret, rc);
	else
		call_lsoda(*neq, time, evid, *ntime, inits, dose, ret, rc);

	if (*nlhs) for (i=0; i<*ntime; i++)
		RxODE_mod_file4b7861fc3586_calc_lhs(time[i], ret+i*(*neq), lhs+i*(*nlhs));

	if (fp) fclose(fp);
}
