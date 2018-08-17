// neldermead.hpp: population PK/PD modeling library
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


#ifndef __NELDERMEAD_HPP__
#define __NELDERMEAD_HPP__
#include <R_ext/Rdynload.h>

// for resid error

typedef void (*fn_ptr) (double *, double *);
typedef void (*nelder_fn)(fn_ptr func, int n, double *start, double *step,
                                   int itmax, double ftol_rel, double rcoef, double ecoef, double ccoef,
                                   int *iconv, int *it, int *nfcall, double *ynewlo, double *xmin,
                                   int *iprint);
void nelder_ (fn_ptr func, int n, double *start, double *step,
		       int itmax, double ftol_rel, double rcoef, double ecoef, double ccoef,
		       int *iconv, int *it, int *nfcall, double *ynewlo, double *xmin,
		       int *iprint){
    static nelder_fn fun=NULL;
    if (fun == NULL) fun = (nelder_fn) R_GetCCallable("nlmixr","nelder_fn");
    fun(func, n, start, step, itmax, ftol_rel, rcoef, ecoef, ccoef, iconv, it, nfcall, ynewlo, xmin, iprint);
}

double *yptr, *fptr;	//CHK
int len;	//CHK

void obj(double *ab, double *fx)
{
	int i;
	double g, sum;
    double xmin = 1.0e-200;

	for (i=0, sum=0; i<len; ++i) {
        // nelder_() does not allow lower bounds; we force ab[] be positive here
		g = ab[0]*ab[0] + ab[1]*ab[1]*fabs(fptr[i]);
		if (g < xmin) g = xmin;
		sum += pow((yptr[i]-fptr[i])/g, 2.0) + 2*log(g);
	}

	*fx = sum;
}

#endif
