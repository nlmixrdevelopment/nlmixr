// neldermead.cpp: population PK/PD modeling library
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

#include <cstdlib>
#include <math.h>
#include <R.h>

typedef void (*fn_ptr) (double *, double *);
#define MXPAR  45

void nelder_(fn_ptr func, int n, double *start, double *step,
  int itmax, double ftol_rel, double rcoef, double ecoef, double ccoef,
  int *iconv, int *it, int *nfcall, double *ynewlo, double *xmin,
  int *iprint)
{
  double fval;
  int i, j, k;
  int ihi, ilo, nn, ibest=0;
  int konvge, kcount;
  double *p, *pbar, *pstar, *p2star, *y;
  double ystar, xlo, xhi, y2star, bignum;
  double ylo, dchk, z, dn, dabit, yoldlo=0;

  //FIXME: check malloc status
  p = (double *) malloc(n*(n+1)*sizeof(double));
  pstar = (double *) malloc(n*sizeof(double));
  p2star = (double *) malloc(n*sizeof(double));
  pbar = (double *) malloc(n*sizeof(double));
  y = (double *) malloc((n+1)*sizeof(double));

  kcount = 1000000;
  *nfcall = 0;
  *it = 0;
  *iconv = 0;

  /* check inputs */
  if (n <= 0 || n > MXPAR) *nfcall += -10;
  if (*nfcall < 0) return;

  /* constants */
  dabit = 2.2204460492503131e-16;
  bignum = 1e38;
  konvge = 5;
  dn = (double) (n);
  nn = n + 1;

  /* initial simplex */
  for (i = 0; i < n; ++i)
    p[i+n*n] = start[i];

  (*func)(start, &fval);
  y[n] = fval;
  ++(*nfcall);

  if (itmax == 0) {
    for (i = 0; i < n; ++i)
      xmin[i] = start[i];

    *ynewlo = fval;
    return;
  }
  else {
    for (j = 0; j < n; ++j) {
      dchk = start[j];
      start[j] = dchk + step[j];
      for (i = 0; i < n; ++i)
        p[i+j*n] = start[i];
      (*func)(start, &fval);  //CHECKME: start or &start[?]
      y[j] = fval;
      ++(*nfcall);
      start[j] = dchk;
    }
  }

  /* hi/lo values */
  for(;;)
  {
    ylo = y[0];
    *ynewlo = ylo;
    ilo = 0;
    ihi = 0;
    for (i = 1; i < nn; ++i) {
      if (y[i] < ylo) {
        ylo = y[i];
        ilo = i;
      }
      if (y[i] > *ynewlo) {
        *ynewlo = y[i];
        ihi = i;
      }
    }

    if (*nfcall <= nn) yoldlo = ylo;
    else if (ylo < yoldlo) {
      yoldlo = ylo;
      ++(*it);
      if (*it >= itmax) break;
    }

    if (*iprint) Rprintf("%d %d %f\n", *it, *nfcall, ylo);

    /* convergence checks*/
    dchk = (*ynewlo + dabit) / (ylo + dabit) - 1.;
    if (fabs(dchk) <= ftol_rel) {
      *iconv = 1;
      break;
    }

    --konvge;
    if (konvge == 0) {
      konvge = 5;
      *iconv = 2;

      /* convergence check of coordinates */
      for (i = 0; i < n; ++i) {
        xlo = p[i];
        xhi = xlo;

        for (j = 1; j < nn; ++j) {
          if (p[i+j*n] < xlo) xlo = p[i+j*n];
          if (p[i+j*n] > xhi) xhi = p[i+j*n];
        }

        dchk = (xhi + dabit) / (xlo + dabit) - 1.;
        if (fabs(dchk) > ftol_rel) {
         *iconv = 0;
         break;
        }
      }

      if (*iconv) {
#ifdef __DEBUG__
    	printf("fval chk: %21.16f %21.16f %21.16f\n", *ynewlo, ylo, (*ynewlo + dabit) / (ylo + dabit));
#endif
	    break;
	  }
    }
    if (*nfcall >= kcount) break;

    /* calculate centroid of simplex */
    for (i = 0; i < n; ++i) {
      z = 0.;
      for (j = 0; j < nn; ++j)
        z += p[i+j*n];
      z -= p[i+ihi*n];
      pbar[i] = z / dn;
    }

    /* reflection */
    for (i = 0; i < n; ++i)
      pstar[i] = pbar[i] + rcoef*(pbar[i] - p[i+ihi*n]);
    (*func)(pstar, &fval);
    ystar = fval;
    ++(*nfcall);

    if (ystar < ylo) {
      if (*nfcall < kcount) {
        for (i = 0; i < n; ++i)
          p2star[i] = pbar[i] + ecoef*(pstar[i] - pbar[i]);
        (*func)(p2star, &fval);
        y2star = fval;
        ++(*nfcall);

        if (y2star < ystar) {
          for (i = 0; i < n; ++i)
            p[i+ihi*n] = p2star[i];
          y[ihi] = y2star;
          continue;
        }
        else {
          for (i = 0; i < n; ++i)
            p[i+ihi*n] = pstar[i];
          y[ihi] = ystar;
          continue;
        }
      }
      else {
        for (i = 0; i < n; ++i)
          p[i+ihi*n] = pstar[i];
        y[ihi] = ystar;
        continue;
      }
    }
    else {
      k = 0;
      for (i = 0; i < nn; ++i) {
        if (y[i] > ystar) ++k;
      }
      if (k > 1) {
        for (i = 0; i < n; ++i)
          p[i+ihi*n] = pstar[i];
        y[ihi] = ystar;
        continue;
      }
      if (k != 0) {
        for (i = 0; i < n; ++i) {
          p[i+ihi*n] = pstar[i];
        }
        y[ihi] = ystar;
      }
    }

    /* contraction */
    if (*nfcall >= kcount) break;
    for (i = 0; i < n; ++i)
      p2star[i] = pbar[i] + ccoef*(p[i+ihi*n] - pbar[i]);
    (*func)(p2star, &fval);
    y2star = fval;
    ++(*nfcall);

    if (y2star < y[ihi]) {
      for (i = 0; i < n; ++i)
        p[i+ihi*n] = p2star[i];
      y[ihi] = y2star;
      continue;
    }

    for (j = 0; j < nn; ++j) {
      for (i = 0; i < n; ++i) {
        p[i+j*n] = (p[i+j*n] + p[i+ilo*n])*.5;
        xmin[i] = p[i+j*n];
      }
      (*func)(xmin, &fval);
      y[j] = fval;
    }
    *nfcall += nn;

    if (*nfcall >= kcount) break;

  }

  /* get optimum & optima */
  for (j = 0; j < nn; ++j) {
    for (i = 0; i < n; ++i)
      xmin[i] = p[i+j*n];
    (*func)(xmin, &fval);
    y[j] = fval;
  }
  *nfcall += nn;

  *ynewlo = bignum;
  for (j = 0; j < nn; ++j) {
    if (y[j] < *ynewlo) {
      *ynewlo = y[j];
      ibest = j;
    }
  }
  for (i = 0; i < n; ++i)
    xmin[i] = p[i+ibest*n];
  return;
}
