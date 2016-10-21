// slice.cpp: population PK/PD modeling library
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


#include <math.h>
#include <stdlib.h>
#include <Rcpp.h>

//#define genunf(a,b) Rcpp::runif(1, (a), (b))[0]
#define genunf(a,b) (a) + ((b) - (a))*R::unif_rand()

double rexp(double beta)
{
	return -log(genunf(0,1))/beta;
}

/*--
C TRANSLATION OF
R FUNCTION FOR PERFORMING UNIVARIATE SLICE SAMPLING.

Radford M. Neal, 17 March 2008.

Implements, with slight modifications and extensions, the algorithm described
in Figures 3 and 5 of the following paper:

  Neal, R. M (2003) "Slice sampling" (with discussion), Annals of Statistics,
     vol. 31, no. 3, pp. 705-767.
--*/


// GLOBAL VARIABLES FOR RECORDING PERFORMANCE.
static int uni_slice_calls = 0;	// Number of calls of the slice sampling function
static int uni_slice_evals = 0;	// Number of density evaluations done in these calls

/*--
UNIVARIATE SLICE SAMPLING WITH STEPPING OUT AND SHRINKAGE.

Performs a slice sampling update from an initial point to a new point that
leaves invariant the distribution with the specified log density function.

Arguments:

  x0    Initial point
  g     Function returning the log of the probability density (plus constant)
  w     Size of the steps for creating interval (default 1)
  m     Limit on steps (default infinite)
  lower Lower bound on support of the distribution (default -Inf)
  upper Upper bound on support of the distribution (default +Inf)

The log density function may return -Inf for points outside the support
of the distribution.  If a lower and/or upper bound is specified for the
support, the log density function will not be called outside such limits.

The value of this function is the new point sampled, with an attribute
of "log.density" giving the value of the log density function, g, at this
point.  Depending on the context, this log density might be passed as the
gx0 argument of a future call of uni.slice.

The global variable uni_slice_calls is incremented by one for each call
of uni.slice.  The global variable uni_slice_evals is incremented by the
number of calls made to the g function passed.

WARNING:  If you provide a value for g(x0), it must of course be correct!
In addition to giving wrong answers, wrong values for gx0 may result in
the uni.slice function going into an infinite loop.
--*/

double uni_slice(double x0, double (*g)(double), double w, int m, double lower, double upper)
{
	int J, K;
	double u, L, R, logy, gx0, x1, gx1;

	// Keep track of the number of calls made to this function.
	uni_slice_calls ++;

	// Find the log density at the initial point, if not already known.
	uni_slice_evals ++;
    gx0 = g(x0);

	// Determine the slice level, in log terms.
	logy = gx0 - rexp(1);

	// Find the initial interval to sample from.
	u = genunf(0,w);
	L = x0 - u;
	R = x0 + (w-u);  // should guarantee that x0 is in [L,R], even with roundoff

	// Expand the interval until its ends are outside the slice, or until
	// the limit on steps is reached.

	if (m<=0)  // no limit on number of steps
	{
		for(;;)
		{
			if (L<=lower) break;
			uni_slice_evals ++;
			if (g(L)<=logy) break;
			L -= w;
		}
		for(;;)
		{
			if (R>=upper) break;
			uni_slice_evals ++;
			if (g(R)<=logy) break;
			R += w;
		}
	}
	else if (m>1)  // limit on steps, bigger than one
	{
		J = floor(genunf(0,m));
		K = (m-1) - J;

		while (J>0)
		{
			if (L<=lower) break;
			uni_slice_evals ++;
			if (g(L)<=logy) break;
			L -= w;
			J--;
		}

		while (K>0)
		{
			if (R>=upper) break;
			uni_slice_evals ++;
			if (g(R)<=logy) break;
			R += w;
			K --;
		}
	}

	// Shrink interval to lower and upper bounds.
	if (L<lower) L = lower;
	if (R>upper) R = upper;

	// Sample from the interval, shrinking it on each rejection.
	for(;;)
	{
		x1 = genunf(L,R);

		uni_slice_evals ++;
		gx1 = g(x1);

		if (gx1>=logy) break;

		if (x1>x0) R = x1;
		else L = x1;
	}

	// Return the point sampled, with its log density attached as an attribute.
	return (x1);
}


/*

x = rnorm(20, 1, 1)
fr <- function(m) {
    sum(log(dnorm(x, m, 1)))
}
rho = environment(fr)

require(Rcpp)
dyn.load("slice.dll")
is.loaded("slice_wrap")

nsim = 200; o = rep(NA, nsim)
x0 = 1
for (i in 1:nsim)
{
	x1 = .Call("slice_wrap", fr, rho, x0)$x1
	x0 = x1
	o[i] = x1
}

require(lattice)
histogram(~o)
*/
