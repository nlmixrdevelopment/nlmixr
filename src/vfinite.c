#include <R.h>
#include <Rinternals.h>

// modified from scilab to use R interface.  
int F77_SUB(vfinite)(int *n, double *v)
{
    int i;
    for (i = 0; i < *n; i++)
      if (!R_FINITE(v[i])) {
            return 0;
      }
    return 1;
}
