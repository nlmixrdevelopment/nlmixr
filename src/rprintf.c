#define STRICT_R_HEADER
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// modified from scilab to use R interface.  Assumes io is alwyas Rprintf; half of it was cut.

int F77_SUB(basout)(int *io, int *lunit, char *string,long int nbcharacters){
  /* bug 3831 */
  if (string) {
      int i = 0;
      for (i = 0; i < nbcharacters - 1; i++)  {
          if (string[i] == 0) string[i] = ' ';
      }
  }
  if (string) {
    if (nbcharacters > 1) {
          /* on linux , q=[] crashs with previous version 
             in printf.f line 102 
             call basout(io,lunit,'     []')
             if we do basout(io,lunit,'     []',7) it works ...
             temp workaround , we returns to old version with a allocation
          */
	  char *buffer = (char *)R_Calloc(nbcharacters+1,char);
          if (buffer) {
              strncpy(buffer,string,nbcharacters);
              buffer[nbcharacters]='\0';
              Rprintf("%s\n",buffer);
              R_Free(buffer);
	  } else {
              Rprintf("\n");
            }
    } else if (nbcharacters == 1) {
          Rprintf("%c\n", string[0]);
    } else {
          Rprintf("\n");
    }
  } else Rprintf("\n");
  return 0;
} 
/*--------------------------------------------------------------------------*/ 

void F77_SUB(rprintf)(char* msg) {
  Rprintf(msg);
  Rprintf("\n");
}

void F77_SUB(rprintflen)(char* msg, int *i) {
  Rprintf("%.%s,",*i,msg);
  Rprintf("\n");
}

// may be redundant
void F77_SUB(rprintf2)(char* msg) {
  Rprintf(msg);
  Rprintf("\n");
}



void F77_SUB(rprintfid)(char* msg, int *i, double *d) {
  Rprintf(msg, *i, *d);
  Rprintf("\n");
}

void F77_SUB(rprintfdi)(char* msg, double *d, int *i) {
  Rprintf(msg, *d, *i);
  Rprintf("\n");
}

void F77_SUB(rprintfdid)(char* msg, double *d1, int *i, double *d2) {
  Rprintf(msg, *d1, *i, *d2);
  Rprintf("\n");
}

void F77_SUB(rprintfd1)(char* msg, double *d) {
  Rprintf(msg, *d);
  Rprintf("\n");
}

void F77_SUB(rprintfd2)(char* msg, double *d1, double *d2) {
  Rprintf(msg, *d1, *d2);
  Rprintf("\n");
}

void F77_SUB(rprintfi1)(char* msg, int *i) {
  Rprintf(msg, *i);
  Rprintf("\n");
}

void F77_SUB(rprintfi2)(char* msg, int *i1, int *i2) {
  Rprintf(msg, *i1, *i2);
  Rprintf("\n");
}

void F77_SUB(rprintfi3)(char* msg, int *i1, int *i2, int* i3) {
  Rprintf(msg, *i1, *i2, *i3);
  Rprintf("\n");
}
