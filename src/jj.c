// jj.c: population PK/PD modeling library
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

#ifndef __STANDALONE__
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>
#define printf Rprintf
#else
#include <stdio.h>
#endif
#include <dparser.h>


#define max(a,b) (a)>(b) ? (a):(b)
#define MXSYM 5000
#define MXDER 500
#define MXLEN 1200
#define MXBUF 2400
#define SBPTR sb.s+sb.o

// Taken from dparser and changed to use R_alloc
int r_buf_read(const char *pathname, char **buf, int *len) {
  struct stat sb;
  int fd;

  *buf = 0;
  *len = 0;
  fd = open(pathname, O_RDONLY);
  if (fd <= 0) 
    return -1;
  memset(&sb, 0, sizeof(sb));
  fstat(fd, &sb);
  *len = sb.st_size;
  *buf = (char*)R_alloc(*len + 2,sizeof(char));
  // MINGW likes to convert cr lf => lf which messes with the size
  size_t real_size = read(fd, *buf, *len);
  (*buf)[real_size] = 0;
  (*buf)[real_size + 1] = 0;
  *len = real_size;
  close(fd);
  return *len;
}

// Taken from dparser and changed to use R_alloc
char * r_sbuf_read(const char *pathname) {
  char *buf;
  int len;
  if (r_buf_read(pathname, &buf, &len) < 0)
    return NULL;
  return buf;
}


// Taken from dparser and changed to use Calloc
char * rc_dup_str(const char *s, const char *e) {
  int l = e ? e-s : strlen(s);
  char *ss = Calloc(l+1,char);
  memcpy(ss, s, l);
  ss[l] = 0;
  return ss;
}


typedef void (print_node_fn_t)(int depth, char *token_name, char *token_value, void *client_data);

extern D_ParserTables parser_tables_gram;


char *ignore_symb_str;
int *ignore_offset;
int nignore;

typedef struct Symtable {
  char *symb_str;       /* symbol string: all vars*/
  int offset[MXSYM];    /* offset of symbols */
  int is_lhs[MXSYM];    /* lhs symbols? =9 if a state var*/
  int der_ix[MXDER];    /* ith of state vars */
  int nvar;             /* nbr of symbols */
  int nder;             /* nbr of dydt */
  int nrhs;             /* nbr of rhs vars*/
  int nlhs;             /* nbr of lhs vars*/
  int is_fn;            /* curr symbol a fn?*/
  int var_ix;           /* curr symbol ix*/
  int ncovar;           /* number of covariates */
} Symtable;
Symtable symtab;

typedef struct Sbuf {
  char s[MXBUF];        /* curr print buffer */
  int o;                /* offset of print buffer */
} Sbuf;
Sbuf sb;                /* buffer w/ current parsed & translated line */
                        /* to be stored in a temp file */

static FILE *fpIO, *fp_inits;
int symstr_offset=0;


/* new symbol? if no, find it's ith */
int new_or_ith(const char *s) {
  int i, len, len_s=strlen(s);

  // not new symbols
  if (symtab.is_fn) return 0;
  if (!strcmp("t", s)) return 0;
  if (!strcmp("podo", s)) return 0;
  if (!strcmp("tlast", s)) return 0;

  // ignore?
  for (i=0; i<nignore; i++) {
    len = ignore_offset[i+1] - ignore_offset[i] - 1;    /* -1 for added ',' */
    //printf("%s\n", ignore_symb_str+ignore_offset[i]);
    if (!strncmp(ignore_symb_str+ignore_offset[i], s, max(len, len_s))) {/* note we need take the max in order not to match a sub-string */
      return 0;
    }
  }

  // starting point
  if (!symtab.nvar) return 1;

  for (i=0; i<symtab.nvar; i++) {
    len = symtab.offset[i+1] - symtab.offset[i] - 1;    /* -1 for added ',' */
    if (!strncmp(symtab.symb_str+symtab.offset[i], s, max(len, len_s))) {/* note we need take the max in order not to match a sub-string */
      symtab.var_ix = i;
      return 0;
    }
  }
  return 1;
}

void wprint_node(int depth, char *name, char *value, void *client_data) {
  sprintf(SBPTR, " %s", value);
  sb.o += strlen(value)+1;
}

void wprint_parsetree(D_ParserTables pt, D_ParseNode *pn, int depth, print_node_fn_t fn, void *client_data) {
  char *name = (char*)pt.symbols[pn->symbol].name;
  int nch = d_get_number_of_children(pn), i;
  char *value = (char*)rc_dup_str(pn->start_loc.s, pn->end);

  if (!strcmp("identifier", name) && new_or_ith(value)) {
    sprintf(symtab.symb_str+symstr_offset, "%s,", value);
    symstr_offset += strlen(value)+1;
    symtab.offset[++symtab.nvar] = symstr_offset;
  }

  if (!strcmp("(", name)) {sprintf(SBPTR, "("); sb.o++;}
  if (!strcmp(")", name)) {sprintf(SBPTR, ")"); sb.o++;}
  if (!strcmp(",", name)) {sprintf(SBPTR, ","); sb.o++;}

  if (
    !strcmp("identifier", name) ||
    !strcmp("constant", name) ||
    !strcmp("+", name) ||
    !strcmp("-", name) ||
    !strcmp("*", name) ||
    !strcmp("/", name) ||

    !strcmp("&&", name) ||
    !strcmp("||", name) ||
    !strcmp("!=", name) ||
    !strcmp("==", name) ||
    !strcmp("<=", name) ||
    !strcmp(">=", name) ||
    !strcmp("!", name) ||
    !strcmp("<", name) ||
    !strcmp(">", name) ||

    !strcmp("=", name)
      ){
    fn(depth, name, value, client_data);
  } else {
    /* Free(value); */
  }
  
  depth++;
  if (nch != 0) {

    if (!strcmp("power_expression", name)) {
      sprintf(SBPTR, " pow(");
      sb.o += 5;
    }

    for (i = 0; i < nch; i++) {
      if (!strcmp("derivative", name) && i< 2) continue;
      if (!strcmp("derivative", name) && i==3) continue;
      if (!strcmp("derivative", name) && i==4) continue;

      symtab.is_fn = (!strcmp("function", name) && i==0) ? 1 : 0;
      D_ParseNode *xpn = d_get_child(pn,i);
      wprint_parsetree(pt, xpn, depth, fn, client_data);

      //inits
      if (!strcmp("selection_statement", name) && i==1) {
        sprintf(sb.s, "if (");
        sb.o = strlen(sb.s);
        continue;
      }
      if (!strcmp("selection_statement", name) && i==3) {
        sprintf(SBPTR, " {");
        sb.o += 2;
        fprintf(fpIO, "%s\n", sb.s);
        continue;
      }
      if (!strcmp("selection_statement__8", name) && i==0) {
        fprintf(fpIO, "}\nelse {\n");
        continue;
      }

      if (!strcmp("power_expression", name) && i==0) {
        sprintf(SBPTR, ",");
        sb.o++;
      }

      if (!strcmp("derivative", name) && i==2) {
        sprintf(sb.s, "dxdt[%d] = InfusionRate[%d] +", symtab.nder, symtab.nder);
        sb.o = strlen(sb.s);

        char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
        new_or_ith(v);
        symtab.is_lhs[symtab.var_ix] = 9;
        symtab.der_ix[symtab.nder] = symtab.var_ix;
        symtab.nder++;
        Free(v);
        continue;
      }

      if (!strcmp("assignment", name) && i==0) {
        char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
        sprintf(sb.s, "%s", v);
        sb.o = strlen(v);

        new_or_ith(v);
        symtab.is_lhs[symtab.var_ix] = 1;
        Free(v);
      }
    } // end for

    if (!strcmp("assignment", name) || !strcmp("derivative", name))
      fprintf(fpIO, "%s;\n", sb.s);

    if (!strcmp("selection_statement", name))
      fprintf(fpIO, "}\n");

    if (!strcmp("power_expression", name)) {
      sprintf(SBPTR, ")");
      sb.o++;
    }
  } // end if
}

//retrieve the ith symbol
void retieve_var(int i, char *buf) {
  int len;

  len = symtab.offset[i+1] - symtab.offset[i] - 1;
  strncpy(buf, symtab.symb_str+symtab.offset[i], len);
  buf[len] = 0;
}

void err_msg(int chk, const char *msg, int code)
{
  if(!chk) {
    if (symtab.symb_str)
      free(symtab.symb_str);
    printf("%s", msg);
    //exit(code);
    return;
  }
}

//output aux_files -- to be read by dvode() in R.
void prnt_aux_files(char *prefix) {
  int i, islhs;
  char buf[512];
  FILE *fp[3];

  sprintf(buf, "%sODE_PARS.txt",   prefix); fp[0] = fopen(buf, "w");
  sprintf(buf, "%sLHS_VARS.txt",   prefix); fp[1] = fopen(buf, "w");
  sprintf(buf, "%sSTATE_VARS.txt", prefix); fp[2] = fopen(buf, "w");
  i =  (intptr_t)fp[0] * (intptr_t)fp[1] *  (intptr_t) fp[2];
  err_msg(i, "Coudln't open file to write.\n", -1);

  for (i=0; i<symtab.nvar; i++) {
    islhs = symtab.is_lhs[i];
    if (islhs>1) continue;    /* is a state var */
    retieve_var(i, buf);
    fprintf(fp[islhs], "%s ", buf);
  }

  for (i=0; i<symtab.nder; i++) {            /* name state vars */
    retieve_var(symtab.der_ix[i], buf);
    fprintf(fp[2], "%s ", buf);
  }

  fclose(fp[0]);
  fclose(fp[1]);
  fclose(fp[2]);
}

void codegen(FILE *outpt, int dosep) {
  int i, j, k, ncov;
  char sLine[MXLEN+1];
  char buf[64];
  FILE *fpIO;

#if 0
  char *hdft[]=
    {
      "#ifndef __STAN__MATH__FUNCTIONS__GENERIC_ODE_INTERFACE_HPP__\n#define __STAN__MATH__FUNCTIONS__GENERIC_ODE_INTERFACE_HPP__\n\n#include <vector>\n#include <stan/math/rev/core.hpp>\n#include <stan/math/prim/mat/fun/Eigen.hpp>\n#include <stan/math/PMXStan/lsoda_tmpl_class.hpp>#include <ctime>\n\n\nnamespace stan {\nnamespace math {\n\ntypedef stan::math::var AVAR;\ntypedef double ADBL;\n\nstatic const size_t neq = %d;\nstatic const size_t npar= %d;\nstatic const size_t ncov=%d;\nstatic std::vector<AVAR> pars_var(npar);\nstatic std::vector<ADBL> pars_dbl(npar, 0.0);\nstatic std::vector<ADBL> ptr_cov(ncov, 0.0);\nstatic std::vector<double> InfusionRate(neq, 0.0);\n\nvoid dydt(double t, AVAR *x, AVAR *dxdt, void *data) {\n",
      "void dydt(double t, ADBL *x, ADBL *dxdt, void *data) {\n",
      "}\n\n",
      "\n} // ns math\n}// ns stan\n#endif\n"
    };
#endif

  for (k=1; k<2; k++) {
    //fprintf(outpt, hdft[k], symtab.nder, symtab.nrhs, symtab.ncovar);

    fprintf(outpt, "  %s\n", "double");
    for (i=0, j=0; i<symtab.nvar; i++) {
      j++;
      retieve_var(i, buf);
      fprintf(outpt, j<symtab.nvar ? "  %s,\n":"  %s;\n\n", buf);
    }

    if (dosep) fprintf(outpt, "/*=================*/\n");

    for (i=0, j=0; i<symtab.nvar; i++) {
      if (symtab.is_lhs[i]) continue;
      retieve_var(i, buf);
      fprintf(outpt, "  %s = _phi(_i, %d);\n", buf, j++);
    }

    if (dosep) fprintf(outpt, "/*=================*/\n");

    for (i=0; i<symtab.nder; i++) {            /* name state vars */
      retieve_var(symtab.der_ix[i], buf);
      fprintf(outpt, "  %s = _x[%d];\n", buf, i);
    }
    fprintf(outpt,"\n");

    ncov = 0;
    fpIO = fopen("out2.txt", "r");
    err_msg((intptr_t)fpIO, "Coudln't access out2.txt.\n", -1);
    while(fgets(sLine, MXLEN, fpIO)) {        /* parsed eqns */
	  char *s;
	  s = strstr(sLine, "= 9999.999 + - 9999.999;");
	  if (s) {
        s[0] = '\0';
        fprintf(outpt, "  %s = _ptr_cov[%d];\n", sLine, ncov++);
	  }
	  else {
        fprintf(outpt, "  %s", sLine);
      }
    }
    fclose(fpIO);
    //fprintf(outpt, "%s", hdft[2]);
  }

  //fprintf(outpt, "%s", hdft[3]);
}


void codegen_stanode(FILE *outpt, const char* template_file) {
  int i, j, k, ncov;
  char sLine[MXLEN+1];
  char buf[64];
  FILE *fpIO;

  char *vartype[]={"AVAR", "ADBL"};
  char *partype[]={"var", "dbl"};

  char *hdft[]=
    {
      "#ifndef __STAN__MATH__FUNCTIONS__GENERIC_ODE_INTERFACE_HPP__\n#define __STAN__MATH__FUNCTIONS__GENERIC_ODE_INTERFACE_HPP__\n\n#include <vector>\n#include <stan/math/rev/core.hpp>\n#include <stan/math/prim/mat/fun/Eigen.hpp>\n#include <stan/math/PMXStan/lsoda_tmpl_class.hpp>\n#include <ctime>\n\n\nnamespace stan {\nnamespace math {\n\ntypedef stan::math::var AVAR;\ntypedef double ADBL;\n\nstatic const size_t neq = %d;\nstatic const size_t npar= %d;\nstatic const size_t ncov=%d;\nstatic std::vector<AVAR> pars_var(npar);\nstatic std::vector<ADBL> pars_dbl(npar, 0.0);\nstatic std::vector<ADBL> ptr_cov(ncov, 0.0);\nstatic std::vector<double> InfusionRate(neq, 0.0);\n\nvoid dydt(double t, AVAR *x, AVAR *dxdt, void *data) {\n",
      "void dydt(double t, ADBL *x, ADBL *dxdt, void *data) {\n",
      "}\n\n",
      "\n} // ns math\n}// ns stan\n#endif\n"
    };

  for (k=0; k<2; k++) {
    fprintf(outpt, hdft[k], symtab.nder, symtab.nrhs, symtab.ncovar);

    fprintf(outpt, "  %s\n", vartype[k]);
    for (i=0, j=0; i<symtab.nvar; i++) {
      j++;
      retieve_var(i, buf);
      fprintf(outpt, j<symtab.nvar ? "  %s,\n":"  %s;\n\n", buf);
    }

    for (i=0, j=0; i<symtab.nvar; i++) {
      if (symtab.is_lhs[i]) continue;
      retieve_var(i, buf);
      fprintf(outpt, "  %s = pars_%s[%d];\n", buf, partype[k], j++);
    }

    for (i=0; i<symtab.nder; i++) {            /* name state vars */
      retieve_var(symtab.der_ix[i], buf);
      fprintf(outpt, "  %s = x[%d];\n", buf, i);
    }
    fprintf(outpt,"\n");

    ncov = 0;
    fpIO = fopen("out2.txt", "r");
    err_msg((intptr_t)
	    fpIO, "Coudln't access out2.txt.\n", -1);
    while(fgets(sLine, MXLEN, fpIO)) {        /* parsed eqns */
	  char *s;
	  s = strstr(sLine, "= 9999.999 + - 9999.999;");
	  if (s) {
        s[0] = '\0';
        fprintf(outpt, "  %s = ptr_cov[%d];\n", sLine, ncov++);
	  }
	  else {
        fprintf(outpt, "  %s", sLine);
      }
    }
    fclose(fpIO);
    fprintf(outpt, "%s", hdft[2]);
  }

  fpIO = fopen(template_file, "r");
  err_msg((intptr_t)fpIO, "Couldn't access template file.\n", -1);
  while(fgets(sLine, MXLEN, fpIO)) {
    fprintf(outpt, "%s", sLine);
  }
  fclose(fpIO);
  fprintf(outpt, "%s", hdft[3]);
}



void inits() {
  symtab.symb_str = (char *) malloc(64*MXSYM*sizeof(int));
  err_msg((intptr_t)symtab.symb_str, "error allocating vars", 1);

  symtab.offset[0]=0;
  memset(symtab.is_lhs, 0, MXSYM*sizeof(int));
  symtab.nvar=0;
  symtab.nder=0;
  symtab.is_fn=0;
  symstr_offset=0;
}

#ifdef __STANDALONE_PARS__
int main(int argc, char *argv[]) {
#else
void parse_pars(char **model_file, char **result_file, int *nrhs, int *dosep, char **ign_str, int *offset, int *nign) {
  int argc = 3;
  char *argv[] = {"", *model_file, *result_file, *ign_str};
#endif

  int i;
  char *buf;

  ignore_symb_str = argv[3];
  ignore_offset = offset;
  nignore = *nign;
  //printf("%d %s\n", strlen(argv[3]), argv[3]);

  D_ParseNode *pn;
  /* any number greater than sizeof(D_ParseNode_User) will do;
     below 1024 is used */
  D_Parser *p = new_D_Parser(&parser_tables_gram, 1024);
  p->save_parse_tree = 1;

  if (argc<3) {
    printf("Usage: %s FILE_to_parse c_FILE [aux file direcory]\n",argv[0]);
    return;
  }
  else {
    buf = r_sbuf_read(argv[1]);
    err_msg((intptr_t)buf, "error: empty buf\n", -2);
  }

  if ((pn=dparse(p, buf, strlen(buf))) && !p->syntax_errors) {
    inits();
    fpIO = fopen( "out2.txt", "w" );
    

    err_msg((intptr_t)fpIO, "error opening out2.txt\n", -2);
    wprint_parsetree(parser_tables_gram, pn, 0, wprint_node, NULL);
    fclose(fpIO);
    if (fp_inits)
      fclose(fp_inits);

    // count number of rhs vars
    for (i=0, symtab.nrhs=0; i<symtab.nvar; i++) {
      if (symtab.is_lhs[i]==0) {
        symtab.nrhs ++;
	  }
	}
    symtab.nlhs = symtab.nvar - symtab.nrhs;
    //printf("symtab.nlhs=%d\n", symtab.nlhs);
    symtab.ncovar = 0;
    *nrhs = symtab.nrhs;

    fpIO = fopen(argv[2], "w");
    codegen(fpIO, *dosep);
    fclose(fpIO);
    prnt_aux_files("");
    remove("out2.txt");
    free(symtab.symb_str);
  }
  else {
    printf("\nfailure\n");
  }
  return;
}


#ifdef __STANDALONE_ODE__
int main(int argc, char *argv[]) {
#else
void parse_ode(char **tmplt_file, char **model_file, char **result_file, char **ncovar) {
  int argc = 5;
  char *argv[] = {"", *tmplt_file, *model_file, *result_file, *ncovar, ""};
#endif

  int i;
  char *buf;

  D_ParseNode *pn;
  /* any number greater than sizeof(D_ParseNode_User) will do;
     below 1024 is used */
  D_Parser *p = new_D_Parser(&parser_tables_gram, 1024);
  p->save_parse_tree = 1;

  if (argc<5) {
    printf("Usage: %s TEMPLATE_FILE FILE_to_parse c_FILE NCOVAR [aux file direcory]\n",argv[0]);
    return;
  }
  else {
    buf = r_sbuf_read(argv[2]);
    err_msg((intptr_t)buf, "error: empty buf\n", -2);
  }

  if ((pn=dparse(p, buf, strlen(buf))) && !p->syntax_errors) {
    inits();
    fpIO = fopen( "out2.txt", "w" );
    err_msg((intptr_t)fpIO, "error opening out2.txt\n", -2);
    wprint_parsetree(parser_tables_gram, pn, 0, wprint_node, NULL);
    fclose(fpIO);
    if (fp_inits)
      fclose(fp_inits);

    // count number of rhs vars
    for (i=0, symtab.nrhs=0; i<symtab.nvar; i++) {
      if (symtab.is_lhs[i]==0) {
        symtab.nrhs ++;
	  }
	}
    symtab.nlhs = symtab.nvar - symtab.nrhs;
    //printf("symtab.nlhs=%d\n", symtab.nlhs);
    symtab.ncovar = atoi(argv[4]);

    fpIO = fopen(argv[3], "w");
    codegen_stanode(fpIO, argv[1]);
    fclose(fpIO);
    prnt_aux_files(argc<6 ? "" : argv[5]);
    remove("out2.txt");
    free(symtab.symb_str);
  }
  else {
    printf("\nfailure\n");
  }
  return;
}

