#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void cdgqf ( int nt, int kind, double alpha, double beta, double t[], 
  double wts[] );
void cgqf ( int nt, int kind, double alpha, double beta, double a, double b, 
  int lo, double t[], double wts[] );
void chebyshev_set ( int n, double x[], double w[] );
void chebyshev1_compute ( int n, double x[], double w[] );
double chebyshev1_integral ( int expon );
void chebyshev1_set ( int n, double x[], double w[] );
void chebyshev2_compute ( int n, double x[], double w[] );
double chebyshev2_integral ( int expon );
void chebyshev2_set ( int n, double x[], double w[] );
void chebyshev3_compute ( int n, double xtab[], double weight[] );
double chebyshev3_integral ( int expon );
void chebyshev3_set ( int n, double x[], double w[] );
double class_matrix ( int kind, int m, double alpha, double beta, double aj[], 
  double bj[] );
void clenshaw_curtis_compute ( int n, double x[], double w[] );
void clenshaw_curtis_set ( int n, double xtab[], double weight[] );
void fejer1_compute ( int n, double x[], double w[] );
void fejer1_set ( int n, double xtab[], double weight[] );
void fejer2_compute ( int n, double x[], double w[] );
void fejer2_set ( int n, double xtab[], double weight[] );
double gegenbauer_integral ( int expon, double alpha );
void gegenbauer_ek_compute ( int n, double alpha, double a, double b,
  double xtab[], double weight[] );
void gegenbauer_ss_compute ( int n, double alpha, double xtab[], 
  double weight[] );
void gegenbauer_ss_recur ( double *p2, double *dp2, double *p1, double x, 
  int order, double alpha, double c[] );
void gegenbauer_ss_root ( double *x, int order, double alpha,  double *dp2, 
  double *p1, double c[] );
void gen_hermite_dr_compute ( int order, double alpha, double x[], double w[] );
void gen_hermite_ek_compute ( int order, double alpha, double x[], double w[] );
double gen_hermite_integral ( int expon, double alpha );
void gen_laguerre_ek_compute ( int order, double alpha, double x[], 
  double w[] );
double gen_laguerre_integral ( int expon, double alpha );
void gen_laguerre_ss_compute ( int order, double alpha, double x[], 
  double w[] );
void gen_laguerre_ss_recur ( double *p2, double *dp2, double *p1, double x, 
  int order, double alpha, double b[], double c[] );
void gen_laguerre_ss_root ( double *x, int order, double alpha, double *dp2, 
  double *p1, double b[], double c[] );
void hermite_ek_compute ( int n, double x[], double w[] );
void hermite_gk16_set ( int n, double x[], double w[] );
void hermite_gk18_set ( int n, double x[], double w[] );
void hermite_gk22_set ( int n, double x[], double w[] );
void hermite_gk24_set ( int n, double x[], double w[] );
double hermite_integral ( int n );
double hermite_integral2 ( double a );
void hermite_probabilist_set ( int order, double xtab[], double weight[] );
void hermite_set ( int order, double xtab[], double weight[] );
void hermite_1_set ( int order, double xtab[], double weight[] );
void hermite_ss_compute ( int order, double xtab[], double weight[] );
void hermite_ss_recur ( double *p2, double *dp2, double *p1, double x, 
  int order );
void hermite_ss_root ( double *x, int order, double *dp2, double *p1 );
int i4_factorial2 ( int n );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_power ( int i, int j );
void imtqlx ( int n, double d[], double e[], double z[] );
void jacobi_ek_compute ( int n, double alpha, double beta, double x[], 
  double w[] );
double jacobi_integral ( int expon, double alpha, double beta );
void jacobi_ss_compute ( int n, double alpha, double beta, double x[], 
  double w[] );
void jacobi_ss_root ( double *x, int order, double alpha, double beta, 
  double *dp2, double *p1, double b[], double c[] );
void jacobi_ss_recur ( double *p2, double *dp2, double *p1, double x, int order, 
  double alpha, double beta, double b[], double c[] );
void kronrod_set ( int order, double xtab[], double weight[] );
void laguerre_ek_compute ( int n, double x[], double w[] );
double laguerre_integral ( int expon );
void laguerre_ss_compute ( int n, double x[], double w[] );
void laguerre_ss_recur ( double *p2, double *dp2, double *p1, double x, 
  int order, double b[], double c[] );
void laguerre_ss_root ( double *x, int order, double *dp2, 
  double *p1, double b[], double c[] );
void laguerre_set ( int order, double xtab[], double weight[] );
void laguerre_1_set ( int order, double xtab[], double weight[] );
double laguerre_sum ( double func ( double x ), double a, int order, 
  double xtab[], double weight[] );
void legendre_dr_compute ( int order, double xtab[], double weight[] );
void legendre_ek_compute ( int n, double x[], double w[] );
double legendre_integral ( int expon );
void legendre_recur ( double *p2, double *dp2, double *p1, double x, 
  int order );
void legendre_set ( int order, double xtab[], double weight[] );
void lobatto_compute ( int n, double x[], double w[] );
void lobatto_set ( int order, double xtab[], double weight[] );
double log_gamma ( double x );
void nc_compute_weights ( int order, double a, double b, double xtab[], 
  double weight[] );
void ncc_compute ( int order, double xtab[], double weight[] );
void ncc_set ( int order, double xtab[], double weight[] );
void nco_compute ( int n, double x[], double w[] );
void nco_set ( int order, double xtab[], double weight[] );
void ncoh_compute ( int order, double xtab[], double weight[] );
void ncoh_set ( int order, double xtab[], double weight[] );
void parchk ( int kind, int m, double alpha, double beta );
void patterson_set ( int order, double xtab[], double weight[] );
void psi_values ( int *n_data, double *x, double *fx );
double r8_abs ( double x );
double r8_epsilon ( );
double r8_factorial ( int n );
double r8_factorial2 ( int n );
double r8_gamma ( double x );
double r8_huge ( );
double r8_hyper_2f1 ( double a, double b, double c, double x );
double r8_max ( double x, double y );
double r8_psi ( double xx );
double r8_sign ( double x );
void r8vec_copy ( int n, double a1[], double a2[] );
double r8vec_dot_product ( int n, double a1[], double a2[] );
void r8vec_linspace ( int n, double a, double b, double x[] );
void r8vec_print ( int n, double a[], char *title );
void r8vec_reverse ( int n, double x[] );
void radau_compute ( int n, double x[], double w[] );
void radau_set ( int order, double xtab[], double weight[] );
void scqf ( int nt, double t[], int mlt[], double wts[], int nwts, int ndx[], 
  double swts[], double st[], int kind, double alpha, double beta, double a, 
  double b );
void sgqf ( int nt, double aj[], double bj[], double zemu, double t[], 
  double wts[] );
void timestamp ( );
