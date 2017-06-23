#include <math.h>
#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))
extern long dadt_counter;
extern double InfusionRate[99];
extern double *par_ptr;
extern double podo;
extern double tlast;

// prj-specific differential eqns
void RxODE_mod_file4b7868a723c2_dydt(unsigned int neq, double t, double *__zzStateVar__, double *__DDtStateVar__)
{
double
	centr,
	K21,
	periph,
	K12,
	VM,
	V,
	KM;

	K21 = par_ptr[0];
	K12 = par_ptr[1];
	VM = par_ptr[2];
	V = par_ptr[3];
	KM = par_ptr[4];

	centr = __zzStateVar__[0];
	periph = __zzStateVar__[1];

	__DDtStateVar__[0] = InfusionRate[0] + K21 * periph - K12 * centr -( VM * centr / V) /( KM + centr / V);
	__DDtStateVar__[1] = InfusionRate[1] + - K21 * periph + K12 * centr;
    dadt_counter++;
}

// prj-specific derived vars
void RxODE_mod_file4b7868a723c2_calc_lhs(double t, double *__zzStateVar__, double *lhs) {
double
	centr,
	K21,
	periph,
	K12,
	VM,
	V,
	KM;

	K21 = par_ptr[0];
	K12 = par_ptr[1];
	VM = par_ptr[2];
	V = par_ptr[3];
	KM = par_ptr[4];

	centr = __zzStateVar__[0];
	periph = __zzStateVar__[1];


}
