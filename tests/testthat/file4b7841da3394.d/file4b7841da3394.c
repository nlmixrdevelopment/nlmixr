#include <math.h>
#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))
extern long dadt_counter;
extern double InfusionRate[99];
extern double *par_ptr;
extern double podo;
extern double tlast;

// prj-specific differential eqns
void RxODE_mod_file4b7841da3394_dydt(unsigned int neq, double t, double *__zzStateVar__, double *__DDtStateVar__)
{
double
	abs,
	KA,
	centr,
	K21,
	periph,
	K12,
	VM,
	V,
	KM;

	KA = par_ptr[0];
	K21 = par_ptr[1];
	K12 = par_ptr[2];
	VM = par_ptr[3];
	V = par_ptr[4];
	KM = par_ptr[5];

	abs = __zzStateVar__[0];
	centr = __zzStateVar__[1];
	periph = __zzStateVar__[2];

	__DDtStateVar__[0] = InfusionRate[0] + - KA * abs;
	__DDtStateVar__[1] = InfusionRate[1] + KA * abs + K21 * periph - K12 * centr -( VM * centr / V) /( KM + centr / V);
	__DDtStateVar__[2] = InfusionRate[2] + - K21 * periph + K12 * centr;
    dadt_counter++;
}

// prj-specific derived vars
void RxODE_mod_file4b7841da3394_calc_lhs(double t, double *__zzStateVar__, double *lhs) {
double
	abs,
	KA,
	centr,
	K21,
	periph,
	K12,
	VM,
	V,
	KM;

	KA = par_ptr[0];
	K21 = par_ptr[1];
	K12 = par_ptr[2];
	VM = par_ptr[3];
	V = par_ptr[4];
	KM = par_ptr[5];

	abs = __zzStateVar__[0];
	centr = __zzStateVar__[1];
	periph = __zzStateVar__[2];


}
