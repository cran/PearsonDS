#include <float.h>
#include "config.h"
#include "xpfpa.h"

// #ifndef _FPU_GETCW
// #define _FPU_GETCW(x) asm volatile ("fnstcw %0":"=m" (x));
// #endif

// #ifndef _FPU_SETCW
// #define _FPU_SETCW(x) asm volatile ("fldcw %0": :"m" (x));
// #endif

// #ifndef _FPU_EXTENDED
// #define _FPU_EXTENDED 0x0300
// #endif

// #ifndef _FPU_DOUBLE
// #define _FPU_DOUBLE 0x0200
// #endif

typedef struct {
	double r;
	double e;
} flerr;

typedef struct {
	double r;
	double e;
} ddouble;

typedef struct {
	double s;
	double u;
	double v;
} daccres;

typedef struct {
	double hi;
	double lo;
} splitr;

typedef struct {
	double v[4];
} qdouble;

typedef struct {
	double v[5];
} unqdouble;

typedef struct {
	double v[3];
} tsares;

typedef struct {
	double v[2];
} tsbres;

typedef struct {
	qdouble r;
	qdouble i;
} QDcomplex;

typedef struct {
	ddouble r;
	ddouble i;
} DDcomplex;

// QuadDouble.c
flerr qtwosum(double a, double b);
flerr twosum(double a, double b);
splitr split(double a);
flerr twoprod(double a, double b);
qdouble renormalize(unqdouble a);
double qd2d (qdouble a);
qdouble qdplusd(qdouble a, double b);
daccres doubleacc(double u, double v, double x);
qdouble qdadd(qdouble a, qdouble b);
tsares threesumsa(double x, double y, double z);
tsbres threesumsb(double x, double y, double z);
qdouble qdadd2 (qdouble a, qdouble b);
qdouble qdtimesd (qdouble a, double b);
tsares sixthreesum(double a, double b, double c, double d, double e, double f);
flerr doubledoublesum(flerr a, flerr b);
flerr ninetwosum_old(double a, double b, double c, double d, double e, double f,
                     double g, double h, double i);
flerr ninetwosum(double a, double b, double c, double d, double e, double f,
                 double g, double h, double i);
qdouble qdmult (qdouble a, qdouble b);
qdouble qddiv (qdouble a, qdouble b);
qdouble qdneg (qdouble a);
QDcomplex QDMultC(QDcomplex A, Rcomplex B);
QDcomplex QDMult(QDcomplex A, QDcomplex B);
QDcomplex QDMultR(QDcomplex A, double B);
QDcomplex QDDivC(QDcomplex A, Rcomplex B);
QDcomplex QDDivR(QDcomplex A, double B);
QDcomplex QDAdd(QDcomplex A, QDcomplex B);
QDcomplex QDAddR(QDcomplex A, double B);
QDcomplex QDAdd1(QDcomplex A);
double QDnormInf(QDcomplex A);
void FPUcheck();

// DoubleDouble.c
double dd2d (ddouble a);
DDcomplex DDMultC(DDcomplex A, Rcomplex B);
DDcomplex DDMult(DDcomplex A, DDcomplex B);
DDcomplex DDMultR(DDcomplex A, double B);
DDcomplex DDDivC(DDcomplex A, Rcomplex B);
DDcomplex DDDivR(DDcomplex A, double B);
DDcomplex DDAdd(DDcomplex A, DDcomplex B);
DDcomplex DDAddR(DDcomplex A, double B);
DDcomplex DDAdd1(DDcomplex A);
double DDnormInf(DDcomplex A);

// Complex.c
Rcomplex CMult(Rcomplex A, Rcomplex B);
Rcomplex CMultR(Rcomplex A, double B);
Rcomplex CDiv(Rcomplex A, Rcomplex B);
Rcomplex CDivR(Rcomplex A, double B);
Rcomplex CAdd(Rcomplex A, Rcomplex B);
Rcomplex CAddR(Rcomplex A, double B);
Rcomplex CAdd1(Rcomplex A);
double Cabs2(Rcomplex A);
double CnormInf(Rcomplex A);
