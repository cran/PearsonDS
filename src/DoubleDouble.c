/* Code (partially) based on
     'Algorithms for Quad-Double Precision Floating Point Arithmetic'
   by
     Yozo Hida, Xiaoye S. Li, and David H. Bailey
*/      

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>
#include <Rmath.h>
#include "pearson.h"

ddouble qtwosumd(double a, double b) {
  volatile ddouble res;
  res.r = a + b;
  res.e = b - (res.r - a);
//  Rprintf("  Quick-Two-Sum; a=%g, b=%g, sum=%g, error=%g\n",a,b,res.r,res.e);
  return(res);
}

ddouble twosumd(double a, double b) {
  volatile double v;
  volatile ddouble res;
  res.r = a + b;
  v     = res.r - a;
  res.e = (a - (res.r - v)) + (b - v);
//  Rprintf("        Two-Sum; a=%g, b=%g, sum=%g, error=%g\n",a,b,res.r,res.e);
  return(res);
}

ddouble twoprodd(double a, double b) {
  volatile ddouble res;
  res.r = a * b;
  volatile splitr as, bs;
  as    = split(a);
  bs    = split(b);
  res.e = ((as.hi * bs.hi - res.r)+as.hi*bs.lo+as.lo*bs.hi)+as.lo*bs.lo;
  return(res);
}

double dd2d (ddouble a) {
  return(a.r+a.e);
}

ddouble ddplusd(ddouble a, double b) {
  ddouble tmp;
  tmp = twosumd(a.r,b);
  tmp = qtwosumd(tmp.r,tmp.e+a.e);
  return(tmp);
}

ddouble ddadd(ddouble a, ddouble b) {
  ddouble tmp1,tmp2;
  tmp1 = twosumd(a.r,b.r);
  tmp2 = twosumd(a.e,b.e);
  tmp1 = qtwosumd(tmp1.r,tmp1.e+tmp2.r);
  tmp1 = qtwosumd(tmp1.r,tmp1.e+tmp2.e);
  return(tmp1);
}

ddouble ddtimesd (ddouble a, double b) {
  volatile ddouble tmp1;
  tmp1 = twoprodd(a.r,b);
  tmp1 = qtwosumd(tmp1.r,tmp1.e + a.e*b);
  return(tmp1);
}

ddouble ddmult (ddouble a, ddouble b) {
  volatile ddouble tmp1;
  tmp1 = twoprodd(a.r,b.r);
  tmp1 = qtwosumd(tmp1.r,tmp1.e+(a.r*b.e + a.e*b.r));
  return(tmp1);
}

ddouble dddivd (ddouble a, double b) {
  volatile double q0, q1, err;
  volatile ddouble tmp1,tmp2,res;
  
  q0   = a.r / b;
  tmp1 = twoprodd(q0,b);
  tmp2 = twosumd(a.r,-tmp1.r);
  err  = a.e;
  err -= tmp1.e;
  q1   = (tmp2.r + err) / b;
  res  = qtwosumd(q0,q1);
  return(res);
}  

ddouble dddiv (ddouble a, ddouble b) {
  volatile double q0, q1, q2;
  volatile ddouble r;
  
  q0  = a.r / b.r;
  r   = ddadd(a,ddtimesd(b,-q0));
  q1  = r.r / b.r;
  r   = ddadd(r,ddtimesd(b,-q1));
  q2  = r.r / b.r;
  r   = qtwosumd(q0,q1);
  r   = ddplusd(r,q2);
  return(r);
}

ddouble ddneg (ddouble a) {
  ddouble res;
  res.r = -a.r;
  res.e = -a.e;
  return(res);
}

SEXP DDplusD (SEXP DD, SEXP D) {
  SEXP Res;
  PROTECT (Res=allocVector(REALSXP,2));
  double d = REAL(D)[0];
  double *res = REAL(Res);
  ddouble dd,tmp;
  dd.r = REAL(DD)[0]; dd.e = REAL(DD)[1];
  tmp = ddplusd(dd,d);
  res[0] = tmp.r; res[1] = tmp.e;
  UNPROTECT(1);
  return(Res);
}

SEXP DDmalD (SEXP DD, SEXP D) {
  SEXP Res;
  PROTECT (Res=allocVector(REALSXP,2));
  double d = REAL(D)[0];
  double *res = REAL(Res);
  ddouble dd,tmp;
  dd.r = REAL(DD)[0]; dd.e = REAL(DD)[1];
  tmp = ddtimesd(dd,d);
  res[0] = tmp.r; res[1] = tmp.e;
  UNPROTECT(1);
  return(Res);
}

SEXP DDplusDD (SEXP DD1, SEXP DD2) {
  SEXP Res;
  PROTECT (Res=allocVector(REALSXP,2));
  double *res = REAL(Res);
  ddouble tmp,dd1,dd2;
  dd1.r = REAL(DD1)[0]; dd1.e = REAL(DD1)[1];
  dd2.r = REAL(DD2)[0]; dd2.e = REAL(DD2)[1];
  tmp = ddadd(dd1,dd2);
  res[0] = tmp.r; res[1] = tmp.e;
  UNPROTECT(1);
  return(Res);
}

SEXP DDmalDD (SEXP DD1, SEXP DD2) {
  SEXP Res;
  PROTECT (Res=allocVector(REALSXP,2));
  double *res = REAL(Res);
  ddouble tmp,dd1,dd2;
  dd1.r = REAL(DD1)[0]; dd1.e = REAL(DD1)[1];
  dd2.r = REAL(DD2)[0]; dd2.e = REAL(DD2)[1];
  tmp = ddmult(dd1,dd2);
  res[0] = tmp.r; res[1] = tmp.e;
  UNPROTECT(1);
  return(Res);
}

SEXP DDdurchDD (SEXP DD1, SEXP DD2) {
  SEXP Res;
  PROTECT (Res=allocVector(REALSXP,2));
  double *res = REAL(Res);
  ddouble tmp,dd1,dd2;
  dd1.r = REAL(DD1)[0]; dd1.e = REAL(DD1)[1];
  dd2.r = REAL(DD2)[0]; dd2.e = REAL(DD2)[1];
  tmp = dddiv(dd1,dd2);
  res[0] = tmp.r; res[1] = tmp.e;
  UNPROTECT(1);
  return(Res);
}

DDcomplex DDMultC(DDcomplex A, Rcomplex B) {
  DDcomplex C;
  C.r = ddadd(ddtimesd(A.r,B.r),ddtimesd(A.i,-B.i));
  C.i = ddadd(ddtimesd(A.r,B.i),ddtimesd(A.i,B.r));
  return(C);
}

DDcomplex DDMult(DDcomplex A, DDcomplex B) {
  DDcomplex C;
  C.r = ddadd(ddmult(A.r,B.r),ddmult(A.i,ddneg(B.i)));
  C.i = ddadd(ddmult(A.r,B.i),ddmult(A.i,B.r));
  return(C);
}

DDcomplex DDMultR(DDcomplex A, double B) {
  DDcomplex C;
  C.r = ddtimesd(A.r,B);
  C.i = ddtimesd(A.i,B);
  return(C);
}

DDcomplex DDDivC(DDcomplex A, Rcomplex B) {                                     // need Smith's formula here? TODO!
  DDcomplex C;
  ddouble nenner;
  ddouble zaehler;
  ddouble c;
  ddouble d;
  c.r  = B.r; c.e = 0.;
  d.r  = B.i; d.e = 0.;
  nenner  = ddadd(ddtimesd(c,B.r),ddtimesd(d,B.i));
  zaehler = ddadd(ddtimesd(A.r,B.r),ddtimesd(A.i,B.i));
  C.r     = dddiv(zaehler,nenner);
  zaehler = ddadd(ddtimesd(A.i,B.r),ddtimesd(A.r,-B.i));
  C.i     = dddiv(zaehler,nenner);
  return(C);
}

DDcomplex DDDivR(DDcomplex A, double B) {                                       // improvable!!!
  DDcomplex C;
  ddouble b;
  b.r = B; b.e = 0.;
  C.r     = dddiv(A.r,b);
  C.i     = dddiv(A.i,b);
  return(C);
}

DDcomplex DDAdd(DDcomplex A, DDcomplex B) {
  DDcomplex C;
  C.r = ddadd(A.r,B.r);
  C.i = ddadd(A.i,B.i);
  return(C);
}

DDcomplex DDAddR(DDcomplex A, double B) {
  DDcomplex C;
  C.r = ddplusd(A.r,B);
  C.i = A.i;
  return(C);
}

DDcomplex DDAdd1(DDcomplex A) {
  DDcomplex C;
  C.r = ddplusd(A.r , 1.);
  C.i = A.i;
  return(C);
}

double DDnormInf(DDcomplex A) {
  return(fmax2(fabs(dd2d(A.r)),fabs(dd2d(A.i))));
}
