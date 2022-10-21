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

flerr qtwosum(double a, double b) {
  volatile flerr res;
  res.r = a + b;
  res.e = b - (res.r - a);
//  Rprintf("  Quick-Two-Sum; a=%g, b=%g, sum=%g, error=%g\n",a,b,res.r,res.e);
  return(res);
}

flerr twosum(double a, double b) {
  volatile double v;
  volatile flerr res;
  res.r = a + b;
  v     = res.r - a;
  res.e = (a - (res.r - v)) + (b - v);
//  Rprintf("        Two-Sum; a=%g, b=%g, sum=%g, error=%g\n",a,b,res.r,res.e);
  return(res);
}

splitr split(double a) {
  volatile double t;
  volatile splitr res;
  t = 134217729. * a;
  res.hi = t - (t - a);
  res.lo = a - res.hi;
  return(res);
}

SEXP Split(SEXP A) {
  splitr tmp;
  SEXP Res;
  PROTECT (Res=allocVector(REALSXP,2));
  double *res = REAL(Res);
  tmp = split(REAL(A)[0]);
  res[0] = tmp.hi;
  res[1] = tmp.lo;
  UNPROTECT(1);
  return(Res);
}

flerr twoprod(double a, double b) {
  volatile flerr res;
  res.r = a * b;
  volatile splitr as, bs;
  as    = split(a);
  bs    = split(b);
  res.e = ((as.hi * bs.hi - res.r)+as.hi*bs.lo+as.lo*bs.hi)+as.lo*bs.lo;
  return(res);
}

qdouble renormalize(unqdouble a) {
  int k,i;
  volatile double t[5];
  volatile flerr tmp;
  qdouble b;
  for (i=0;i<4;i++) b.v[i] = 0.;
  double s,e;
  tmp  = qtwosum(a.v[3],a.v[4]);
  t[4] = tmp.e;
  tmp  = qtwosum(a.v[2],tmp.r);
  t[3] = tmp.e;
  tmp  = qtwosum(a.v[1],tmp.r);
  t[2] = tmp.e;
  tmp  = qtwosum(a.v[0],tmp.r);
  t[0] = tmp.r;
  t[1] = tmp.e;
  s    = t[0];
  k    = 0;
//  Rprintf("  Intermediate t: %g %g %g %g %g\n",t[0],t[1],t[2],t[3],t[4]);
  for (i=1;i<=4;i++) {
    tmp = qtwosum(s,t[i]);
    s   = tmp.r;
    e   = tmp.e;
    b.v[k] = s;
    if (e!=0.) {
      s      = e;
      k      = k+1;
    }
  }
//  Rprintf("  Renormalize; a=%g+%g+%g+%g+%g\n",a.v[0],a.v[1],a.v[2],a.v[3],a.v[4]);
//  Rprintf("       result: b=%g+%g+%g+%g\n",b.v[0],b.v[1],b.v[2],b.v[3]);
  return(b);
}

double qd2d (qdouble a) {
  return(a.v[0]+a.v[1]+a.v[2]+a.v[3]);
}

qdouble qdplusd(qdouble a, double b) {
  unqdouble us;
  flerr tmp;
  tmp = twosum(a.v[0],b);
  us.v[0] = tmp.r;
  tmp = twosum(a.v[1],tmp.e);
  us.v[1] = tmp.r;
  tmp = twosum(a.v[2],tmp.e);
  us.v[2] = tmp.r;
  tmp = twosum(a.v[3],tmp.e);
  us.v[3] = tmp.r;
  us.v[4] = tmp.e;
//  Rprintf(" sum of %g + %g + %g + %g and %g\n",a.v[0],a.v[1],a.v[2],a.v[3],b);
  return(renormalize(us));
}

daccres doubleacc(double u, double v, double x) {                               
//  Rprintf("DACC In: %20g, %20g, %20g\n",u,v,x);
  daccres res;
  flerr tmp;
  double vint,uint,s;
  tmp  = twosum(v,x);
  vint = tmp.e;
  tmp  = twosum(u,tmp.r);
  uint = tmp.e;
  s    = tmp.r;
  if (uint==0.) {
    uint = s; s = 0.0;
  }
  if (vint==0.) {
    vint = uint; uint = s; s = 0.0;
  }
  res.s = s; res.u = uint; res.v = vint;
//  Rprintf("DACC Out: %20g, %20g, %20g\n",res.s,res.u,res.v);  
  return(res);
}

qdouble qdadd(qdouble a, qdouble b) {
  double x[8];
  int i=0,j=0,k;
  daccres   tmp;
  unqdouble us;
  for (k=0;k<5;k++) us.v[k] = 0.;
  for (k=0;k<8;k++) {                                                            // merge - sort
    if (i>3) { x[k] = b.v[j]; j++; } else
    if (j>3) { x[k] = a.v[i]; i++; } else
    if (fabs(a.v[i])>=fabs(b.v[j])) { x[k] = a.v[i]; i++; } else { x[k] = b.v[j]; j++; }
  }
//  Rprintf("%g %g %g %g %g %g %g %g\n",x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7]);
  double s=0., u=0., v=0.;
  k = 0; i = 0;
  while ((k<4)&&(i<8)) {
    tmp = doubleacc(u,v,x[i]);
//    Rprintf("-   i=%d, k=%d, u=%f, v=%f, x[i]=%f\n",i,k,u,v,x[i]);
//    Rprintf("+   i=%d, k=%d, s=%f, u=%f, v[i]=%f\n",i,k,tmp.s,tmp.u,tmp.v);
    s = tmp.s; u = tmp.u; v = tmp.v;
    if (s!=0.) {
//      Rprintf("   i=%d, k=%d, v[k]=%f\n",i,k,s);
      us.v[k] = s; k++;
    }
    i++;
  }
  if (k<3) us.v[k+1] = v;
  if (k<4) us.v[k] = u;
//  Rprintf("%g %g %g %g %g\n",us.v[0],us.v[1],us.v[2],us.v[3],us.v[4]);
  return(renormalize(us));
}

tsares threesumsa(double x, double y, double z) {
  volatile tsares res;
  volatile flerr tmp1,tmp2;
  tmp1     = twosum(x,y);
  tmp2     = twosum(tmp1.r,z);
  res.v[0] = tmp2.r;
  tmp1     = twosum(tmp1.e,tmp2.e);
  res.v[1] = tmp1.r;
  res.v[2] = tmp1.e;
//  Rprintf("      Three-Sum; a=%g, b=%g, c=%g, sum=%g, error1=%g, error2=%g\n",x,y,z,res.v[0],res.v[1],res.v[2]);
  return(res);
}

tsbres threesumsb(double x, double y, double z) {
  volatile tsbres res;
  volatile flerr tmp1,tmp2;
  tmp1     = twosum(x,y);
  tmp2     = twosum(tmp1.r,z);
  res.v[0] = tmp2.r;
  res.v[1] = tmp1.e+tmp2.e;
//  Rprintf("      Three-Sum; a=%g, b=%g, c=%g, sum=%g, error1=%g\n",x,y,z,res.v[0],res.v[1],res.v[2]);
  return(res);
}

qdouble qdadd2 (qdouble a, qdouble b) {
  unqdouble us;
  flerr tmp1, tmp2, tmp3, tmp4;
  tsares tmpA;
  tsbres tmpB;
  tmp1 = twosum(a.v[0],b.v[0]);
  us.v[0] = tmp1.r;
  tmp2 = twosum(a.v[1],b.v[1]);
  tmp1 = twosum(tmp2.r,tmp1.e);
  us.v[1] = tmp1.r;
  tmp3 = twosum(a.v[2],b.v[2]);
  tmpA = threesumsa(tmp3.r,tmp1.e,tmp2.e);
  us.v[2] = tmpA.v[0];
  tmp4 = twosum(a.v[3],b.v[3]);
  tmpB = threesumsb(tmp4.r,tmp3.e,tmpA.v[1]);
  us.v[3] = tmpB.v[0];
  us.v[4] = tmpB.v[1]+tmpA.v[2]+tmp4.e;
  return(renormalize(us));
}

qdouble qdtimesd (qdouble a, double b) {
  unqdouble us;
  volatile flerr tmp1, tmp2, tmp3;
  tmp1 = twoprod(a.v[0],b);
  tmp2 = twoprod(a.v[1],b);
  tmp3 = twoprod(a.v[2],b);
  us.v[0] = tmp1.r;
  tmp1 = twosum(tmp2.r,tmp1.e);
  us.v[1] = tmp1.r;
  volatile tsares tmp4 = threesumsa(tmp3.r,tmp2.e,tmp1.e);
  us.v[2] = tmp4.v[0];
  volatile tsbres tmp5 = threesumsb(a.v[3]*b,tmp3.e,tmp4.v[1]);
  us.v[3] = tmp5.v[0];
  us.v[4] = tmp5.v[1]+tmp4.v[2];
  return(renormalize(us));
}

tsares sixthreesum(double a, double b, double c, double d, double e, double f) {
  volatile tsares res,tmp1,tmp2;
  volatile flerr tmp3, tmp4, tmp5;
  tmp1 = threesumsa(a,b,c);
  tmp2 = threesumsa(d,e,f);
  tmp3 = twosum(tmp1.v[0],tmp2.v[0]);
  res.v[0] = tmp3.r;
  tmp4 = twosum(tmp1.v[1],tmp2.v[1]);
  tmp5 = twosum(tmp4.r,tmp3.e);
  res.v[1] = tmp5.r;
  res.v[2] = tmp1.v[2]+tmp2.v[2]+tmp4.e+tmp5.e;
  return(res);
}

flerr doubledoublesum(flerr a, flerr b) {
  volatile flerr res,tmp;
  tmp   = twosum(a.r,b.r);
  res.r = tmp.r;
  res.e = tmp.e+a.e+b.e;
  return(res);
}

flerr ninetwosum_old(double a, double b, double c, double d, double e, double f,
                     double g, double h, double i) {
  volatile flerr res;
  volatile flerr tmp1,tmp2,tmp3,tmp4;
  tmp1  = twosum(a,b);
  tmp2  = twosum(c,d);
  tmp3  = twosum(e,f);
  tmp4  = twosum(g,h);
  tmp1  = doubledoublesum(tmp1,tmp2);
  tmp2  = doubledoublesum(tmp3,tmp4);
  tmp3  = doubledoublesum(tmp1,tmp2);
  tmp4  = twosum(tmp3.r,i);
  res.r = tmp4.r;
  res.e = tmp4.e+tmp3.e;
  return(res);
}

flerr ninetwosum(double a, double b, double c, double d, double e, double f,
                 double g, double h, double i) {
  volatile flerr res;
  volatile tsbres tmp1,tmp2,tmp3,tmp4;
  tmp1  = threesumsb(a,b,c);
  tmp2  = threesumsb(d,e,f);
  tmp3  = threesumsb(g,h,i);
  tmp4  = threesumsb(tmp1.v[0],tmp2.v[0],tmp3.v[0]);
  res.r = tmp4.v[0];
  res.e = tmp4.v[1]+tmp1.v[1]+tmp2.v[1]+tmp3.v[1];
  return(res);
}

qdouble qdmult (qdouble a, qdouble b) {
  unqdouble us;
  volatile flerr tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9;
  volatile tsares tmpA,tmpB;
  volatile flerr tmpC;
  tmp0    = twoprod(a.v[0],b.v[0]);
  us.v[0] = tmp0.r;
  tmp1    = twoprod(a.v[0],b.v[1]);
  tmp2    = twoprod(a.v[1],b.v[0]);
  tmpA    = threesumsa(tmp0.e,tmp1.r,tmp2.r);
  us.v[1] = tmpA.v[0];
  tmp3    = twoprod(a.v[0],b.v[2]);
  tmp4    = twoprod(a.v[1],b.v[1]);
  tmp5    = twoprod(a.v[2],b.v[0]);
  tmpB    = sixthreesum(tmpA.v[1],tmp1.e,tmp2.e,tmp3.r,tmp4.r,tmp5.r);
  us.v[2] = tmpB.v[0];
  tmp6    = twoprod(a.v[0],b.v[3]);
  tmp7    = twoprod(a.v[1],b.v[2]);
  tmp8    = twoprod(a.v[2],b.v[1]);
  tmp9    = twoprod(a.v[3],b.v[0]);
  tmpC    = ninetwosum(tmpA.v[2],tmpB.v[1],tmp3.e,tmp4.e,tmp5.e,
                       tmp6.r,tmp7.r,tmp8.r,tmp9.r);
  us.v[3] = tmpC.r;
  us.v[4] = tmpB.v[2] + tmpC.e + tmp6.e + tmp7.e + tmp8.e + tmp9.e +
            a.v[1]*b.v[3] + a.v[2]*b.v[2] + a.v[3]*b.v[1];
  return(renormalize(us));
}

qdouble qddiv (qdouble a, qdouble b) {
  volatile double q0, q1, q2, q3, q4;
  unqdouble us;
  volatile qdouble r;
  
  q0 = a.v[0] / b.v[0];
  r  = qdadd(a,qdtimesd(b,-q0));
  q1 = r.v[0] / b.v[0];
  r  = qdadd(r,qdtimesd(b,-q1));
  q2 = r.v[0] / b.v[0];
  r  = qdadd(r,qdtimesd(b,-q2));
  q3 = r.v[0] / b.v[0];
  r  = qdadd(r,qdtimesd(b,-q3));
  q4 = r.v[0] / b.v[0];
  us.v[0] = q0; us.v[1] = q1; us.v[2] = q2; us.v[3] = q3; us.v[4] = q4;
  return(renormalize(us));
}

qdouble qdneg (qdouble a) {
  qdouble res;
  for (int i=0;i<4;i++) res.v[i] = -a.v[i];
  return(res);
}

SEXP RENORM1 (SEXP QD) {
  SEXP Res;
  PROTECT (Res=allocVector(REALSXP,4));
  double *res = REAL(Res);
  qdouble tmp;
  unqdouble qd1;
  for (int i=0;i<4;i++) qd1.v[i] = REAL(QD)[i];
  qd1.v[4] = 0.0;
  tmp = renormalize(qd1);
  for (int i=0;i<4;i++) res[i] = tmp.v[i];
  UNPROTECT(1);
  return(Res);
}

SEXP QDplusD (SEXP QD, SEXP D) {
  SEXP Res;
  PROTECT (Res=allocVector(REALSXP,4));
  double d = REAL(D)[0];
  double *res = REAL(Res);
  qdouble q,tmp;
  for (int i=0;i<4;i++) q.v[i] = REAL(QD)[i];
  tmp = qdplusd(q,d);
  for (int i=0;i<4;i++) res[i] = tmp.v[i];
  UNPROTECT(1);
  return(Res);
}

SEXP QDmalD (SEXP QD, SEXP D) {
  SEXP Res;
  PROTECT (Res=allocVector(REALSXP,4));
  double d = REAL(D)[0];
  double *res = REAL(Res);
  qdouble q,tmp;
  for (int i=0;i<4;i++) q.v[i] = REAL(QD)[i];
  tmp = qdtimesd(q,d);
  for (int i=0;i<4;i++) res[i] = tmp.v[i];
  UNPROTECT(1);
  return(Res);
}

SEXP QDplusQD (SEXP QD1, SEXP QD2) {
  SEXP Res;
  PROTECT (Res=allocVector(REALSXP,4));
  double *res = REAL(Res);
  qdouble tmp,qd1,qd2;
  for (int i=0;i<4;i++) qd1.v[i] = REAL(QD1)[i];
  for (int i=0;i<4;i++) qd2.v[i] = REAL(QD2)[i];
  tmp = qdadd(qd1,qd2);
  for (int i=0;i<4;i++) res[i] = tmp.v[i];
  UNPROTECT(1);
  return(Res);
}

SEXP QDplus2QD (SEXP QD1, SEXP QD2) {
  SEXP Res;
  PROTECT (Res=allocVector(REALSXP,4));
  double *res = REAL(Res);
  qdouble tmp,qd1,qd2;
  for (int i=0;i<4;i++) qd1.v[i] = REAL(QD1)[i];
  for (int i=0;i<4;i++) qd2.v[i] = REAL(QD2)[i];
  tmp = qdadd2(qd1,qd2);
  for (int i=0;i<4;i++) res[i] = tmp.v[i];
  UNPROTECT(1);
  return(Res);
}

SEXP QDmalQD (SEXP QD1, SEXP QD2) {
  SEXP Res;
  PROTECT (Res=allocVector(REALSXP,4));
  double *res = REAL(Res);
  qdouble tmp,qd1,qd2;
  for (int i=0;i<4;i++) qd1.v[i] = REAL(QD1)[i];
  for (int i=0;i<4;i++) qd2.v[i] = REAL(QD2)[i];
  tmp = qdmult(qd1,qd2);
  for (int i=0;i<4;i++) res[i] = tmp.v[i];
  UNPROTECT(1);
  return(Res);
}

SEXP QDdurchQD (SEXP QD1, SEXP QD2) {
  SEXP Res;
  PROTECT (Res=allocVector(REALSXP,4));
  double *res = REAL(Res);
  qdouble tmp,qd1,qd2;
  for (int i=0;i<4;i++) qd1.v[i] = REAL(QD1)[i];
  for (int i=0;i<4;i++) qd2.v[i] = REAL(QD2)[i];
  tmp = qddiv(qd1,qd2);
  for (int i=0;i<4;i++) res[i] = tmp.v[i];
  UNPROTECT(1);
  return(Res);
}

QDcomplex QDMultC(QDcomplex A, Rcomplex B) {
  QDcomplex C;
  C.r = qdadd(qdtimesd(A.r,B.r),qdtimesd(A.i,-B.i));
  C.i = qdadd(qdtimesd(A.r,B.i),qdtimesd(A.i,B.r));
  return(C);
}

QDcomplex QDMult(QDcomplex A, QDcomplex B) {
  QDcomplex C;
  C.r = qdadd(qdmult(A.r,B.r),qdmult(A.i,qdneg(B.i)));
  C.i = qdadd(qdmult(A.r,B.i),qdmult(A.i,B.r));
  return(C);
}

QDcomplex QDMultR(QDcomplex A, double B) {
  QDcomplex C;
  C.r = qdtimesd(A.r,B);
  C.i = qdtimesd(A.i,B);
  return(C);
}

QDcomplex QDDivC(QDcomplex A, Rcomplex B) {                                     // need Smith's formula here? TODO!
  QDcomplex C;
  qdouble nenner;
  qdouble zaehler;
  qdouble c;
  qdouble d;
  c.v[0]  = B.r; c.v[1] = 0.; c.v[2] = 0.; c.v[3] = 0.;
  d.v[0]  = B.i; d.v[1] = 0.; d.v[2] = 0.; d.v[3] = 0.;
  nenner  = qdadd(qdtimesd(c,B.r),qdtimesd(d,B.i));
  zaehler = qdadd(qdtimesd(A.r,B.r),qdtimesd(A.i,B.i));
  C.r     = qddiv(zaehler,nenner);
  zaehler = qdadd(qdtimesd(A.i,B.r),qdtimesd(A.r,-B.i));
  C.i     = qddiv(zaehler,nenner);
  return(C);
}

QDcomplex QDDivR(QDcomplex A, double B) {                                       // improvable!!!
  QDcomplex C;
  qdouble b;
  b.v[0] = B; b.v[1] = 0.; b.v[2] = 0.; b.v[3] = 0.;
  C.r     = qddiv(A.r,b);
  C.i     = qddiv(A.i,b);
  return(C);
}

QDcomplex QDAdd(QDcomplex A, QDcomplex B) {
  QDcomplex C;
  C.r = qdadd(A.r,B.r);
  C.i = qdadd(A.i,B.i);
  return(C);
}

QDcomplex QDAddR(QDcomplex A, double B) {
  QDcomplex C;
  C.r = qdplusd(A.r,B);
  C.i = A.i;
  return(C);
}

QDcomplex QDAdd1(QDcomplex A) {
  QDcomplex C;
  C.r = qdplusd(A.r , 1.);
  C.i = A.i;
  return(C);
}

double QDnormInf(QDcomplex A) {
  return(fmax2(fabs(qd2d(A.r)),fabs(qd2d(A.i))));
}

double FPUcheck(void) {
  qdouble qda,qdb;
  qda.v[0] = 1.0; qda.v[1] = 0.0; qda.v[2] = 0.0; qda.v[3] = 0.0;
  qdb.v[0] = M_PI; qdb.v[1] = 0.0; qdb.v[2] = 0.0; qdb.v[3] = 0.0;

  // Turn off internal 80-bit precision on x86 CPU/FPU
/*  volatile unsigned short cw, new_cw;
  _FPU_GETCW(cw);
  new_cw = (cw & ~_FPU_EXTENDED) | _FPU_DOUBLE;
  _FPU_SETCW(new_cw);
*/  
  XPFPA_DECLARE()
  XPFPA_SWITCH_DOUBLE()

  qdouble erg = qddiv(qda,qdb);
//  Rprintf("%23.20e %23.20e %23.20e %23.20e\n",erg.v[0],erg.v[1],erg.v[2],erg.v[3]);
  if((erg.v[0]!=  3.18309886183790702319e-01)||
     (erg.v[1]!= -7.27041035184147510506e-18)||
     (erg.v[2]!=  1.14870604126263509137e-34)||
     (erg.v[3]!= -7.01498859463270463266e-51)) 
     warning("CPU does not follow IEEE double precision arithmetics (although we tried to force this).\n  ppearson/qpearson NOT reliable for type IV distributions!\n  Please report this to m@rtinbecker.de!\n"); 
  // Restore internal precision on x86 CPU/FPU
//  _FPU_SETCW(cw);
  XPFPA_RESTORE()
  return(0.0);
}
