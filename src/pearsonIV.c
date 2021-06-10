/* Code (partially) adapted from
   CDF/MEMO/STATISTICS/PUBLIC/6820
   A Guide to the Pearson Type IV Distribution
   Joel Heinrich—University of Pennsylvania
   December 21, 2004
*/      

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>
#include <Rmath.h>
#include "pearson.h"

double gammar2(double x,double y) {
  /* returns abs(gamma(x+iy)/gamma(x))^2 */
  const double y2=y*y, xmin = (2*y2>10.0) ? 2*y2 : 10.0;
  double r=1, s=1, p=1, f=0;
  while(x<xmin) {
    const double t = y/x++;
    r *= 1 + t*t;
  }
  while (p > s*DOUBLE_EPS) {
    p *= y2 + f*f;
    p /= x++ * ++f;
    s += p;
  }
  return 1.0/(r*s);
}

double loggammar2(double x,double y) {
  /* returns log(abs(gamma(x+iy)/gamma(x))^2) */
  const double y2=y*y, xmin = (2*y2>10.0) ? 2*y2 : 10.0;
  double r=0, s=1, p=1, f=0;
  while(x<xmin) {
    const double t = y/x++;
    r += log(1 + t*t);
  }
  while (p > s*DOUBLE_EPS) {
    p *= y2 + f*f;
    p /= x++ * ++f;
    s += p;
  }
  return -r-log(s);
}

double type4norm(double m,double nu,double a) {
  /* returns k */
  return 0.5*M_2_SQRTPI*gammar2(m,0.5*nu)*exp(lgammafn(m)-lgammafn(m-0.5))/a;
}

double logtype4norm(double m,double nu,double a) {
  /* returns log(k) */
  return -M_LN_SQRT_PI+loggammar2(m,0.5*nu)+lgammafn(m)-lgammafn(m-0.5)-log(a);
}

double rpears4(double m,double nu,double a,double lam) {
  /* returns random Pearson IV deviate */
  const double k=type4norm(m,nu,1.0), b=2*m-2, M=atan(-nu/b);
  const double cosM=b/sqrt(b*b+nu*nu), r=b*log(cosM)-nu*M, rc=exp(-r)/k;
  double x,z;
  do {
    int s=0;
    z = 0;
    if( (x=4*unif_rand()) > 2 ) {
      x -= 2;
      s = 1;
    } 
    if (x > 1) x = 1 - (z=log(x-1)) ;
    x = (s) ? M + rc*x : M - rc*x;
  } while (fabs(x) >= M_PI_2 || z + log(unif_rand()) > b*log(cos(x)) - nu*x - r);
  return a*tan(x) + lam;
}

double rpears4k(double m,double nu,double a,double lam, double k) {             // not used any longer...
  /* returns random Pearson IV deviate */
  const double b=2*m-2, M=atan(-nu/b);
  const double cosM=b/sqrt(b*b+nu*nu), r=b*log(cosM)-nu*M, rc=exp(-r)/k;
  double x,z;
  do {
    int s=0;
    z = 0;
    if( (x=4*unif_rand()) > 2 ) {
      x -= 2;
      s = 1;
    } 
    if (x > 1) x = 1 - (z=log(x-1)) ;
    x = (s) ? M + rc*x : M - rc*x;
  } while (fabs(x) >= M_PI_2 || z + log(unif_rand()) > b*log(cos(x)) - nu*x - r);
  return a*tan(x) + lam;
}

double rpears4logk(double m,double nu,double a,double lam, double logk) {
  /* returns random Pearson IV deviate */
  const double b=2*m-2, M=atan(-nu/b);
  const double cosM=b/sqrt(b*b+nu*nu), r=b*log(cosM)-nu*M, rc=exp(-r-logk);
  double x,z;
  do {
    int s=0;
    z = 0;
    if( (x=4*unif_rand()) > 2 ) {
      x -= 2;
      s = 1;
    } 
    if (x > 1) x = 1 - (z=log(x-1)) ;
    x = (s) ? M + rc*x : M - rc*x;
  } while (fabs(x) >= M_PI_2 || z + log(unif_rand()) > b*log(cos(x)) - nu*x - r);
  return a*tan(x) + lam;
}

SEXP PearsonIVnorm(SEXP M, SEXP Nu, SEXP A) {
  SEXP Res;
  PROTECT(Res = allocVector (REALSXP, 1));
  REAL(Res)[0] = type4norm(REAL(M)[0],REAL(Nu)[0],REAL(A)[0]);
  UNPROTECT(1);
  return(Res);
}

SEXP logPearsonIVnorm(SEXP M, SEXP Nu, SEXP A) {
  SEXP Res;
  PROTECT(Res = allocVector (REALSXP, 1));
  REAL(Res)[0] = logtype4norm(REAL(M)[0],REAL(Nu)[0],REAL(A)[0]);
  UNPROTECT(1);
  return(Res);
}

SEXP rPearsonIV(SEXP N, SEXP M, SEXP Nu, SEXP A, SEXP Lambda) {
  int n = INTEGER(AS_INTEGER(N))[0];
  SEXP Res;
  PROTECT(Res = allocVector(REALSXP, n));
  double *res;
  res = REAL(Res);
  double m      = REAL(M)[0];
  double nu     = REAL(Nu)[0];
  double a      = REAL(A)[0];
  double lambda = REAL(Lambda)[0];
  double k      = type4norm(m,nu,1.0);
  GetRNGstate();
  for(int i=0; i<n; i++) {
    res[i] = rpears4k(m,nu,a,lambda,k);
  }
  PutRNGstate();
  UNPROTECT(1);
  return(Res);
}

SEXP rPearsonIVk(SEXP N, SEXP M, SEXP Nu, SEXP A, SEXP Lambda, SEXP K) {
  int n = INTEGER(AS_INTEGER(N))[0];
  SEXP Res;
  PROTECT(Res = allocVector(REALSXP, n));
  double *res;
  res = REAL(Res);
  double m      = REAL(M)[0];
  double nu     = REAL(Nu)[0];
  double a      = REAL(A)[0];
  double lambda = REAL(Lambda)[0];
  double k      = REAL(K)[0];
  GetRNGstate();
  for(int i=0; i<n; i++) {
    res[i] = rpears4k(m,nu,a,lambda,k);
  }
  PutRNGstate();
  UNPROTECT(1);
  return(Res);
}

SEXP rPearsonIVlog(SEXP N, SEXP M, SEXP Nu, SEXP A, SEXP Lambda) {              // rng using logspace
  int n = INTEGER(AS_INTEGER(N))[0];
  SEXP Res;
  PROTECT(Res = allocVector(REALSXP, n));
  double *res;
  res = REAL(Res);
  double m      = REAL(M)[0];
  double nu     = REAL(Nu)[0];
  double a      = REAL(A)[0];
  double lambda = REAL(Lambda)[0];
  double logk   = logtype4norm(m,nu,1.0);
  GetRNGstate();
  for(int i=0; i<n; i++) {
    res[i] = rpears4logk(m,nu,a,lambda,logk);
  }
  PutRNGstate();
  UNPROTECT(1);
  return(Res);
}

SEXP rPearsonIVlogK(SEXP N, SEXP M, SEXP Nu, SEXP A, SEXP Lambda, SEXP LogK) {  // rng using logspace, given logk
  int n = INTEGER(AS_INTEGER(N))[0];
  SEXP Res;
  PROTECT(Res = allocVector(REALSXP, n));
  double *res;
  res = REAL(Res);
  const double m      = REAL(M)[0];
  const double nu     = REAL(Nu)[0];
  const double a      = REAL(A)[0];
  const double lam    = REAL(Lambda)[0];
  const double logk   = REAL(LogK)[0];
  const double b      = 2*m-2;
  const double Mode   = atan(-nu/b);
  const double cosM   = b/sqrt(b*b+nu*nu);
  const double r      = b*log(cosM)-nu*Mode;
  const double rc     = exp(-r-logk);
  double x,z;
  int s;
  GetRNGstate();
  for(int i=0; i<n; i++) {
    do {
      s = 0;
      z = 0;
      if( (x=4*unif_rand()) > 2 ) {
        x -= 2;
        s = 1;
      } 
      if (x > 1) x = 1 - (z=log(x-1)) ;
      x = (s) ? Mode + rc*x : Mode - rc*x;
    } while (fabs(x) >= M_PI_2 || z + log(unif_rand()) > b*log(cos(x)) - nu*x - r);
    res[i] = a*tan(x) + lam;
  }
  PutRNGstate();
  UNPROTECT(1);
  return(Res);
}

/* Old Code; working, but slower */
/*
SEXP rPearsonIVlogK(SEXP N, SEXP M, SEXP Nu, SEXP A, SEXP Lambda, SEXP LogK) {  // rng using logspace, given logk
  int n = INTEGER(AS_INTEGER(N))[0];
  SEXP Res;
  PROTECT(Res = allocVector(REALSXP, n));
  double *res;
  res = REAL(Res);
  double m      = REAL(M)[0];
  double nu     = REAL(Nu)[0];
  double a      = REAL(A)[0];
  double lambda = REAL(Lambda)[0];
  double logk   = REAL(LogK)[0];
  GetRNGstate();
  for(int i=0; i<n; i++) {
    res[i] = rpears4logk(m,nu,a,lambda,logk);
  }
  PutRNGstate();
  UNPROTECT(1);
  return(Res);
}
*/
