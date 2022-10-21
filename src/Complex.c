#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>
#include <Rmath.h>
#include "pearson.h"

Rcomplex CMult(Rcomplex A, Rcomplex B) {
  Rcomplex C;
  C.r = A.r * B.r - A.i * B.i;
  C.i = A.r * B.i + A.i * B.r;
  return(C);
}

Rcomplex CMultR(Rcomplex A, double B) {
  Rcomplex C;
  C.r = A.r * B;
  C.i = A.i * B;
  return(C);
}

Rcomplex CDivR(Rcomplex A, double B) {
  Rcomplex C;
  C.r = A.r / B;
  C.i = A.i / B;
  return(C);
}

Rcomplex CDiv(Rcomplex A, Rcomplex B) {
  Rcomplex C;
  double ratio, den;
  double abr, abi;
  if( (abr = B.r) < 0) abr = - abr;
  if( (abi = B.i) < 0) abi = - abi;
  if( abr <= abi ) {
  	ratio = B.r / B.i;
  	den   = B.i * (1. + ratio*ratio);
  	C.r   = (A.r*ratio + A.i) / den;
   	C.i   = (A.i*ratio - A.r) / den;
  } else {
  	ratio = B.i / B.r ;
  	den  = B.r * (1. + ratio*ratio);
  	C.r  = (A.r + A.i*ratio) / den;
  	C.i  = (A.i - A.r*ratio) / den;
  }
  return(C);
}

Rcomplex CAdd(Rcomplex A, Rcomplex B) {
  Rcomplex C;
  C.r = A.r + B.r;
  C.i = A.i + B.i;
  return(C);
}

Rcomplex CAddR(Rcomplex A, double B) {
  Rcomplex C;
  C.r = A.r + B;
  C.i = A.i;
  return(C);
}

Rcomplex CAdd1(Rcomplex A) {
  Rcomplex C;
  C.r = A.r + 1.;
  C.i = A.i;
  return(C);
}

double Cabs2(Rcomplex A) {
  return(A.r*A.r+A.i*A.i);
}

double CnormInf(Rcomplex A) {
  return(fmax2(fabs(A.r),fabs(A.i)));
}
