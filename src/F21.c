#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>
#include <Rmath.h>
#include "pearson.h"

double StopCritD (Rcomplex qd1, Rcomplex qd2) {
  return(fmax2(fabs(qd1.r),fabs(qd1.i))/
         fmax2(fabs(qd2.r),fabs(qd2.i)));
}

SEXP F21D(SEXP A, SEXP B, SEXP C, SEXP Z, SEXP Minit, SEXP Maxit) {
  int    n     = LENGTH(Z);
  double maxit = REAL(Maxit)[0];
  double minit = REAL(Minit)[0];
  double f, maxsum;
  Rcomplex  a  = COMPLEX(AS_COMPLEX(A))[0];
  Rcomplex  b  = COMPLEX(AS_COMPLEX(B))[0];
  Rcomplex  c  = COMPLEX(AS_COMPLEX(C))[0];
  Rcomplex *z  = COMPLEX(Z);
  Rcomplex currc,currb,curra,currsum,tres;
  SEXP LRes, LNames, Res, Rel;
  PROTECT (LRes   = allocVector(VECSXP, 2));
  PROTECT (LNames = allocVector(STRSXP, 2));
  PROTECT (Res    = allocVector(CPLXSXP, n));
  PROTECT (Rel    = allocVector(REALSXP, n));
  Rcomplex *res = COMPLEX(Res);
  double   *rel = REAL(Rel);
  for (int i=0; i<n; i++) {
    curra = a; currb = b; currc = c; currsum.r = 1.; currsum.i = 0.;
    tres  = currsum; maxsum = 1.;
    for (f = 1.; (f<minit)||((f<maxit)&&(StopCritD(currsum,tres)>DOUBLE_EPS)); f=f+1.) {
      R_CheckUserInterrupt();
      currsum = CMult(currsum,curra);
      currsum = CMult(currsum,currb);
      currsum = CDiv(currsum,currc);
      currsum = CMult(currsum,z[i]);
      currsum = CMultR(currsum,1./f);
      tres    = CAdd(tres,currsum);
      curra   = CAdd1(curra);
      currb   = CAdd1(currb);
      currc   = CAdd1(currc);
//      Rprintf("%f: %g + %g i\n",f,currsum.r,currsum.i);
      maxsum  = fmax2(maxsum,Cabs2(currsum));
    }
    if (f>=maxit) {
//      Rprintf("D:Appr: %f - Z: %f + %f i, Currsum; %f + %f i, Rel: %g\n",f,z[i].r,z[i].i,currsum.r,currsum.i,StopCritD(currsum,tres));
      warning("approximation of hypergeometric function inexact");
    }  
    res[i] = tres;
    rel[i] = sqrt(Cabs2(res[i])/maxsum);
//    Rprintf("Iterations: %f, Result: %g+%g i\n",f,res[i].r,res[i].i);
  }
  
  SET_VECTOR_ELT(LRes, 0, Res);
  SET_STRING_ELT(LNames, 0, mkChar("value"));
  SET_VECTOR_ELT(LRes, 1, Rel);
  SET_STRING_ELT(LNames, 1, mkChar("rel"));
  setAttrib(LRes, R_NamesSymbol, LNames);
  UNPROTECT(4);
  return(LRes);
}

SEXP F21Da1cR(SEXP A, SEXP B, SEXP C, SEXP Z, SEXP Minit, SEXP Maxit) {
  int    n     = LENGTH(Z);
  double maxit = REAL(Maxit)[0];
  double minit = REAL(Minit)[0];
  double f, maxsum;
  Rcomplex  b  = COMPLEX(AS_COMPLEX(B))[0];
  double    c  = REAL(C)[0];
  Rcomplex *z  = COMPLEX(Z);
  double   currc; 
  Rcomplex currb,currsum,tres;
  SEXP LRes, LNames, Res, Rel;
  PROTECT (LRes   = allocVector(VECSXP, 2));
  PROTECT (LNames = allocVector(STRSXP, 2));
  PROTECT (Res    = allocVector(CPLXSXP, n));
  PROTECT (Rel    = allocVector(REALSXP, n));
  Rcomplex *res = COMPLEX(Res);
  double   *rel = REAL(Rel);
  for (int i=0; i<n; i++) {
    currb = b; currc = c; currsum.r = 1.; currsum.i = 0.;
    tres  = currsum; maxsum = 1.;
    for (f = 1.; (f<minit)||((f<maxit)&&(StopCritD(currsum,tres)>DOUBLE_EPS)); f=f+1.) {
      R_CheckUserInterrupt();
      currsum = CMult(currsum,currb);
      currsum = CDivR(currsum,currc);
      currsum = CMult(currsum,z[i]);
      tres    = CAdd(tres,currsum);
      currb   = CAdd1(currb);
      currc   = currc+1.;
//      Rprintf("%f: %g + %g i\n",f,currsum.r,currsum.i);
      maxsum  = fmax2(maxsum,Cabs2(currsum));
    }
    if (f>=maxit) {
//      Rprintf("D:Appr: %f - Z: %f + %f i, Currsum; %f + %f i, Rel: %g\n",f,z[i].r,z[i].i,currsum.r,currsum.i,StopCritD(currsum,tres));
      warning("approximation of hypergeometric function inexact");
    }  
    res[i] = tres;
    rel[i] = sqrt(Cabs2(res[i])/maxsum);
//    Rprintf("Iterations: %f, Result: %g+%g i\n",f,res[i].r,res[i].i);
  }
  
  SET_VECTOR_ELT(LRes, 0, Res);
  SET_STRING_ELT(LNames, 0, mkChar("value"));
  SET_VECTOR_ELT(LRes, 1, Rel);
  SET_STRING_ELT(LNames, 1, mkChar("rel"));
  setAttrib(LRes, R_NamesSymbol, LNames);
  UNPROTECT(4);
  return(LRes);
}

SEXP F21Da1bR(SEXP A, SEXP B, SEXP C, SEXP Z, SEXP Minit, SEXP Maxit) {
  int    n     = LENGTH(Z);
  double maxit = REAL(Maxit)[0];
  double minit = REAL(Minit)[0];
  double f, maxsum;
  double    b  = REAL(B)[0];
  Rcomplex  c  = COMPLEX(AS_COMPLEX(C))[0];
  Rcomplex *z  = COMPLEX(Z);
  double   currb;
  Rcomplex currc,currsum,tres;
  SEXP LRes, LNames, Res, Rel;
  PROTECT (LRes   = allocVector(VECSXP, 2));
  PROTECT (LNames = allocVector(STRSXP, 2));
  PROTECT (Res    = allocVector(CPLXSXP, n));
  PROTECT (Rel    = allocVector(REALSXP, n));
  Rcomplex *res = COMPLEX(Res);
  double   *rel = REAL(Rel);
  for (int i=0; i<n; i++) {
    currb = b; currc = c; currsum.r = 1.; currsum.i = 0.;
    tres  = currsum; maxsum = 1.;
    for (f = 1.; (f<minit)||((f<maxit)&&(StopCritD(currsum,tres)>DOUBLE_EPS)); f=f+1.) {
      R_CheckUserInterrupt();
      currsum = CMultR(currsum,currb);
      currsum = CDiv(currsum,currc);
      currsum = CMult(currsum,z[i]);
      tres    = CAdd(tres,currsum);
      currb   = currb+1.;
      currc   = CAdd1(currc);
//      Rprintf("%f: %g + %g i\n",f,currsum.r,currsum.i);
      maxsum  = fmax2(maxsum,Cabs2(currsum));
    }
    if (f>=maxit) {
//      Rprintf("D:Appr: %f - Z: %f + %f i, Currsum; %f + %f i, Rel: %g\n",f,z[i].r,z[i].i,currsum.r,currsum.i,StopCritD(currsum,tres));
      warning("approximation of hypergeometric function inexact");
    }  
    res[i] = tres;
    rel[i] = sqrt(Cabs2(res[i])/maxsum);
//    Rprintf("Iterations: %f, Result: %g+%g i\n",f,res[i].r,res[i].i);
  }
  
  SET_VECTOR_ELT(LRes, 0, Res);
  SET_STRING_ELT(LNames, 0, mkChar("value"));
  SET_VECTOR_ELT(LRes, 1, Rel);
  SET_STRING_ELT(LNames, 1, mkChar("rel"));
  setAttrib(LRes, R_NamesSymbol, LNames);
  UNPROTECT(4);
  return(LRes);
}

SEXP F21DaR(SEXP A, SEXP B, SEXP C, SEXP Z, SEXP Minit, SEXP Maxit) {
  int    n     = LENGTH(Z);
  double maxit = REAL(Maxit)[0];
  double minit = REAL(Minit)[0];
  double f, maxsum;
  double    a  = REAL(A)[0];
  Rcomplex  b  = COMPLEX(AS_COMPLEX(B))[0];
  Rcomplex  c  = COMPLEX(AS_COMPLEX(C))[0];
  Rcomplex *z  = COMPLEX(Z);
  double   curra;
  Rcomplex currc,currb,currsum,tres;
  SEXP LRes, LNames, Res, Rel;
  PROTECT (LRes   = allocVector(VECSXP, 2));
  PROTECT (LNames = allocVector(STRSXP, 2));
  PROTECT (Res    = allocVector(CPLXSXP, n));
  PROTECT (Rel    = allocVector(REALSXP, n));
  Rcomplex *res = COMPLEX(Res);
  double   *rel = REAL(Rel);
  for (int i=0; i<n; i++) {
    curra = a; currb = b; currc = c; currsum.r = 1.; currsum.i = 0.;
    tres  = currsum; maxsum = 1.;
    for (f = 1.; (f<minit)||((f<maxit)&&(StopCritD(currsum,tres)>DOUBLE_EPS)); f=f+1.) {
      R_CheckUserInterrupt();
      currsum = CMultR(currsum,curra);
      currsum = CMult(currsum,currb);
      currsum = CDiv(currsum,currc);
      currsum = CMult(currsum,z[i]);
      currsum = CDivR(currsum,f);
      tres    = CAdd(tres,currsum);
      curra   = curra+1.;
      currb   = CAdd1(currb);
      currc   = CAdd1(currc);
//      Rprintf("%f: %g + %g i\n",f,currsum.r,currsum.i);
      maxsum  = fmax2(maxsum,Cabs2(currsum));
    }
    if (f>=maxit) {
//      Rprintf("D:Appr: %f - Z: %f + %f i, Currsum; %f + %f i, Rel: %g\n",f,z[i].r,z[i].i,currsum.r,currsum.i,StopCritD(currsum,tres));
      warning("approximation of hypergeometric function inexact");
    }  
    res[i] = tres;
    rel[i] = sqrt(Cabs2(res[i])/maxsum);
//    Rprintf("Iterations: %f, Result: %g+%g i\n",f,res[i].r,res[i].i);
  }
  
  SET_VECTOR_ELT(LRes, 0, Res);
  SET_STRING_ELT(LNames, 0, mkChar("value"));
  SET_VECTOR_ELT(LRes, 1, Rel);
  SET_STRING_ELT(LNames, 1, mkChar("rel"));
  setAttrib(LRes, R_NamesSymbol, LNames);
  UNPROTECT(4);
  return(LRes);
}
