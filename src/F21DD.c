#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>
#include <Rmath.h>
#include "pearson.h"

double StopCritDD (DDcomplex dd1, DDcomplex dd2) {
  return(fmax2(fabs(dd2d(dd1.r)),fabs(dd2d(dd1.i)))/
         fmax2(fabs(dd2d(dd2.r)),fabs(dd2d(dd2.i))));
}

double DDnorm2 (DDcomplex dd1) {
  double re = dd2d(dd1.r);
  double im = dd2d(dd1.i); 
  return(re*re+im*im);
}

SEXP F21DD(SEXP A, SEXP B, SEXP C, SEXP Z, SEXP Minit, SEXP Maxit) {
  // Turn off internal 80-bit precision on x86 CPU/FPU
/*  volatile unsigned short cw, new_cw;
  _FPU_GETCW(cw);
  new_cw = (cw & ~_FPU_EXTENDED) | _FPU_DOUBLE;
  _FPU_SETCW(new_cw); */
  XPFPA_DECLARE()
  XPFPA_SWITCH_DOUBLE()
  int    n     = LENGTH(Z);
  double maxit = REAL(Maxit)[0];
  double minit = REAL(Minit)[0];
  double f, maxsum;
  Rcomplex  a  = COMPLEX(AS_COMPLEX(A))[0];
  Rcomplex  b  = COMPLEX(AS_COMPLEX(B))[0];
  Rcomplex  c  = COMPLEX(AS_COMPLEX(C))[0];
  Rcomplex *z  = COMPLEX(Z);
  Rcomplex currc,currb,curra;
  DDcomplex currsum,tres;
  SEXP LRes, LNames, Res, Rel;
  PROTECT (LRes   = allocVector(VECSXP, 2));
  PROTECT (LNames = allocVector(STRSXP, 2));
  PROTECT (Res    = allocVector(CPLXSXP, n));
  PROTECT (Rel    = allocVector(REALSXP, n));
  Rcomplex *res = COMPLEX(Res);
  double   *rel = REAL(Rel);
  for (int i=0; i<n; i++) {
    curra = a; currb = b; currc = c;
    tres.r.r = 1.; tres.r.e = 0.;
    tres.i.r = 0.; tres.i.e = 0.;
    currsum = tres;
    maxsum  = 1.;
    for (f = 1.; (f<minit)||((f<maxit)&&(StopCritDD(currsum,tres)>DOUBLE_EPS)); f=f+1.) {
      R_CheckUserInterrupt();
      currsum = DDMultC(currsum,curra);
      currsum = DDMultC(currsum,currb);
      currsum = DDDivC(currsum,currc);
      currsum = DDMultC(currsum,z[i]);
      currsum = DDDivR(currsum,f);
      tres    = DDAdd(tres,currsum); 
      curra   = CAdd1(curra);
      currb   = CAdd1(currb);
      currc   = CAdd1(currc);
//      Rprintf("%f: %g + %g i\n",f,qd2d(currsum.r),qd2d(currsum.i));
      maxsum  = fmax2(maxsum,DDnorm2(currsum));
    }
    if (f>=maxit) {
//      Rprintf("Appr: %f - Z: %f + %f i, Currsum; %f + %f i, Rel: %g\n",f,(z[i].r),(z[i].i),qd2d(currsum.r),qd2d(currsum.i),StopCritQD(currsum,tres));
      warning("approximation of hypergeometric function inexact");
    }  
    res[i].r = dd2d(tres.r);
    res[i].i = dd2d(tres.i);
    rel[i]   = sqrt(Cabs2(res[i])/maxsum);
//    Rprintf("Iterations: %f, Result: %g+%g i\n",f,res[i].r,res[i].i);
  }
  
  // Restore internal precision on x86 CPU/FPU
//  _FPU_SETCW(cw);
  XPFPA_RESTORE()
  SET_VECTOR_ELT(LRes, 0, Res);
  SET_STRING_ELT(LNames, 0, mkChar("value"));
  SET_VECTOR_ELT(LRes, 1, Rel);
  SET_STRING_ELT(LNames, 1, mkChar("rel"));
  setAttrib(LRes, R_NamesSymbol, LNames);
  UNPROTECT(4);
  return(LRes);
}

SEXP F21DDa1cR(SEXP A, SEXP B, SEXP C, SEXP Z, SEXP Minit, SEXP Maxit) {
  // Turn off internal 80-bit precision on x86 CPU/FPU
/*  volatile unsigned short cw, new_cw;
  _FPU_GETCW(cw);
  new_cw = (cw & ~_FPU_EXTENDED) | _FPU_DOUBLE;
  _FPU_SETCW(new_cw);
*/
  XPFPA_DECLARE()
  XPFPA_SWITCH_DOUBLE()
  int    n     = LENGTH(Z);
  double maxit = REAL(Maxit)[0];
  double minit = REAL(Minit)[0];
  double f, maxsum;
  Rcomplex  b  = COMPLEX(AS_COMPLEX(B))[0];
  double    c  = REAL(C)[0];
  Rcomplex *z  = COMPLEX(Z);
  double   currc;
  Rcomplex currb;
  DDcomplex currsum,tres;
  SEXP LRes, LNames, Res, Rel;
  PROTECT (LRes   = allocVector(VECSXP, 2));
  PROTECT (LNames = allocVector(STRSXP, 2));
  PROTECT (Res    = allocVector(CPLXSXP, n));
  PROTECT (Rel    = allocVector(REALSXP, n));
  Rcomplex *res = COMPLEX(Res);
  double   *rel = REAL(Rel);
  for (int i=0; i<n; i++) {
    currb = b; currc = c;
    tres.r.r = 1.; tres.r.e = 0.;
    tres.i.r = 0.; tres.i.e = 0.;
    currsum = tres;
    maxsum  = 1.;
    for (f = 1.; (f<minit)||((f<maxit)&&(StopCritDD(currsum,tres)>DOUBLE_EPS)); f=f+1.) {
      R_CheckUserInterrupt();
      currsum = DDMultC(currsum,currb);
      currsum = DDDivR(currsum,currc);
      currsum = DDMultC(currsum,z[i]);
      tres    = DDAdd(tres,currsum); 
      currb   = CAdd1(currb);
      currc   = currc+1.;
//      Rprintf("%f: %g + %g i\n",f,qd2d(currsum.r),qd2d(currsum.i));
      maxsum  = fmax2(maxsum,DDnorm2(currsum));
    }
    if (f>=maxit) {
//      Rprintf("Appr: %f - Z: %f + %f i, Currsum; %f + %f i, Rel: %g\n",f,(z[i].r),(z[i].i),qd2d(currsum.r),qd2d(currsum.i),StopCritQD(currsum,tres));
      warning("approximation of hypergeometric function inexact");
    }  
    res[i].r = dd2d(tres.r);
    res[i].i = dd2d(tres.i);
    rel[i]   = sqrt(Cabs2(res[i])/maxsum);
//    Rprintf("Iterations: %f, Result: %g+%g i\n",f,res[i].r,res[i].i);
  }
  
  // Restore internal precision on x86 CPU/FPU
//  _FPU_SETCW(cw);
  XPFPA_RESTORE()
  SET_VECTOR_ELT(LRes, 0, Res);
  SET_STRING_ELT(LNames, 0, mkChar("value"));
  SET_VECTOR_ELT(LRes, 1, Rel);
  SET_STRING_ELT(LNames, 1, mkChar("rel"));
  setAttrib(LRes, R_NamesSymbol, LNames);
  UNPROTECT(4);
  return(LRes);
}

SEXP F21DDa1bR(SEXP A, SEXP B, SEXP C, SEXP Z, SEXP Minit, SEXP Maxit) {
  // Turn off internal 80-bit precision on x86 CPU/FPU
/*  volatile unsigned short cw, new_cw;
  _FPU_GETCW(cw);
  new_cw = (cw & ~_FPU_EXTENDED) | _FPU_DOUBLE;
  _FPU_SETCW(new_cw);
*/
  XPFPA_DECLARE()
  XPFPA_SWITCH_DOUBLE()
  int    n     = LENGTH(Z);
  double maxit = REAL(Maxit)[0];
  double minit = REAL(Minit)[0];
  double f, maxsum;
  double    b  = REAL(B)[0];
  Rcomplex  c  = COMPLEX(AS_COMPLEX(C))[0];
  Rcomplex *z  = COMPLEX(Z);
  double   currb;
  Rcomplex currc;
  DDcomplex currsum,tres;
  SEXP LRes, LNames, Res, Rel;
  PROTECT (LRes   = allocVector(VECSXP, 2));
  PROTECT (LNames = allocVector(STRSXP, 2));
  PROTECT (Res    = allocVector(CPLXSXP, n));
  PROTECT (Rel    = allocVector(REALSXP, n));
  Rcomplex *res = COMPLEX(Res);
  double   *rel = REAL(Rel);
  for (int i=0; i<n; i++) {
    currb = b; currc = c;
    tres.r.r = 1.; tres.r.e = 0.;
    tres.i.r = 0.; tres.i.e = 0.;
    currsum = tres;
    maxsum  = 1.;
    for (f = 1.; (f<minit)||((f<maxit)&&(StopCritDD(currsum,tres)>DOUBLE_EPS)); f=f+1.) {
      R_CheckUserInterrupt();
      currsum = DDMultR(currsum,currb);
      currsum = DDDivC(currsum,currc);
      currsum = DDMultC(currsum,z[i]);
      tres    = DDAdd(tres,currsum); 
      currb   = currb+1.;
      currc   = CAdd1(currc);
//      Rprintf("%f: %g + %g i\n",f,qd2d(currsum.r),qd2d(currsum.i));
      maxsum  = fmax2(maxsum,DDnorm2(currsum));
    }
    if (f>=maxit) {
//      Rprintf("Appr: %f - Z: %f + %f i, Currsum; %f + %f i, Rel: %g\n",f,(z[i].r),(z[i].i),qd2d(currsum.r),qd2d(currsum.i),StopCritQD(currsum,tres));
      warning("approximation of hypergeometric function inexact");
    }  
    res[i].r = dd2d(tres.r);
    res[i].i = dd2d(tres.i);
    rel[i]   = sqrt(Cabs2(res[i])/maxsum);
//    Rprintf("Iterations: %f, Result: %g+%g i\n",f,res[i].r,res[i].i);
  }
  
  // Restore internal precision on x86 CPU/FPU
//  _FPU_SETCW(cw);
  XPFPA_RESTORE()
  SET_VECTOR_ELT(LRes, 0, Res);
  SET_STRING_ELT(LNames, 0, mkChar("value"));
  SET_VECTOR_ELT(LRes, 1, Rel);
  SET_STRING_ELT(LNames, 1, mkChar("rel"));
  setAttrib(LRes, R_NamesSymbol, LNames);
  UNPROTECT(4);
  return(LRes);
}

SEXP F21DDaR(SEXP A, SEXP B, SEXP C, SEXP Z, SEXP Minit, SEXP Maxit) {
  // Turn off internal 80-bit precision on x86 CPU/FPU
/*  volatile unsigned short cw, new_cw;
  _FPU_GETCW(cw);
  new_cw = (cw & ~_FPU_EXTENDED) | _FPU_DOUBLE;
  _FPU_SETCW(new_cw);
*/
  XPFPA_DECLARE()
  XPFPA_SWITCH_DOUBLE()
  int    n     = LENGTH(Z);
  double maxit = REAL(Maxit)[0];
  double minit = REAL(Minit)[0];
  double f, maxsum;
  double    a  = REAL(A)[0];
  Rcomplex  b  = COMPLEX(AS_COMPLEX(B))[0];
  Rcomplex  c  = COMPLEX(AS_COMPLEX(C))[0];
  Rcomplex *z  = COMPLEX(Z); 
  double   curra;
  Rcomplex currc,currb;
  DDcomplex currsum,tres;
  SEXP LRes, LNames, Res, Rel;
  PROTECT (LRes   = allocVector(VECSXP, 2));
  PROTECT (LNames = allocVector(STRSXP, 2));
  PROTECT (Res    = allocVector(CPLXSXP, n));
  PROTECT (Rel    = allocVector(REALSXP, n));
  Rcomplex *res = COMPLEX(Res);
  double   *rel = REAL(Rel);
  for (int i=0; i<n; i++) {
    curra = a; currb = b; currc = c;
    tres.r.r = 1.; tres.r.e = 0.;
    tres.i.r = 0.; tres.i.e = 0.;
    currsum = tres;
    maxsum  = 1.;
    for (f = 1.; (f<minit)||((f<maxit)&&(StopCritDD(currsum,tres)>DOUBLE_EPS)); f=f+1.) {
      R_CheckUserInterrupt();
      currsum = DDMultR(currsum,curra);
      currsum = DDMultC(currsum,currb);
      currsum = DDDivC(currsum,currc);
      currsum = DDMultC(currsum,z[i]);
      currsum = DDDivR(currsum,f);
      tres    = DDAdd(tres,currsum); 
      curra   = curra+1.;
      currb   = CAdd1(currb);
      currc   = CAdd1(currc);
//      Rprintf("%f: %g + %g i\n",f,qd2d(currsum.r),qd2d(currsum.i));
      maxsum  = fmax2(maxsum,DDnorm2(currsum));
    }
    if (f>=maxit) {
//      Rprintf("Appr: %f - Z: %f + %f i, Currsum; %f + %f i, Rel: %g\n",f,(z[i].r),(z[i].i),qd2d(currsum.r),qd2d(currsum.i),StopCritQD(currsum,tres));
      warning("approximation of hypergeometric function inexact");
    }  
    res[i].r = dd2d(tres.r);
    res[i].i = dd2d(tres.i);
    rel[i]   = sqrt(Cabs2(res[i])/maxsum);
//    Rprintf("Iterations: %f, Result: %g+%g i\n",f,res[i].r,res[i].i);
  }
  
  // Restore internal precision on x86 CPU/FPU
//  _FPU_SETCW(cw);
  XPFPA_RESTORE()
  SET_VECTOR_ELT(LRes, 0, Res);
  SET_STRING_ELT(LNames, 0, mkChar("value"));
  SET_VECTOR_ELT(LRes, 1, Rel);
  SET_STRING_ELT(LNames, 1, mkChar("rel"));
  setAttrib(LRes, R_NamesSymbol, LNames);
  UNPROTECT(4);
  return(LRes);
}
