#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>
#include <Rmath.h>
#include "pearson.h"

double StopCritQD (QDcomplex qd1, QDcomplex qd2) {
  return(fmax2(fabs(qd2d(qd1.r)),fabs(qd2d(qd1.i)))/
         fmax2(fabs(qd2d(qd2.r)),fabs(qd2d(qd2.i))));
}

double QDnorm2 (QDcomplex qd1) {
  double re = qd2d(qd1.r);
  double im = qd2d(qd1.i); 
  return(re*re+im*im);
}

SEXP F21QD(SEXP A, SEXP B, SEXP C, SEXP Z, SEXP Minit, SEXP Maxit) {
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
  QDcomplex currsum,tres;
  SEXP LRes, LNames, Res, Rel;
  PROTECT (LRes   = allocVector(VECSXP, 2));
  PROTECT (LNames = allocVector(STRSXP, 2));
  PROTECT (Res    = allocVector(CPLXSXP, n));
  PROTECT (Rel    = allocVector(REALSXP, n));
  Rcomplex *res = COMPLEX(Res);
  double   *rel = REAL(Rel);
  for (int i=0; i<n; i++) {
    curra = a; currb = b; currc = c;
    tres.r.v[0] = 1.; tres.r.v[1] = 0.; tres.r.v[2] = 0.; tres.r.v[3] = 0.;
    tres.i.v[0] = 0.; tres.i.v[1] = 0.; tres.i.v[2] = 0.; tres.i.v[3] = 0.;
    currsum = tres;
    maxsum  = 1.;
    for (f = 1.; (f<minit)||((f<maxit)&&(StopCritQD(currsum,tres)>DOUBLE_EPS)); f=f+1.) {
      R_CheckUserInterrupt();
      currsum = QDMultC(currsum,curra);
      currsum = QDMultC(currsum,currb);
      currsum = QDDivC(currsum,currc);
      currsum = QDMultC(currsum,z[i]);
      currsum = QDDivR(currsum,f);
      tres    = QDAdd(tres,currsum); 
      curra   = CAdd1(curra);
      currb   = CAdd1(currb);
      currc   = CAdd1(currc);
//      Rprintf("%f: %g + %g i\n",f,qd2d(currsum.r),qd2d(currsum.i));
      maxsum  = fmax2(maxsum,QDnorm2(currsum));
    }
    if (f>=maxit) {
//      Rprintf("Appr: %f - Z: %f + %f i, Currsum; %f + %f i, Rel: %g\n",f,(z[i].r),(z[i].i),qd2d(currsum.r),qd2d(currsum.i),StopCritQD(currsum,tres));
      warning("approximation of hypergeometric function inexact");
    }  
    res[i].r = qd2d(tres.r);
    res[i].i = qd2d(tres.i);
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

SEXP F21QDa1cR(SEXP A, SEXP B, SEXP C, SEXP Z, SEXP Minit, SEXP Maxit) {
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
  QDcomplex currsum,tres;
  SEXP LRes, LNames, Res, Rel;
  PROTECT (LRes   = allocVector(VECSXP, 2));
  PROTECT (LNames = allocVector(STRSXP, 2));
  PROTECT (Res    = allocVector(CPLXSXP, n));
  PROTECT (Rel    = allocVector(REALSXP, n));
  Rcomplex *res = COMPLEX(Res);
  double   *rel = REAL(Rel);
  for (int i=0; i<n; i++) {
    currb = b; currc = c;
    tres.r.v[0] = 1.; tres.r.v[1] = 0.; tres.r.v[2] = 0.; tres.r.v[3] = 0.;
    tres.i.v[0] = 0.; tres.i.v[1] = 0.; tres.i.v[2] = 0.; tres.i.v[3] = 0.;
    currsum = tres;
    maxsum  = 1.;
    for (f = 1.; (f<minit)||((f<maxit)&&(StopCritQD(currsum,tres)>DOUBLE_EPS)); f=f+1.) {
      R_CheckUserInterrupt();
      currsum = QDMultC(currsum,currb);
      currsum = QDDivR(currsum,currc);
      currsum = QDMultC(currsum,z[i]);
      tres    = QDAdd(tres,currsum); 
      currb   = CAdd1(currb);
      currc   = currc+1.;
//      Rprintf("%f: %g + %g i\n",f,qd2d(currsum.r),qd2d(currsum.i));
      maxsum  = fmax2(maxsum,QDnorm2(currsum));
    }
    if (f>=maxit) {
//      Rprintf("Appr: %f - Z: %f + %f i, Currsum; %f + %f i, Rel: %g\n",f,(z[i].r),(z[i].i),qd2d(currsum.r),qd2d(currsum.i),StopCritQD(currsum,tres));
      warning("approximation of hypergeometric function inexact");
    }  
    res[i].r = qd2d(tres.r);
    res[i].i = qd2d(tres.i);
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

SEXP F21QDa1bR(SEXP A, SEXP B, SEXP C, SEXP Z, SEXP Minit, SEXP Maxit) {
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
  QDcomplex currsum,tres;
  SEXP LRes, LNames, Res, Rel;
  PROTECT (LRes   = allocVector(VECSXP, 2));
  PROTECT (LNames = allocVector(STRSXP, 2));
  PROTECT (Res    = allocVector(CPLXSXP, n));
  PROTECT (Rel    = allocVector(REALSXP, n));
  Rcomplex *res = COMPLEX(Res);
  double   *rel = REAL(Rel);
  for (int i=0; i<n; i++) {
    currb = b; currc = c;
    tres.r.v[0] = 1.; tres.r.v[1] = 0.; tres.r.v[2] = 0.; tres.r.v[3] = 0.;
    tres.i.v[0] = 0.; tres.i.v[1] = 0.; tres.i.v[2] = 0.; tres.i.v[3] = 0.;
    currsum = tres;
    maxsum  = 1.;
    for (f = 1.; (f<minit)||((f<maxit)&&(StopCritQD(currsum,tres)>DOUBLE_EPS)); f=f+1.) {
      R_CheckUserInterrupt();
      currsum = QDMultR(currsum,currb);
      currsum = QDDivC(currsum,currc);
      currsum = QDMultC(currsum,z[i]);
      tres    = QDAdd(tres,currsum); 
      currb   = currb+1.;
      currc   = CAdd1(currc);
//      Rprintf("%f: %g + %g i\n",f,qd2d(currsum.r),qd2d(currsum.i));
      maxsum  = fmax2(maxsum,QDnorm2(currsum));
    }
    if (f>=maxit) {
//      Rprintf("Appr: %f - Z: %f + %f i, Currsum; %f + %f i, Rel: %g\n",f,(z[i].r),(z[i].i),qd2d(currsum.r),qd2d(currsum.i),StopCritQD(currsum,tres));
      warning("approximation of hypergeometric function inexact");
    }  
    res[i].r = qd2d(tres.r);
    res[i].i = qd2d(tres.i);
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

SEXP F21QDaR(SEXP A, SEXP B, SEXP C, SEXP Z, SEXP Minit, SEXP Maxit) {
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
  QDcomplex currsum,tres;
  SEXP LRes, LNames, Res, Rel;
  PROTECT (LRes   = allocVector(VECSXP, 2));
  PROTECT (LNames = allocVector(STRSXP, 2));
  PROTECT (Res    = allocVector(CPLXSXP, n));
  PROTECT (Rel    = allocVector(REALSXP, n));
  Rcomplex *res = COMPLEX(Res);
  double   *rel = REAL(Rel);
  for (int i=0; i<n; i++) {
    curra = a; currb = b; currc = c;
    tres.r.v[0] = 1.; tres.r.v[1] = 0.; tres.r.v[2] = 0.; tres.r.v[3] = 0.;
    tres.i.v[0] = 0.; tres.i.v[1] = 0.; tres.i.v[2] = 0.; tres.i.v[3] = 0.;
    currsum = tres;
    maxsum  = 1.;
    for (f = 1.; (f<minit)||((f<maxit)&&(StopCritQD(currsum,tres)>DOUBLE_EPS)); f=f+1.) {
      R_CheckUserInterrupt();
      currsum = QDMultR(currsum,curra);
      currsum = QDMultC(currsum,currb);
      currsum = QDDivC(currsum,currc);
      currsum = QDMultC(currsum,z[i]);
      currsum = QDDivR(currsum,f);
      tres    = QDAdd(tres,currsum); 
      curra   = curra+1.;
      currb   = CAdd1(currb);
      currc   = CAdd1(currc);
//      Rprintf("%f: %g + %g i\n",f,qd2d(currsum.r),qd2d(currsum.i));
      maxsum  = fmax2(maxsum,QDnorm2(currsum));
    }
    if (f>=maxit) {
//      Rprintf("Appr: %f - Z: %f + %f i, Currsum; %f + %f i, Rel: %g\n",f,(z[i].r),(z[i].i),qd2d(currsum.r),qd2d(currsum.i),StopCritQD(currsum,tres));
      warning("approximation of hypergeometric function inexact");
    }  
    res[i].r = qd2d(tres.r);
    res[i].i = qd2d(tres.i);
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
