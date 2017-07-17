AC_PREREQ(2.59)
AC_INIT(PearsonDS, 0.91, m@rtinbecker.de)
AC_CONFIG_SRCDIR([src/QuadDouble.c]) ## representative source file
AC_CONFIG_HEADER([src/config.h])     ## Autoheader is used

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
echo "could not determine R_HOME"
exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`

host=`echo 'cat(R.version$platform,"\n")' | "${R_HOME}/bin/R" --vanilla --slave`

dnl --------------------------------------------------------------------
dnl Source of this part: http://www.christian-seiler.de/projekte/fpmath/
dnl --------------------------------------------------------------------

dnl Floating point precision checks
dnl
dnl This file contains floating point precision checks
dnl
dnl This file is released under public domain - or - in countries where this is
dnl not possible under the following license:
dnl
dnl    Permission is hereby granted, free of charge, to any person obtaining a
dnl    copy of this software, to deal in the software without restriction,
dnl    including without limitation the rights to use, copy, modify, merge,
dnl    publish, distribute, sublicense, and/or sell copies of the software,
dnl    and to permit persons to whom the software is furnished to do so, subject
dnl    to no condition whatsoever.
dnl
dnl    This software is provided AS IS, without warranty of any kind, express or
dnl    implied.

dnl FIXME:
dnl  This file was only built for x86 and x86_64 platforms but it does not check
dnl  for the platform since CMake does not provide a viable variable.

AC_DEFUN([CHECK_FLOAT_PRECISION],[
  AC_MSG_CHECKING([for usable _FPU_SETCW])
  AC_TRY_RUN([
   #include <stdio.h>
   #include <string.h>
   #include <fpu_control.h>
  
   double div (double a, double b) {
     fpu_control_t fpu_oldcw, fpu_cw;
     volatile double result;
  
     _FPU_GETCW(fpu_oldcw);
     fpu_cw = (fpu_oldcw & ~_FPU_EXTENDED & ~_FPU_SINGLE) | _FPU_DOUBLE;
     _FPU_SETCW(fpu_cw);
     result = a / b;
     _FPU_SETCW(fpu_oldcw);
     return result;
   }
  
   int main (int argc, char **argv) {
     double d = div (2877.0, 1000000.0);
     char buf[255];
     sprintf(buf, "%.30f", d);
     // see if the result is actually in double precision
     return strncmp(buf, "0.00287699", 10) == 0 ? 0 : 1;
   }
  ], [ac_cfp_have__fpu_setcw=yes], [ac_cfp_have__fpu_setcw=no], [ac_cfp_have__fpu_setcw=no])
  if test "$ac_cfp_have__fpu_setcw" = "yes" ; then
    AC_DEFINE(HAVE__FPU_SETCW, 1, [whether _FPU_SETCW is present and usable])
    AC_MSG_RESULT(yes)
  else
    AC_MSG_RESULT(no)
  fi
  
  AC_MSG_CHECKING([for usable fpsetprec])
  AC_TRY_RUN([
   #include <stdio.h>
   #include <string.h>
   #include <machine/ieeefp.h>
  
   double div (double a, double b) {
     fp_prec_t fpu_oldprec;
     volatile double result;
  
     fpu_oldprec = fpgetprec();
     fpsetprec(FP_PD);
     result = a / b;
     fpsetprec(fpu_oldprec);
     return result;
   }
  
   int main (int argc, char **argv) {
     double d = div (2877.0, 1000000.0);
     char buf[255];
     sprintf(buf, "%.30f", d);
     // see if the result is actually in double precision
     return strncmp(buf, "0.00287699", 10) == 0 ? 0 : 1;
   }
  ], [ac_cfp_have_fpsetprec=yes], [ac_cfp_have_fpsetprec=no], [ac_cfp_have_fpsetprec=no])
  if test "$ac_cfp_have_fpsetprec" = "yes" ; then
    AC_DEFINE(HAVE_FPSETPREC, 1, [whether fpsetprec is present and usable])
    AC_MSG_RESULT(yes)
  else
    AC_MSG_RESULT(no)
  fi

  AC_MSG_CHECKING([for usable _controlfp])
  AC_TRY_RUN([
   #include <stdio.h>
   #include <string.h>
   #include <float.h>
  
   double div (double a, double b) {
     unsigned int fpu_oldcw;
     volatile double result;
  
     fpu_oldcw = _controlfp(0, 0);
     _controlfp(_PC_53, _MCW_PC);
     result = a / b;
     _controlfp(fpu_oldcw, _MCW_PC);
     return result;
   }
  
   int main (int argc, char **argv) {
     double d = div (2877.0, 1000000.0);
     char buf[255];
     sprintf(buf, "%.30f", d);
     // see if the result is actually in double precision
     return strncmp(buf, "0.00287699", 10) == 0 ? 0 : 1;
   }
  ], [ac_cfp_have__controlfp=yes], [ac_cfp_have__controlfp=no], [ac_cfp_have__controlfp=no])
  if test "$ac_cfp_have__controlfp" = "yes" ; then
    AC_DEFINE(HAVE__CONTROLFP, 1, [whether _controlfp is present usable])
    AC_MSG_RESULT(yes)
  else
    AC_MSG_RESULT(no)
  fi

  AC_MSG_CHECKING([for usable _controlfp_s])
  AC_TRY_RUN([
   #include <stdio.h>
   #include <string.h>
   #include <float.h>
  
   double div (double a, double b) {
     unsigned int fpu_oldcw, fpu_cw;
     volatile double result;
  
     _controlfp_s(&fpu_cw, 0, 0);
     fpu_oldcw = fpu_cw;
     _controlfp_s(&fpu_cw, _PC_53, _MCW_PC);
     result = a / b;
     _controlfp_s(&fpu_cw, fpu_oldcw, _MCW_PC);
     return result;
   }
  
   int main (int argc, char **argv) {
     double d = div (2877.0, 1000000.0);
     char buf[255];
     sprintf(buf, "%.30f", d);
     // see if the result is actually in double precision
     return strncmp(buf, "0.00287699", 10) == 0 ? 0 : 1;
   }
  ], [ac_cfp_have__controlfp_s=yes], [ac_cfp_have__controlfp_s=no], [ac_cfp_have__controlfp_s=no])
  if test "$ac_cfp_have__controlfp_s" = "yes" ; then
    AC_DEFINE(HAVE__CONTROLFP_S, 1, [whether _controlfp_s is present and usable])
    AC_MSG_RESULT(yes)
  else
    AC_MSG_RESULT(no)
  fi

  AC_MSG_CHECKING([whether FPU control word can be manipulated by inline assembler])
  AC_TRY_RUN([
   #include <stdio.h>
   #include <string.h>
   
   double div (double a, double b) {
     unsigned int oldcw, cw;
     volatile double result;
     
     __asm__ __volatile__ ("fnstcw %0" : "=m" (*&oldcw));
     cw = (oldcw & ~0x0 & ~0x300) | 0x200;
     __asm__ __volatile__ ("fldcw %0" : : "m" (*&cw));
     
     result = a / b;
     
     __asm__ __volatile__ ("fldcw %0" : : "m" (*&oldcw));
     
     return result;
   }
  
   int main (int argc, char **argv) {
     double d = div (2877.0, 1000000.0);
     char buf[255];
     sprintf(buf, "%.30f", d);
     // see if the result is actually in double precision
     return strncmp(buf, "0.00287699", 10) == 0 ? 0 : 1;
   }
  ], [ac_cfp_have_fpu_inline_asm_x86=yes], [ac_cfp_have_fpu_inline_asm_x86=no], [ac_cfp_have_fpu_inline_asm_x86=no])
  if test "$ac_cfp_have_fpu_inline_asm_x86" = "yes" ; then
    AC_DEFINE(HAVE_FPU_INLINE_ASM_X86, 1, [whether FPU control word can be manipulated by inline assembler])
    AC_MSG_RESULT(yes)
  else
    AC_MSG_RESULT(no)
  fi
])

dnl -----------
dnl End of part
dnl -----------

# Check for x86 FPU fix
AC_MSG_CHECKING([whether x86 FPU control code is needed])
x86_fix="no"
case "$host" in
  i?86-*-* | k?*-*-* | athlon-*-* | x86_64-*-*)
    AC_DEFINE([X86], [1], [Whether to use x86 fpu fix.])
    x86_fix="yes"
  ;;
esac
AC_MSG_RESULT($x86_fix)
if test "$x86_fix" = "yes"; then
  CHECK_FLOAT_PRECISION
fi

AC_OUTPUT
