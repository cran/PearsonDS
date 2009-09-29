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
  sinclude(float_precision.m4)
  CHECK_FLOAT_PRECISION
fi

AC_OUTPUT
