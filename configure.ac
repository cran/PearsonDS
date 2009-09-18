AC_PREREQ(2.59)
AC_INIT(PearsonDS, 0.90, m@rtinbecker.de)
AC_CONFIG_SRCDIR([src/QuadDouble.c]) ## representative source file
AC_CONFIG_HEADER([src/config.h])     ## Autoheader is used

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
