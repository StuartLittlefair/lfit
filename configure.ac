dnl Process this file with autoconf to produce a configure script.

AC_INIT([lfit], [1.0])

AC_CONFIG_SRCDIR([src/lfit.cc])

AM_INIT_AUTOMAKE

dnl Get rid of all the -D stuff
AC_CONFIG_HEADERS([config.h])

AC_PROG_MAKE_SET

dnl Compiler checks
AC_PROG_CXX
AC_PROG_CC

dnl Make C++ the default
AC_LANG_CPLUSPLUS
AC_REQUIRE_CPP

dnl libtool
AC_CONFIG_MACRO_DIR([m4])
AC_PROG_LIBTOOL

dnl check for other trm software

AC_CHECK_HEADERS([trm/subs.h trm/constants.h trm/array1d.h],[],AC_MSG_ERROR([cannot find 'subs' headers]))

AC_CHECK_LIB([subs], [main], [], AC_MSG_ERROR([cannot find 'subs' library]))

dnl Third-party software checks

dnl First headers

AC_CHECK_HEADERS([stdlib.h float.h math.h sstream string fstream iostream], [], AC_MSG_ERROR([missing header; please fix]))

AC_CHECK_HEADERS([slalib.h], [], AC_MSG_ERROR([missing header; please fix]))

dnl Now libraries
AC_CHECK_LIB([csla], [main], [], AC_MSG_ERROR([cannot find the sla C library]))

dnl PGPLOT has its own macro 'cos its a pain
TRM_LIB_PGPLOT

dnl Installation program
AC_PROG_INSTALL

dnl The Makefiles to create
AC_CONFIG_FILES([Makefile include/Makefile src/Makefile])

AC_OUTPUT
