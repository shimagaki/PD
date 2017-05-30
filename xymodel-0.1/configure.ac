#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT([xymodel],[0.1],[tobias.rautenkranz@stud.unibas.ch])
#AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([xymodel.c])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_AWK
AC_PROG_CC
AC_PROG_MAKE_SET

AC_GNU_SOURCE

# Checks for libraries.
AC_CHECK_LIB([curses], [setupterm], [
		HAVE_LIBCURSES=1
		AC_DEFINE(HAVE_LIBCURSES, 1) ],
		AC_MSG_WARN([*** progress: no progress bar fit to line with ***]))
AM_CONDITIONAL([HAVE_LIBCURSES], [test x$HAVE_LIBCURSES = x1])

AC_CHECK_LIB([glut], [glutInit], [
		HAVE_LIBGLUT=1
		AC_DEFINE(HAVE_LIBGLUT,1)],
		AC_MSG_WARN([*** glvec can't be compiled ***]))

AM_CONDITIONAL([HAVE_LIBGLUT], [test x$HAVE_LIBGLUT = x1])

AC_CHECK_LIB([dl], [dlopen])

AC_CHECK_LIB([m], [cos])

# Checks for header files.
	AC_HEADER_STDC
AC_CHECK_HEADERS([fcntl.h stdlib.h string.h sys/time.h unistd.h wchar.h])

# Checks for typedefs, structures, and compiler characteristics.
	AC_C_CONST
	AC_C_INLINE
	AC_TYPE_SIZE_T
	AC_HEADER_TIME

# Checks for library functions.
	AC_FUNC_ERROR_AT_LINE
	AC_FUNC_MALLOC
	AC_FUNC_STRTOD
	AC_CHECK_FUNCS([floor gettimeofday memset pow sqrt strtol])
AC_CHECK_FUNC([getline], , AC_MSG_ERROR([*** getline() function needed (GNU extension)***]))

AC_CONFIG_FILES([Makefile])
	AC_OUTPUT
