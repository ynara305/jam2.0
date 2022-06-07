AC_INIT(jam, 2.0097)
AC_PREFIX_DEFAULT($HOME)
AC_CONFIG_SRCDIR(jam2/main.cxx)
#AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE(dist-bzip2)
AC_CONFIG_HEADERS([config.h])

prevCXX="$CXX"
AC_PROG_CXX
AC_PROG_CC
if test -n "$prevCXX"; then CXX="$prevCXX"; fi
AC_SUBST(CXXFLAGS)
#CXXFLAGS="-O3 -g -Wall -march=corei7 -pipe -fbounds-check\
# -Wuninitialized"
if test "$GCC" = "yes"; then CXXFLAGS="$CXXFLAGS -Wall"; fi

AC_PROG_INSTALL
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE

# Checks for library functions.
AC_CHECK_FUNCS([pow sqrt])

#AC_SUBST(FFLAGS) # => FFLAGS can be specified as an environment variable
#AC_PROG_F77
#AC_F77_LIBRARY_LDFLAGS

#
# libtool
#
AC_CONFIG_MACRO_DIRS([m4])
LT_INIT
#AC_LIBTOOL_DLOPEN
#AM_PROG_LIBTOOL
#AC_SUBST(LOADLIB_FLAG)
#AC_SUBST(LOADLIB_LIB)
#AC_SUBST(LOADLIB_LDFLAGS)
#AC_CHECK_HEADER(ltdl.h,
#	LOADLIB_FLAG="-DUSE_LOADLIB"
#	LOADLIB_LIB="-lltdl"
#	LOADLIB_LDFLAGS="-export-dynamic")
 

# for  pythia8 library
AC_SUBST(PYTHIALIBS)
PYTHIALIBS="-lpythia8"
LDFLAGS="-L$Pythia8/lib -Wl,-rpath,$Pythia8/lib"
CXXFLAGS="$CXXFLAGS -fPIC -I$Pythia8/include -I.."

#AC_SUBST(FORTLIBS)
#AC_SUBST(WITH_PYTHIA8)
#AC_SUBST(WITHOUT_PYTHIA8)
#AC_MSG_CHECKING(whether to use Hijing)
#AC_ARG_WITH(pythia,
#[ --with-pythia[=PATH] Compile with Pythia8 library.],
#[ case "$withval" in
# no)
#   AC_MSG_RESULT(no)
#   ;;
#  *)
#   AC_MSG_RESULT(yes)
##   AC_DEFINE(WITH_PYTHIA8)
#    WITH_PYTHIA8 = 1;
#   if test -d "$withval"; then
#   	PYTHIA8LIBS="-L$withval -lpythia8"
#   else
#   	AC_MSG_RESULT(can not find Hijing library)
#   fi
#   ;;
# esac ],
# AC_MSG_RESULT(no)
## AC_DEFINE(WITHOUT_PYTHIA8)
#  WITHOUT_PYTHIA8 = 1;
#)

#######################################################################
#AC_MSG_CHECKING(whether to use Hijing)
#AC_ARG_ENABLE(hijing,
#[  --enable-hijing       compile with Hijing],
#[case "${enableval}" in
#  yes)  hijing=true
#        AC_MSG_RESULT(yes) ;;
#  no)   hijing=false
#        AC_MSG_RESULT(no) ;;
#  *) AC_MSG_ERROR(bad value ${enableval} for --enable-hijing) ;;
#esac],
#[
#        hijing=false
#        AC_MSG_RESULT(no)
#])
#AM_CONDITIONAL(PYTHIA8BB, test x$hijing = xtrue)
#######################################################################

#if test -n "$PYTHIA8LIBS"; then
#	FORTLIBS="-L/usr/lib/gcc-lib/i386-redhat-linux/2.96/ -lg2c"
#fi



AC_CONFIG_FILES([Makefile jam2/Makefile
                 jam2/hadrons/Makefile
                 jam2/xsection/Makefile
                 jam2/initcond/Makefile
                 jam2/collision/Makefile
                 jam2/interaction/Makefile
                 jam2/fluid/Makefile
                 jam2/meanfield/Makefile
                 ])
AC_CONFIG_SUBDIRS([tmp])
AC_OUTPUT

