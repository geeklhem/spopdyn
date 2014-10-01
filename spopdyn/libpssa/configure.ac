
#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(libpssa,1.0.0RC, incardon@mpi-cbg.de)
LT_INIT()
AC_CANONICAL_SYSTEM
AC_CONFIG_SRCDIR([libpssa/src/PSSA.cpp])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADER([libpssa/src/config.h])
m4_ifdef([AX_CHECK_COMPILE_FLAG],,[m4_include([m4/ax_check_compile_flag.m4])])
m4_ifdef([AX_GCC_ARCHFLAG],,[m4_include([m4/ax_gcc_archflag.m4])])
m4_ifdef([AX_GCC_X86_CPUID],,[m4_include([m4/ax_gcc_x86_cpuid.m4])])
m4_ifdef([AX_BOOST_BASE],,[m4_include([m4/ax_boost_base.m4])])
m4_ifdef([AX_BOOST_PROGRAM_OPTIONS],,[m4_include([m4/ax_boost_program_options.m4])])
m4_ifdef([AX_PATH_GSL],,[m4_include([m4/gsl.m4])])
m4_ifdef([ACX_MPI],,[m4_include([m4/acx_mpi.m4])])

CXXFLAGS="$CXXFLAGS"

# Check for args


###### Check for debug compilation

AC_MSG_CHECKING(whether to build with debug information)
debuger=no
AC_ARG_ENABLE(debug,
	AC_HELP_STRING(
	    [--enable-debug],
	    [enable debug data generation (def=no)]
	),
  	debuger="$enableval"
)


AC_MSG_RESULT($debuger)
if test x"$debuger" = x"yes"; then
	AC_DEFINE([DEBUG_MODE],[],[Debug])
	AC_DEFINE([DEBUG],[],[Debug])
	CXXFLAGS="$CXXFLAGS -g3  -Wall -O0 "
else
	CXXFLAGS="$CXXFLAGS -O2 "
fi

# Check for MPI

AC_ARG_WITH([mpi],
        AS_HELP_STRING([--with-mpi], [build with MPI]),[
        AC_LANG_PUSH([C++])
        ACX_MPI([
                use_mpi=yes
                 ],[AC_MSG_ERROR([No MPI found])])
        AC_LANG_POP([C++])
        AC_DEFINE([HAVE_MPI],[],[Use Mpi])

],[use_mpi=no])

# check for header

AC_PROG_CXX
AC_LANG_PUSH([C++])
AC_CHECK_HEADER([time.h],,[AC_MSG_ERROR([*** This package needs cmath header])])
AC_CHECK_HEADER([cmath],,[AC_MSG_ERROR([*** This package needs cmath header])])
AC_CHECK_HEADER([fstream],,[AC_MSG_ERROR([*** This package needs set header])])
AC_CHECK_HEADER([iostream],,[AC_MSG_ERROR([*** This package needs set header])])
AC_CHECK_HEADER([sys/stat.h],,[AC_MSG_ERROR([*** This package needs set header])])
AC_CHECK_HEADER([vector],,[AC_MSG_ERROR([*** This package needs set header])])
AC_CHECK_HEADER([string.h],,[AC_MSG_ERROR([*** This package needs set header])])
AC_CHECK_HEADER([set],,[AC_MSG_ERROR([*** This package needs set header])])
AC_CHECK_HEADER([stdlib.h],,[AC_MSG_ERROR([*** This package needs stdlib.h])])
AC_LANG_POP([C++])

# Checks for libraries.

AX_PATH_GSL("1.15",
                [],
                [AC_MSG_ERROR("Could not find GSL package, use the following command, Red Hat, Fedora, CentOS: "yum install gsl-devel" or Ubuntu: "apt-get install gsl-dev", in order to solve the problem")],[],[])

AX_BOOST_BASE([1.33],, AC_MSG_ERROR([*** Boost library not found. Please install it.]))
AX_BOOST_PROGRAM_OPTIONS

# Checks g++ flags

AC_CANONICAL_HOST

if test x"$is_custom_arch" = x"no"; then 
  case $build_cpu in
  i[[3456]]86*) CXXFLAGS="$CXXFLAGS"
                no_64=yes;;
  x86_64*) CXXFLAGS="$CXXFLAGS -m64"
                   no_64=no;;
  powerpc64*) CXXFLAGS="$CXXFLAGS"
                     no_64=no;;
  powerpc*) CXXFLAGS="$CXXFLAGS"
            LDFLAGS="$LDFLAGS  -qnostaticlink -qnostaticlink=libgcc"
                     use_bgp=yes
                     no_64=yes;;
  *) AC_MSG_ERROR([*** $build_cpu ,Unsupported build cpu])
  esac
else
  no_64=yes
fi

# Check target architetture

is_custom_arch=no
if test x"$is_custom_arch" = x"no"; then
    AX_GCC_ARCHFLAG([], [CXXFLAGS="$CXXFLAGS"], [])
else
    ax_cv_gcc_archflag="$custom_arch"
fi

no_avx=no
no_sse42=no
no_sse41=no
no_sse3=no
no_sse2=no
no_sse=no
no_mmx=no

if test x"$no_64" = x"no"; then
	AX_CHECK_COMPILE_FLAG([-m64],[CXXFLAGS="$CXXFLAGS -m64"],[no_64=yes])
fi

case "${host}" in
 *-*-darwin*)
  AC_MSG_CHECKING([for a avx processor])
  if [ test `sysctl -n hw.optional.avx1_0` == 1]
  then
	AC_MSG_RESULT([yes])
	no_avx=no
  else
	AC_MSG_RESULT([no])
	no_avx=yes
  fi

AC_MSG_CHECKING([for a sse4_2 processor])
if [ test `sysctl -n hw.optional.sse4_2` == 1]
then
	AC_MSG_RESULT([yes])
	no_sse42=no
else
	AC_MSG_RESULT([no])
	no_sse42=yes
fi

AC_MSG_CHECKING([for a sse4_1 processor])
if [ test `sysctl -n hw.optional.sse4_1` == 1]
then
	AC_MSG_RESULT([yes])
	no_sse41=no
else
	AC_MSG_RESULT([no])
	no_sse41=yes
fi

AC_MSG_CHECKING([for a sse3 processor])
if [ test `sysctl -n hw.optional.sse3` == 1]
then
	AC_MSG_RESULT([yes])
	no_sse3=no
else
	AC_MSG_RESULT([no])
	no_sse3=yes
fi

AC_MSG_CHECKING([for a sse2 processor])
if [ test `sysctl -n hw.optional.sse2` == 1]
then
	AC_MSG_RESULT([yes])
	no_sse2=no
else
	AC_MSG_RESULT([no])
	no_sse2=yes
fi

AC_MSG_CHECKING([for a sse processor])
if [ test `sysctl -n hw.optional.sse` == 1]
then
	AC_MSG_RESULT([yes])
	no_sse=no
else
	AC_MSG_RESULT([no])
	no_sse=yes
fi
;;
*)
if grep -q "avx" /proc/cpuinfo
then
	AC_MSG_RESULT([yes])
	no_avx=no
else
	AC_MSG_RESULT([no])
	no_avx=yes
fi

AC_MSG_CHECKING([for a sse4_2 processor])
if grep -q "sse4_2" /proc/cpuinfo
then
	AC_MSG_RESULT([yes])
	no_sse42=no
else
	AC_MSG_RESULT([no])
	no_sse42=yes
fi

AC_MSG_CHECKING([for a sse4_1 processor])
if grep -q "sse4_1" /proc/cpuinfo
then
	AC_MSG_RESULT([yes])
	no_sse41=no
else
	AC_MSG_RESULT([no])
	no_sse41=yes
fi

AC_MSG_CHECKING([for a sse3 processor])
if grep -q "sse3" /proc/cpuinfo
then
	AC_MSG_RESULT([yes])
	no_sse3=no
else
	AC_MSG_RESULT([no])
	no_sse3=yes
fi

AC_MSG_CHECKING([for a sse2 processor])
if grep -q "sse2" /proc/cpuinfo
then
	AC_MSG_RESULT([yes])
	no_sse2=no
else
	AC_MSG_RESULT([no])
	no_sse2=yes
fi

AC_MSG_CHECKING([for a sse processor])
if grep -q "sse" /proc/cpuinfo
then
	AC_MSG_RESULT([yes])
	no_sse=no
else
	AC_MSG_RESULT([no])
	no_sse=yes
fi
;;
esac

if test x"$no_avx" = x"no"; then
     AX_CHECK_COMPILE_FLAG([-mavx],[CXXFLAGS="$CXXFLAGS -mavx"],[])
fi

if test x"$no_sse42" = x"no"; then
     AX_CHECK_COMPILE_FLAG([-msse4.2],[CXXFLAGS="$CXXFLAGS -msse4.2"],[])
fi

if test x"$no_avx" = x"no"; then
     AX_CHECK_COMPILE_FLAG([-msse4.1],[CXXFLAGS="$CXXFLAGS -msse4.1"],[])
fi

if test x"$no_avx" = x"no"; then
     AX_CHECK_COMPILE_FLAG([-msse3],[CXXFLAGS="$CXXFLAGS -msse3"],[])
fi

if test x"$no_avx" = x"no"; then
     AX_CHECK_COMPILE_FLAG([-msse2],[CXXFLAGS="$CXXFLAGS -msse2"],[])
fi

if test x"$no_avx" = x"no"; then
     AX_CHECK_COMPILE_FLAG([-msse],[CXXFLAGS="$CXXFLAGS -msse"],[])
fi

AX_CHECK_COMPILE_FLAG([-m64],[CXXFLAGS="$CXXFLAGS -m64"],[])
AX_CHECK_COMPILE_FLAG([-mmmx],[CXXFLAGS="$CXXFLAGS -mmmx"],[no_mmx=yes])
AX_CHECK_COMPILE_FLAG([-Wno-unused-but-set-variable],[CXXFLAGS="$CXXFLAGS -Wno-unused-but-set-variable"],[])

#AX_CHECK_COMPILE_FLAG([-frounding-math],[CXXFLAGS="$CXXFLAGS -frounding-math"],[])

# Define avx or sse

if test x"$no_avx" = x"no"; then
	AC_DEFINE([HAVE_AVX],[],[Enable avx])
elif test x"$no_sse41" = x"no"; then
	AC_DEFINE([HAVE_SSE],[],[Enable sse])
fi

# Checks for typedefs, structures, and compiler characteristics.

# Checks for programs.
if [ test x"$use_mpi" = x"yes" ]; then
  CXX="$MPICXX"
fi

# Checks for library functions.

AC_CHECK_LIB([sbml],[zipOpen],[],AC_MSG_ERROR("Install libSMBL"),[])
case "${host_os}" in
  darwin*)
    # Nothing to check
    ;;
  *)
    AC_CHECK_LIB([rt],[clock_gettime],[],AC_MSG_ERROR("Install librt"),[])
    OPT_LIBS="$OPT_LIBS -lrt"
    ;;
esac

AC_CONFIG_FILES([Makefile libpssa/Makefile libpssa/src/Makefile pssa_cli/Makefile pssa_test/Makefile examples/Makefile examples/Gray_Scott_2D/Makefile])
AC_OUTPUT
echo ""
echo   "***********************************"
echo   "*                                 *"
#arch_str="${ax_cv_gcc_archflag#-march=#-mtune=}"
#arch_str="${arch_str#-mtune=}"
#n_arch_str=${#arch_str}
#for (( X=0; X<23-n_arch_str; X++ ))
#do
#   arch_str="$arch_str "
#done
#echo "*    arch: $arch_str*"
if [ test x"$no_avx" = x"no" ]; then
    echo "*    avx   : yes                  *"
else
    echo "*    avx   : no                   *"
fi
if [ test x"$no_sse42" = x"no" ]; then
    echo "*    sse4.2: yes                  *"
else
    echo "*    sse4.2: no                   *"
fi
if [ test x"$no_sse41" = x"no" ]; then
    echo "*    sse4.1: yes                  *"
else
    echo "*    sse4.1: no                   *"
fi
if [ test x"$no_sse3" = x"no" ]; then
    echo "*    sse3: yes                    *"
else
    echo "*    sse3: no                     *"
fi
if [ test x"$no_sse2" = x"no" ]; then
    echo "*    sse2: yes                    *"
else
    echo "*    sse2: no                     *"
fi
if [ test x"$no_sse" = x"no" ]; then
    echo "*    sse: yes                     *"
else
    echo "*    sse: no                      *"
fi
if [ test x"$no_mmx" = x"no" ]; then
  echo "*    mmx: yes                     *"
else
  echo "*    mmx: no                      *"
fi
if [ test x"$use_mpi" = x"no" ]; then
  echo "*    MPI: no                      *"
else
  echo "*    MPI: yes                     *"
fi
if [ test x"$debuger" = x"yes" ]; then
  echo "*    debug: yes                   *"
else
  echo "*    debug: no                    *"
fi
if [ test x"$no_64" = x"no" ]; then
  echo "*    64 bit: yes                  *"
else
  echo "*    64 bit: no                   *"
fi
echo   "*                                 *"
echo   "***********************************"
echo   ""
echo   "In order to compile:"
echo   "make"
echo   ""
echo   "Getting started documentation: doc/pSSALib.hml"
echo   ""