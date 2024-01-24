# Scons "makefile" for imfit and related programs
# 
# To use this, type the following on the command line:
#    $ scons <target-name>
# where <target-name> is, e.g., "imft", "makeimage", etc.
# "scons" by itself will build *all* targets
# 
# To run unit tests (if you have the full Github-based distribution)
#    $ scons unit
#
# To clean up files associated with a given target:
#    $ scons -c <target-name>
# To clean up all targets (including programs):
#    $ scons -c


# *** SPECIAL STUFF ***
# To build a version with extra, experimental functions
#    $ scons --extra-funcs <target-name>
#
# To build a version *without* OpenMP enabled
#    $ scons --no-openmp <target-name>
#
# To build a version *without* FFTW threading:
#    $ scons --no-threading <target-name>
#
#
# To build a version using non-default compiler:
#    $ scons --cc=<C_COMPILER> --cpp=<C++_COMPILE> <target-name>
# e.g.
#    $ scons --cc=gcc-4.9 --cpp=g++-4.9 <target-name>
# shorthand for using GCC 9
#    $ scons --use-gcc <target-name>
#
#
# To build a version with full debugging printouts:
#    $ scons define=DEBUG <target-name>
#
# To build export version (all non-system libraries statically linked):
#    $ scons --static <target-name>
#

# *** EXPORT CONFIGURATIONS ***
# macOS static binaries (including static linking of libgcc and libstdc++)
# $ scons --allstatic --mac-distribution

# To add one or more directories to the header or library search paths:
#    $ scons --header-path=/path/to/header/dir
# OR $ scons --header-path=/path/to/header/dir:/alt/path:/another/path
#    $ scons --lib-path=/path/to/lib/dir
# etc.


# Copyright 2010--2020 by Peter Erwin.
# 
# This file is part of Imfit.
# 
# Imfit is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your
# option) any later version.
# 
# Imfit is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
# 
# You should have received a copy of the GNU General Public License along
# with Imfit.  If not, see <http://www.gnu.org/licenses/>.


import os, subprocess, platform, getpass, pwd


def GetLinuxType():
    txt = subprocess.check_output("hostnamectl")
    txt = txt.decode('utf-8')
    lines = txt.splitlines()

    linux_name = None
    for line in lines:
        pp = line.split(":")
        if pp[0].strip() == "Operating System":
            linux_name = pp[1].split()[0]
    return linux_name

# possible values: "Darwin", "Linux"
os_type = platform.system()
linux_type = None
# possible values: "Ubuntu", "CentOS", ...
if os_type == "Linux":
    linux_type = GetLinuxType()



# LIBRARIES:
# m
# pthread [Linux]
# dl [Linux, if using loguru for logging]
# cfitsio
#   -- if static, then on macOS we must link with curl [part of system]
#   -- for Linux, we use our compiled static-library version of cfitsio
#      (not Ubuntu's), so we don't need any extra libraries
# fftw3, fftw3_threads
# [gsl, gslcblas]
# [nlopt]
#   -- note that static-library version of this (libnlopt.a) should be the 
#      *non*-C++ version (build with cmake ... -DNLOPT_CXX=OFF ...),
#      otherwise there can be problems when compiler for these programs is
#      different from compiler used to build libnlopt.a.


# *** Hard-coded paths to static libraries (to ensure static linking when that's what
# we want to do; compilers sometimes insist on linking to shared library even when you 
# don't specify that)

# NOTE TO USERS: For normal compilation, don't use static linking, since it's
# harder to set up and ensure the right library paths. The static-library stuff
# here is primarily for use in preparing the binary distributions of Imfit.

# We assume that FFTW library is static-only (since that's the default installation).

MAC_STATIC_LIBS_PATH = "/usr/local/lib/"
# Debian/Ubuntu standard x86-64 package installation path
LINUX_UBUNTU_STATIC_LIBS_PATH = "/usr/local/lib/"
libDirs = {"Darwin": MAC_STATIC_LIBS_PATH, "Linux": LINUX_UBUNTU_STATIC_LIBS_PATH}
extraSharedLibs_static_cfitsio = {"Darwin": ["curl"], "Linux": []}

BASE_SHARED_LIBS = ["m"]
if os_type == "Linux":
    BASE_SHARED_LIBS.append("pthread")


# The following definitions are for when we want to force static linking to
# the cfitsio, FFTW3, GSL, and NLopt libraries.
# (Change these if the locations are different on your system)

# cfitsio weirdness: if we link to the static-library version, then we
# must also link to the (shared library) libcurl on Mac;
# for Linux, we link to our pre-compiled static library in /usr/local/lib
# Also, static-library path is different than for other static libraries on Linux
# VM installation.
STATIC_CFITSIO_LIBRARY_FILE = File(libDirs[os_type] + "libcfitsio.a")
BASE_SHARED_LIBS += extraSharedLibs_static_cfitsio[os_type]

STATIC_FFTW_LIBRARY_FILE = File(libDirs[os_type] + "libfftw3.a")
STATIC_FFTW_THREADED_LIBRARY_FILE = File(libDirs[os_type] + "libfftw3_threads.a")

# The following is for when we want to force static linking to the GSL library
# (Change these if the locations are different on your system)
STATIC_GSL_LIBRARY_FILE1 = File(libDirs[os_type] + "libgsl.a")
STATIC_GSL_LIBRARY_FILE2 = File(libDirs[os_type] + "libgslcblas.a")

# the following is for when we want to force static linking to the NLopt library
# (Change these if the locations are different on your system)
STATIC_NLOPT_LIBRARY_FILE = File(libDirs[os_type] + "libnlopt.a")
STATIC_NLOPT_LIBRARY_FILE_MACOSX_NOTHREADLOCAL = File("/Users/erwin/coding/imfit/static_libs/nlopt_nothreadlocal/libnlopt.a")



# locations of source-code files (including our header files)
CORE_SUBDIR = "core/"
FUNCTION_SUBDIR = "function_objects/"
FUNCTION_1D_SUBDIR = "function_objects_1d/"
SOLVER_SUBDIR = "solvers/"
CDREAM_SUBDIR = "cdream/"
CDREAM_INCLUDE_SUBDIR = CDREAM_SUBDIR + "include/"
PROFILEFIT_SUBDIR = "profile_fitting/"

BAD_OPENMP_COMPILER_WARNING = """
*** WARNING: You appear to be trying to compile with OpenMP support under macOS 
while using Apple's (Xcode) clang++, which does NOT support OpenMP.
(Note that clang++ may be aliased on a macOS system as /usr/bin/g++.)

Either go ahead and compile without OpenMP by calling scons with the --no-openmp 
flag, or (ideally) use a different compiler (e.g., GCC, or a version of Clang/LLVM
with OpenMP support).
"""


# *** Set up compiler flags, library lists, include paths

# Note on sse2, etc.:
#    SSE2 is supported on Intel processors from xxx onward
#    AVX  is supported on Intel "Core i3/i5/i7" processors from 2011 onward
#    AVX2  is supported on Intel Haswell and later processors (mostly 2014 onward)
#    AVX-512  is supported only on "Knights Landing" Xeon Phi processors (2106 onward)

cflags_opt = ["-O3", "-g0", "-fPIC", "-msse2", "-std=c++11"]
cflags_db = ["-Wall", "-g3", "-O0", "-fPIC", "-std=c++11", "-Wshadow", 
                "-Wredundant-decls", "-Wpointer-arith"]

base_defines = ["ANSI", "USING_SCONS"]

# libraries needed for imfit, makeimage, psfconvolve, & other 2D programs
lib_list = BASE_SHARED_LIBS
# libraries needed for profilefit and psfconvolve1d compilation
lib_list_1d = ["fftw3", "fftw3_threads", "m"]


include_path = [".", "/usr/local/include", CORE_SUBDIR, SOLVER_SUBDIR, CDREAM_SUBDIR,
                CDREAM_INCLUDE_SUBDIR, FUNCTION_SUBDIR, FUNCTION_1D_SUBDIR, PROFILEFIT_SUBDIR]
lib_path = ["/usr/local/lib"]
link_flags = []


# find out what the default compilers (according to SCons) are
default_environ = DefaultEnvironment()
cc_default = default_environ["CC"]
cpp_default = default_environ["CXX"]
CC_COMPILER = cc_default
CPP_COMPILER = cpp_default
c_compiler_changed = False
cpp_compiler_changed = False
usingGCC = True


# ** Special setup for compilation by P.E. on Mac (assumes GCC v9 is installed and
# callable via gcc-9 and g++-9)
# Comment this out otherwise!
# Note that the following way of determining the username seems to be a bit more 
# portable than "getpass.getuser()", which fails for "Ubuntu on Windows" (acc. to 
# Lee Kelvin, who contributed the new version)
userName = pwd.getpwuid(os.getuid())[0]
if (os_type == "Darwin") and (userName == "erwin"): 
    CC_COMPILER = "gcc-9"
    CPP_COMPILER = "g++-9"
    c_compiler_changed = True
    cpp_compiler_changed = True
    usingGCC = True


# *** System-specific setup
# if (os_type == "Darwin"):   # OK, we're compiling on macOS (a.k.a. Mac OS X)
#   # Note: if for some reason you need to compile to 32-bit -- e.g., because
#   # your machine is 32-bit only, or because the fftw3 and cfitsio libraries
#   # are 32-bit, use the following
#   cflags_db = ["-Wall", "-Wshadow", "-Wredundant-decls", "-Wpointer-arith", "-g3"]
if (os_type == "Linux"):
    # change the following path definitions as needed
    include_path.append("/usr/include")
    if userName == "erwin":
        include_path.append("/home/erwin/include")
        lib_path.append("/home/erwin/lib")
defines_opt = base_defines
defines_db = base_defines

extra_defines = []



# *** Set up default settings, check for user-requested changes

# Default settings for compilation
#useGSL = True
useNLopt = True
useOpenMP = True
usingClangOpenMP = False
useExtraFuncs = False
useStaticLibs = False
totalStaticLinking = False
buildFatBinary = False
build32bit = False
buildForOldMacOS = False
scanBuild = False
addressSanitize = False
allSanitize = False
setOptToDebug = False
useLogging = False

# Define some user options
AddOption("--lib-path", dest="libraryPath", type="string", action="store", default=None,
    help="colon-separated list of additional paths to search for libraries")
AddOption("--header-path", dest="headerPath", type="string", action="store", default=None,
    help="colon-separated list of additional paths to search for header files")
AddOption("--no-threading", dest="fftwThreading", action="store_false", 
    default=True, help="compile programs *without* FFTW threading")
# AddOption("--no-gsl", dest="useGSL", action="store_false", 
#     default=True, help="do *not* use GNU Scientific Library")
AddOption("--no-nlopt", dest="useNLopt", action="store_false", 
    default=True, help="do *not* use NLopt library")
AddOption("--no-openmp", dest="noOpenMP", action="store_true", 
    default=False, help="compile *without* OpenMP support")
AddOption("--clang-openmp", dest="useClangOpenMP", action="store_true", 
    default=False, help="compile with clang++ on macOS")
AddOption("--extra-funcs", dest="useExtraFuncs", action="store_true", 
    default=False, help="compile additional FunctionObject classes for testing")
AddOption("--extra-checks", dest="doExtraChecks", action="store_true", 
    default=False, help="turn on additional error-checking and warning flags during compilation")
# options to specify use of non-default compilers, extra checks, loggin
AddOption("--cc", dest="cc_compiler", type="string", action="store", default=None,
    help="C compiler to use instead of system default")
AddOption("--cpp", dest="cpp_compiler", type="string", action="store", default=None,
    help="C++ compiler to use instead of system default")
AddOption("--use-gcc", dest="useGCC", action="store_true", 
    default=False, help="use gcc and g++ v9 compilers")
AddOption("--scan-build", dest="doingScanBuild", action="store_true", 
    default=False, help="set this when using scan-build (only for imfit_db and makeimage_db)")
AddOption("--sanitize", dest="useAllSanitize", action="store_true", 
    default=False, help="set this to generate binaries with -fsanitize-address, -fsanitize-undefined, and -fsanitize=leak")
AddOption("--address-sanitize", dest="useAddressSanitize", action="store_true", 
    default=False, help="set this to generate binaries with -fsanitize-address")
AddOption("--logging", dest="useLogging", action="store_true", 
    default=False, help="compile with support for logging via loguru")

# Define some more arcane options (e.g., for making binaries for distribution)
AddOption("--static", dest="useStaticLibs", action="store_true", 
    default=False, help="force static library linking")
AddOption("--allstatic", dest="useTotalStaticLinking", action="store_true", 
    default=False, help="force static library linking, *including* system libraries if possible")
AddOption("--mac-distribution", dest="compileForMacDistribution", action="store_true", 
    default=False, help="use this to make macOS binaries for public distribution")



# * Check to see if user actually specified something, and implement it
# special stuff for compiling Mac binary distribution
if GetOption("compileForMacDistribution"):
    print("DOING MAC DISTRIBUTION COMPILE!")
    usingClangOpenMP = True
    usingGCC = False
    CC_COMPILER = "clang"
    CPP_COMPILER = "clang++"
    useStaticLibs = True
#   totalStaticLinking = True
    # reset lib_path so linker only looks where we want it to -- i.e., in the 
    # directory with static library files
    lib_path = ["/Users/erwin/coding/imfit/static_libs"]

if GetOption("headerPath") is not None:
    extraPaths = GetOption("headerPath").split(":")
    print("extra header search paths: ", extraPaths)
    include_path += extraPaths
if GetOption("libraryPath") is not None:
    extraPaths = GetOption("libraryPath").split(":")
    print("extra library search paths: ", extraPaths)
    lib_path += extraPaths
# if GetOption("useGSL") is False:
#     useGSL = False
if GetOption("useNLopt") is False:
    useNLopt = False
if GetOption("noOpenMP"):
    useOpenMP = False
if GetOption("useClangOpenMP"):
    usingClangOpenMP = True
    usingGCC = False
    CC_COMPILER = "clang"
    CPP_COMPILER = "clang++"
if GetOption("useExtraFuncs"):
    useExtraFuncs = True
doExtraChecks = False
if GetOption("doExtraChecks"):
    doExtraChecks = True

# change the compilers if user requests it
if GetOption("cc_compiler") is not None:
    CC_COMPILER = GetOption("cc_compiler")
    print("using %s for C compiler" % CC_COMPILER)
    c_compiler_changed = True
if GetOption("cpp_compiler") is not None:
    CPP_COMPILER = GetOption("cpp_compiler")
    print("using %s for C++ compiler" % CPP_COMPILER)
    cpp_compiler_changed = True
if GetOption("useGCC"):
    if usingClangOpenMP:
        print("ERROR: You cannot specify both Clang and GCC as the compiler!")
        Exit(2)
    usingGCC = True
    CC_COMPILER = "gcc-9"
    CPP_COMPILER = "g++-9"
    print("using %s for C compiler" % CC_COMPILER)
    print("using %s for C++ compiler" % CPP_COMPILER)
    c_compiler_changed = True
    cpp_compiler_changed = True

if GetOption("doingScanBuild"):
    scanBuild = True
    useOpenMP = False   # scan-build uses clang, which doesn't have OpenMP

if GetOption("useAddressSanitize"):
    addressSanitize = True
    useOpenMP = False
    setOptToDebug = True
if GetOption("useAllSanitize"):
    allSanitize = True
    useOpenMP = False
    setOptToDebug = True
if GetOption("useLogging"):
    useLogging = True

if GetOption("useStaticLibs"):
    useStaticLibs = True
if GetOption("useTotalStaticLinking"):
    useStaticLibs = True
    if usingGCC:
        totalStaticLinking = True



# *** Setup for various options (either default, or user-altered)

if setOptToDebug:
    print("** Turning off optimizations!")
    cflags_opt = cflags_db


if useStaticLibs:
    lib_list.append(STATIC_CFITSIO_LIBRARY_FILE)
    lib_list.append(STATIC_FFTW_LIBRARY_FILE)
    lib_list.append(STATIC_FFTW_THREADED_LIBRARY_FILE)
else:
    lib_list += ["cfitsio", "fftw3", "fftw3_threads"]
extra_defines.append("FFTW_THREADING")

#if useGSL:   # true by default
# GNU Scientific Library
if useStaticLibs:
    lib_list.append(STATIC_GSL_LIBRARY_FILE1)
    lib_list.append(STATIC_GSL_LIBRARY_FILE2)
    lib_list_1d.append(STATIC_GSL_LIBRARY_FILE1)
    lib_list_1d.append(STATIC_GSL_LIBRARY_FILE2)
else:
    lib_list.append("gsl")
    lib_list.append("gslcblas")     
    lib_list_1d.append("gsl")
    lib_list_1d.append("gslcblas")      
# else:
#     extra_defines.append("NO_GSL")

if useNLopt:   # default is to do this
    if useStaticLibs:
        lib_list.append(STATIC_NLOPT_LIBRARY_FILE)
        lib_list_1d.append(STATIC_NLOPT_LIBRARY_FILE)
    else:
        if (os_type == "Linux") and (linux_type == "CentOS"):
            # CentOS package versions of NLopt tend to call the library file
            # "libnlopt_cxx" instead of just "libnlopt"
            lib_list.append("nlopt_cxx")    
            lib_list_1d.append("nlopt_cxx") 
        else:
            lib_list.append("nlopt")    
            lib_list_1d.append("nlopt") 
else:
    extra_defines.append("NO_NLOPT")



# Special case where we try to link libstdc++, libgomp, and libgcc_s statically
# (probably *won't* work with clang/clang++, only useful with GCC)
# Note that "-static-libstdc++" doesn't seem to accomplish anything; only
# "-static-libgcc" does (and that statically links in libstdc++ and libgomp
# in addition to libgcc_s). But we'll leave it in because it seems to be
# the standard.
if totalStaticLinking:
    link_flags.append("-static-libgcc")
    link_flags.append("-static-libstdc++")

if useOpenMP:   # default is to do this (turn this off with "--no-openmp")
    if usingClangOpenMP:
        # special flags for Apple clang++ (assumes libomp is installed)
        cflags_opt.append("-Xpreprocessor")
        cflags_opt.append("-fopenmp")
        cflags_db.append("-Xpreprocessor")
        cflags_db.append("-fopenmp")
        link_flags.append("-Xpreprocessor")
        link_flags.append("-fopenmp")
        if useStaticLibs:
            link_flags.append("/Users/erwin/coding/imfit/static_libs/libomp.a")
        else:
            # dynamic linking to libomp
            link_flags.append("-lomp")
    else:
        # flags for (real) g++
        cflags_opt.append("-fopenmp")
        cflags_db.append("-fopenmp")
        link_flags.append("-fopenmp")
    extra_defines.append("USE_OPENMP")

if useExtraFuncs:   # default is to NOT do this; user must specify with "--extra-funcs"
    extra_defines.append("USE_EXTRA_FUNCS")

if doExtraChecks:   # default is to NOT do this; user must specify with "--extra-checks"
    cflags_opt.append(["-Wall", "-Wshadow", "-Wredundant-decls", "-Wpointer-arith",
                    "-Wextra", "-pedantic"])

if useLogging:
    extra_defines.append(["-DUSE_LOGGING"])
    if os_type == "Linux":
        lib_list.append("dl")
    
# Add any additional, user-specified preprocessor definitions (e.g., "define=DEBUG")
for key, value in ARGLIST:
    if key == 'define':
        extra_defines.append(value)


if addressSanitize:
    cflags_opt.append("-fsanitize=address")
    cflags_opt.append("-fno-omit-frame-pointer")
    cflags_db.append("-fsanitize=address")
    cflags_db.append("-fno-omit-frame-pointer")
    link_flags.append("-fsanitize=address")
    link_flags.append("-fno-omit-frame-pointer")

if allSanitize:
    cflags_opt += ["-fsanitize=address", "-fsanitize=undefined", "-fsanitize=leak"]
    cflags_opt.append("-fno-omit-frame-pointer")
    cflags_db += ["-fsanitize=address", "-fsanitize=undefined", "-fsanitize=leak"]
    cflags_db.append("-fno-omit-frame-pointer")
    link_flags += ["-fsanitize=address", "-fsanitize=undefined", "-fsanitize=leak"]
    link_flags.append("-fno-omit-frame-pointer")


# * Collect together all the updated preprocessor definitions
defines_db = defines_db + extra_defines
defines_opt = defines_opt + extra_defines




# *** Create Environments for compilation:
# "env_debug" is environment with debugging options turned on
# "env" is an environment for optimized compiling

env = Environment( CC=CC_COMPILER, CXX=CPP_COMPILER, CPPPATH=include_path, LIBS=lib_list, 
                    LIBPATH=lib_path, CCFLAGS=cflags_opt, LINKFLAGS=link_flags, 
                    CPPDEFINES=defines_opt )
env_debug = Environment( CC=CC_COMPILER, CXX=CPP_COMPILER, CPPPATH=include_path, LIBS=lib_list, 
                    LIBPATH=lib_path, CCFLAGS=cflags_db, LINKFLAGS=link_flags, 
                    CPPDEFINES=defines_db )


# Checks for libraries and headers -- if we're not doing scons -c:
# WARNING: This is NOT a good idea for us at the moment, because
#    1. It fails to work on our Linux VM installation
#    2. It automatically inserts "-l<libname>" if it finds the libraries, which ends
#       up forcing the linking of dynamic-library versions even if we're trying to
#       do static compilation
# if not env.GetOption('clean'):
#   conf_opt = Configure(env)
#   cfitsioFound = conf_opt.CheckLibWithHeader('cfitsio', 'fitsio.h', 'c')
#   fftwFound = conf_opt.CheckLibWithHeader('fftw3', 'fftw3.h', 'c')
#   fftwThreadsFound = conf_opt.CheckLib('fftw3_threads')
#   nloptFound = conf_opt.CheckLibWithHeader('nlopt', 'nlopt.h', 'c')
#   gslFound = conf_opt.CheckLib('gsl')
#   libsOK = False
#   if cfitsioFound and fftwFound:
#       libsOK = True
#   else:
#       print("ERROR: Failed to find one or more required libraries and/or header files (cfitsio and/or fftw3)!")
#       print("\tMake sure they are installed; if necessary, include correct path to library with --lib-path option")
#       print("\tand correct path to header with --header-path option")
#       exit(1)
#   if useFFTWThreading and not fftwThreadsFound:
#       print("ERROR: Failed to find fftw3_threading library!")
#       print("\tSuggestion: include correct path to library with --lib-path option")
#       print("\tOR run SCons with --no-threading option")
#       exit(1)
#   if useGSL and not gslFound:
#       print("ERROR: Failed to find gsl library!")
#       print("\tSuggestion: include correct path to library with --lib-path option")
#       print("\tOR run SCons with --no-gsl option")
#       exit(1)
#   if useNLopt and not nloptFound:
#       print("ERROR: Failed to find nlopt library!")
#       print("\tSuggestion: include correct path to library with --lib-path option")
#       print("\tOR run SCons with --no-nlopt option")
#       exit(1)
#   env = conf_opt.Finish()



# We have separate lists of object names (what we want the .o files to be called) and
# source names (.cpp) so that we can specify separate debugging and optimized compilations.

# ModelObject and related classes/files:
modelobject_obj_string = """model_object convolver oversampled_region downsample
        psf_oversampling_info setup_model_object"""
modelobject_objs = [ CORE_SUBDIR + name for name in modelobject_obj_string.split() ]
modelobject_sources = [name + ".cpp" for name in modelobject_objs]

# Function objects:
functionobject_obj_string = """function_object func_gaussian func_exp func_gen-exp  
        func_sersic func_gen-sersic func_core-sersic func_broken-exp
        func_broken-exp2d func_moffat func_flatsky func_tilted-sky-plane 
        func_flatbar func_gaussian-ring func_spiral func_spiral_broken func_spiral_0b func_spiral_1b func_spiral_2b
        func_gaussian-ring2side func_gaussian-ring-az func_edge-on-disk_n4762 
        func_edge-on-disk_n4762v2 func_edge-on-ring func_edge-on-ring2side 
        func_king func_king2 func_ferrersbar2d 
        helper_funcs helper_funcs_3d psf_interpolators"""
#if useGSL:
# NOTE: the following modules require GSL be present
functionobject_obj_string += " func_edge-on-disk"
functionobject_obj_string += " integrator"
functionobject_obj_string += " func_expdisk3d"  # requires integrator
functionobject_obj_string += " func_brokenexpdisk3d"  # requires integrator
functionobject_obj_string += " func_gaussianring3d"  # requires integrator
functionobject_obj_string += " func_ferrersbar3d"  # requires integrator
functionobject_obj_string += " func_pointsource"
if useExtraFuncs:
    # experimental extra functions for personal testing
    functionobject_obj_string += " func_broken-exp-bar"
    functionobject_obj_string += " func_double-broken-exp"
    functionobject_obj_string += " func_gen-exp2"
#    functionobject_obj_string += " func_flatbar"
    functionobject_obj_string += " func_flatbar_trunc"
    functionobject_obj_string += " func_gen-flatbar"
    functionobject_obj_string += " func_bp-cross-section"
    functionobject_obj_string += " func_lorentzian-ring"
    functionobject_obj_string += " func_n4608disk"
#    if useGSL:
    functionobject_obj_string += " func_brokenexpbar3d"
    functionobject_obj_string += " func_boxytest3d"
    functionobject_obj_string += " func_double-brokenexpdisk3d"
    functionobject_obj_string += " func_expdisk3d_trunc"
    functionobject_obj_string += " func_logspiral"
    functionobject_obj_string += " func_logspiral2"
    functionobject_obj_string += " func_logspiral_gauss"
    functionobject_obj_string += " func_nan"
    functionobject_obj_string += " func_triaxbar3d"
    functionobject_obj_string += " func_triaxbar3d_sq"
    functionobject_obj_string += " func_triaxbar3d_gengauss_sq"
    functionobject_obj_string += " func_exp-higher-mom"
# ADD CODE FOR NEW FUNCTIONS HERE
# (NOTE: be sure to include one or more spaces before the file name!)
# e.g.,
# functionobject_obj_string += " func_mynewfunction"

functionobject_objs = [ FUNCTION_SUBDIR + name for name in functionobject_obj_string.split() ]
functionobject_sources = [name + ".cpp" for name in functionobject_objs]


# Solvers and associated code
solver_obj_string = """levmar_fit mpfit diff_evoln_fit DESolver dispatch_solver solver_results"""
if useNLopt:
    solver_obj_string += " nmsimplex_fit nlopt_fit"
solver_objs = [ SOLVER_SUBDIR + name for name in solver_obj_string.split() ]
solver_sources = [name + ".cpp" for name in solver_objs]


# CDREAM and associated code for MCMC
cdream_obj_string = """check_outliers dream dream_initialize dream_pars gelman_rubin gen_CR
restore_state"""
cdream_objs = [ CDREAM_SUBDIR + name for name in cdream_obj_string.split() ]
cdream_sources = [name + ".cpp" for name in cdream_objs]


# Base files for imfit, makeimage, imfit-mcmc, and libimfit:
base_obj_string = """mp_enorm statistics mersenne_twister commandline_parser utilities 
config_file_parser add_functions"""
base_objs = [ CORE_SUBDIR + name for name in base_obj_string.split() ]
# FITS image-file I/O
image_io_obj_string = "image_io getimages"
image_io_objs = [ CORE_SUBDIR + name for name in image_io_obj_string.split() ]

# Main set of files for imfit
imfit_obj_string = """print_results bootstrap_errors estimate_memory 
imfit_main"""
if useLogging:
    imfit_obj_string += " loguru/loguru"
imfit_base_objs = [ CORE_SUBDIR + name for name in imfit_obj_string.split() ]
imfit_base_objs = base_objs + image_io_objs + imfit_base_objs
imfit_base_sources = [name + ".cpp" for name in imfit_base_objs]

# Main set of files for makeimage
makeimage_base_objs = base_objs + image_io_objs + [CORE_SUBDIR + "makeimage_main"]
if useLogging:
    makeimage_base_objs.append("loguru/loguru")
makeimage_base_sources = [name + ".cpp" for name in makeimage_base_objs]

# Main set of files for imfit-mcmc
mcmc_obj_string = """estimate_memory mcmc_main"""
mcmc_base_objs = [ CORE_SUBDIR + name for name in mcmc_obj_string.split() ]
mcmc_base_objs = mcmc_base_objs + base_objs + image_io_objs + cdream_objs
mcmc_base_sources = [name + ".cpp" for name in mcmc_base_objs]


# imfit: put all the object and source-code lists together
imfit_objs = imfit_base_objs + modelobject_objs + functionobject_objs + solver_objs
imfit_sources = imfit_base_sources + modelobject_sources + functionobject_sources + solver_sources

# makeimage: put all the object and source-code lists together
makeimage_objs = makeimage_base_objs + modelobject_objs + functionobject_objs
makeimage_sources = makeimage_base_sources + modelobject_sources + functionobject_sources

# imfit-mcmc: put all the object and source-code lists together
mcmc_objs = mcmc_base_objs + modelobject_objs + functionobject_objs
mcmc_sources = mcmc_base_sources + modelobject_sources + functionobject_sources


# import environment variables if we're doing scan-build static analysis
if scanBuild:
    env_debug["CC"] = os.getenv("CC")
    env_debug["CXX"] = os.getenv("CXX")
    env_debug["ENV"].update(x for x in os.environ.items() if x[0].startswith("CCC_"))


# *** Finally, define the actual targets for building
# specify ".do" as the suffix for "full-debug" object code
imfit_dbg_objlist = [ env_debug.Object(obj + ".do", src) for (obj,src) in zip(imfit_objs, imfit_sources) ]
env_debug.Program("imfit_db", imfit_dbg_objlist)
imfit_opt_objlist = [ env.Object(obj, src) for (obj,src) in zip(imfit_objs, imfit_sources) ]
env.Program("imfit", imfit_opt_objlist)

makeimage_dbg_objlist = [ env_debug.Object(obj + ".do", src) for (obj,src) in zip(makeimage_objs, makeimage_sources) ]
env_debug.Program("makeimage_db", makeimage_dbg_objlist)
env.Program("makeimage", makeimage_sources)

mcmc_dbg_objlist = [ env_debug.Object(obj + ".do", src) for (obj,src) in zip(mcmc_objs, mcmc_sources) ]
env_debug.Program("imfit-mcmc_db", mcmc_dbg_objlist)
env.Program("imfit-mcmc", mcmc_sources)


# Run tests
# Unit tests:
env.Command("unit", None, "./run_unit_tests.sh")

# All tests:
env.Command("alltests", None, 
    "./run_unit_tests.sh ; ./do_makeimage_tests ; ./do_imfit_tests ; ./do_mcmc_tests")




# Build Imfit library
base_for_lib_objstring = """mp_enorm statistics mersenne_twister utilities 
config_file_parser add_functions bootstrap_errors"""
base_for_lib_objs = [ CORE_SUBDIR + name for name in base_for_lib_objstring.split() ]
libimfit_objs = modelobject_objs + functionobject_objs + solver_objs
libimfit_objs += base_for_lib_objs
libimfit_sourcelist = [name + ".cpp" for name in libimfit_objs]
#print(imfit_lib_sourcelist)
# Note that for some reason we have to give the library name with its
# "lib" prefix: "libimfit"
# invoke the following as "scons libimfit.a" for the static-library version
# and "scons libimfit.dylib" (macOS) or "scons libimfit.so" (Linux) for the
# dynamic-library version
libimfit_objlist = [ env.Object(obj, src) for (obj,src) in zip(libimfit_objs, libimfit_sourcelist) ]
staticlib = env.StaticLibrary(target="libimfit", source=libimfit_objlist)
# THE FOLLOWING CURRENTLY DOES NOT WORK
# ("Source file: core/model_object.o is static and is not compatible with shared target: libimfit.dylib")
# sharedlib = env.SharedLibrary(target="libimfit", source=libimfit_objlist,
#                                   LIBS=lib_list_libimfit)





# *** Other programs (profilefit, psfconvolve, older stuff)
# From here to the end of the file: removed from exported distribution version of SConstruct

env_1d = Environment( CC=CC_COMPILER, CXX=CPP_COMPILER, CPPPATH=include_path, LIBS=lib_list_1d, LIBPATH=lib_path,
                        CCFLAGS=cflags_db, LINKFLAGS=link_flags, CPPDEFINES=defines_db )

# ModelObject1d and related classes:
# (Note that model_object includes references to oversampled_region and downsample,
# so we need to include those in the compilation and link, even though they aren't
# actually used in model_object1d. Similarly, code in image_io is referenced from
# downsample.)
modelobject1d_obj_string = """model_object oversampled_region downsample psf_oversampling_info"""
modelobject1d_objs = [CORE_SUBDIR + name for name in modelobject1d_obj_string.split()]
modelobject1d_sources = [name + ".cpp" for name in modelobject1d_objs]

# 1D FunctionObject classes (note that we have to add a separate entry for function_object.cpp,
# which is in a different subdirectory):
functionobject1d_obj_string = """func1d_gaussian func1d_gaussian_linear func1d_exp func1d_sersic 
        func1d_core-sersic func1d_broken-exp func1d_moffat func1d_delta func1d_sech 
        func1d_sech2 func1d_vdksech func1d_gaussian2side  func1d_nuker func1d_spline
        func1d_n1543majmin_circbulge func1d_n1543majmin func1d_n1543majmin2
        func1d_double-gauss-hermite func1d_gauss-hermite"""
functionobject1d_objs = [ FUNCTION_1D_SUBDIR + name for name in functionobject1d_obj_string.split() ]
functionobject1d_objs.append(FUNCTION_SUBDIR + "function_object")
functionobject1d_sources = [name + ".cpp" for name in functionobject1d_objs]

# Base files for profilefit:
profilefit_base_obj_string = """core/commandline_parser core/utilities profile_fitting/read_profile 
        core/config_file_parser core/print_results profile_fitting/add_functions_1d core/convolver 
        core/mp_enorm core/statistics core/mersenne_twister 
        function_objects/psf_interpolators
        profile_fitting/convolver1d profile_fitting/model_object_1d 
        profile_fitting/bootstrap_errors_1d profile_fitting/profilefit_main"""
profilefit_base_objs = profilefit_base_obj_string.split()
profilefit_base_sources = [name + ".cpp" for name in profilefit_base_objs]

# profilefit: put all the object and source-code lists together
profilefit_objs = profilefit_base_objs + modelobject1d_objs + functionobject1d_objs + solver_objs
profilefit_sources = profilefit_base_sources + modelobject1d_sources + functionobject1d_sources + solver_sources

# psfconvolve1d: put all the object and source-code lists together
psfconvolve1d_objs = ["profile_fitting/psfconvolve1d_main", "core/commandline_parser", "core/utilities",
                    "profile_fitting/read_profile", "profile_fitting/convolver1d"]
psfconvolve1d_sources = [name + ".cpp" for name in psfconvolve1d_objs]



# source+obj lists for older or less-used programs:
# readimage: put all the object and source-code lists together
readimage_sources = ["readimage_main.cpp", "image_io.cpp"]

# psfconvolve: put all the object and source-code lists together
psfconvolve_objs = ["extra/psfconvolve_main", "core/commandline_parser", "core/utilities",
                    "core/image_io", "core/convolver"]
psfconvolve_sources = [name + ".cpp" for name in psfconvolve_objs]

# test_parser: put all the object and source-code lists together
testparser_objs = ["test_parser", "core/config_file_parser", "core/utilities"]
testparser_sources = [name + ".cpp" for name in testparser_objs]


# test_2dspline: put all the object and source-code lists together
spline2dtest_objs = ["spline2dtest_main", "function_objects/psf_interpolators", 
                    "core/commandline_parser", "core/utilities", "core/image_io", 
                    "function_objects/function_object", "function_objects/func_pointsource"]
spline2dtest_sources = [name + ".cpp" for name in spline2dtest_objs]





# profile fit is fast, so we don't really need an "optimized" version
profilefit_dbg_objlist = [ env_debug.Object(obj + ".do", src) for (obj,src) in zip(profilefit_objs, profilefit_sources) ]
env_1d.Program("profilefit", profilefit_dbg_objlist)

psfconvolve_dbg_objlist = [ env_debug.Object(obj + ".do", src) for (obj,src) in zip(psfconvolve_objs, psfconvolve_sources) ]
env_debug.Program("psfconvolve", psfconvolve_dbg_objlist)

psfconvolve1d_dbg_objlist = [ env_debug.Object(obj + ".do", src) for (obj,src) in zip(psfconvolve1d_objs, psfconvolve1d_sources) ]
env_1d.Program("psfconvolve1d", psfconvolve1d_dbg_objlist)


spline2dtest_objlist = [ env_debug.Object(obj + ".do", src) for (obj,src) in zip(spline2dtest_objs, spline2dtest_sources) ]
env_debug.Program("spline2dtest", spline2dtest_objlist)


# timing: variation on makeimage designed to time image-generation and convolution
# Base files for timing:
timing_base_obj_string = """core/commandline_parser core/utilities core/image_io 
            core/config_file_parser core/add_functions core/mp_enorm core/mersenne_twister
            extra/timing_main"""
timing_base_objs = timing_base_obj_string.split()
timing_base_sources = [name + ".cpp" for name in timing_base_objs]

timing_sources = timing_base_sources + modelobject_sources + functionobject_sources

env.Program("timing", timing_sources)


# test harnesses, etc.:
# test_commandline_objlist = [ env_debug.Object(obj + ".do", src) for (obj,src) in zip(test_commandline_objs, test_commandline_sources) ]
# env_debug.Program("test_commandline", test_commandline_objlist)

# older programs
testparser_objlist = [ env_debug.Object(obj + ".do", src) for (obj,src) in zip(testparser_objs, testparser_sources) ]
env_debug.Program("testparser", testparser_objlist)
readimage_test_objs = ["readimage_test", "core/image_io"]
readimage_test_sources = [name + ".cpp" for name in readimage_test_objs]
readimage_test_dbg_objlist = [ env_debug.Object(obj + ".do", src) for (obj,src) in zip(readimage_test_objs, readimage_test_sources) ]
env.Program("readimage_test", readimage_test_dbg_objlist)
