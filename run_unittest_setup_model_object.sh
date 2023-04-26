#! /bin/bash

# load environment-dependent definitions for CXXTESTGEN, CPP, etc.
. ./define_unittest_vars.sh

# Predefine some ANSI color escape codes
RED='\033[0;31m'
GREEN='\033[0;0;32m'
NC='\033[0m' # No Color

# Unit tests for setup_model_object
echo
echo "Generating and compiling unit tests for setup_model_object..."
$CXXTESTGEN --error-printer -o test_runner_setup_modelobj.cpp unit_tests/unittest_setup_model_object.t.h
$CPP -std=c++11 -o test_runner_setup_modelobj test_runner_setup_modelobj.cpp core/model_object.cpp \
core/setup_model_object.cpp core/utilities.cpp core/convolver.cpp core/config_file_parser.cpp \
core/mersenne_twister.cpp core/mp_enorm.cpp core/oversampled_region.cpp core/downsample.cpp \
core/image_io.cpp core/psf_oversampling_info.cpp function_objects/psf_interpolators.cpp \
-I. -Icore -Isolvers -I/usr/local/include -Ifunction_objects -I$CXXTEST \
-L/usr/local/lib -lfftw3_threads -lcfitsio -lfftw3 -lgsl -lgslcblas -lm
if [ $? -eq 0 ]
then
  echo "Running unit tests for setup_model_object:"
  ./test_runner_setup_modelobj NewTestSuite
  exit
else
  echo -e "${RED}Compilation of unit tests for setup_model_object.cpp failed.${NC}"
  exit 1
fi
