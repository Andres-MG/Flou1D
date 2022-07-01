cmake_minimum_required(VERSION 3.10.2)

project(CFDsolver)
enable_language(Fortran)

####################################### Definition of macros #######################################

# Type of build
set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build,
options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release
RelWithDebInfo MinSizeRel.")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS None Debug Release
    RelWithDebInfo MinSizeRel)

# Floats precision
set(PRECISION Double CACHE STRING "PRECISION can be single or double")
set_property(CACHE PRECISION PROPERTY STRINGS Single Double)

# OpenMP support (not usable yet)
# option(OPENMP "Enable OpenMP during compilation (not supported yet)" OFF)

################################### Do not edit below this line ####################################

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/obj")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if (CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    add_compile_options(-cpp)
    string(CONCAT CMAKE_Fortran_FLAGS_DEBUG "-g -Wall -Wextra -pedantic "
                                            "-fimplicit-none "
                                            "-ffpe-trap=invalid,zero,overflow "
                                            "-fcheck=all -fbacktrace")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    add_compile_options(-fpp -mkl)
    string(CONCAT CMAKE_Fortran_FLAGS_DEBUG "-g -warn all -traceback -fpe0 "
                                            "-check all -ftrapuv -debug all "
                                            "-gen-interfaces -warn interfaces")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES NVHPC)
    add_compile_options(-Mpreprocess)
    string(CONCAT CMAKE_Fortran_FLAGS_DEBUG "-g -traceback -Mbounds "
                                            "-Mchkstk -Mchkptr")
endif()

if (PRECISION STREQUAL "Single")
    add_definitions(-DSINGLE_PRECISION)
elseif (PRECISION STREQUAL "Double")
    add_definitions(-DDOUBLE_PRECISION)
endif()

# Debug flags
# https://www.nas.nasa.gov/hecc/support/kb/recommended-intel-compiler-debugging-options_92.html

# Not supported yet
#if (OPENMP)
    #if (CMAKE_Fortran_COMPILER_ID MATCHES GNU)
        #add_compile_options("-fopenmp")
        #set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")
    #elseif (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
        #add_compile_options("-qopenmp")
        #set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -qopenmp")
    #elseif (CMAKE_Fortran_COMPILER_ID MATCHES NVHPC)
        #add_compile_options("-mp")
        #set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -mp")
    #endif()
#endif()

# We can do better than this...
file(GLOB_RECURSE SOURCES src/*.f90)
add_executable(flou1d ${SOURCES})
