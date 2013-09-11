# CMake Project file for isosurface stuffing  
# Changxi Zheng (cxzheng@cs.cornell.edu)
# Copyright @ Jan 2009
project(isostuffer)

cmake_minimum_required(VERSION 2.6)

# General Configuration: find all the required libraries.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# compile option
option(USE_DEBUG "Turn on the debug mode" OFF)
option(USE_OPENMP "Turn on the OpenMP feature when compiling" OFF)
option(USE_64BIT_ARCH "Compile the 64bit executable" OFF)

# Options for building individual tools
option(BUILD_ISOSTUFFER         "Compile the command-line isostuffer" ON)
option(BUILD_ISOSTUFFER_GUI     "Compile the isostuffer GUI" ON)
option(BUILD_TETVIEWER          "Compile tetmesh viewer" ON)
option(BUILD_ELASTICITY_SOLVER  "Compiles elasticity solver" ON)
option(BUILD_EIGENSOLVER        "Compile modal eigensolver" ON)
option(BUILD_SINGLE_TET         "Compile single tet builder" ON)
option(BUILD_FBEM_INPUT_GEN     "Compile FastBEM input file builder" ON)
option(BUILD_RIGIDSIM_GUI       "Compile rigid-body solver GUI" ON)
option(BUILD_PAN_PRECOMPUTE_GUI "Compile GUI for PAN precomputation" ON)
option(BUILD_RUN_PULSE_SOLVE    "Compile command-line PAN precomputation" ON)
option(BUILD_MULTITERM_PAN      "Compile multi-term PAN system builder" ON)
option(BUILD_COMPRESSED_PAN     "Compile compressed PAN system builder" ON)
option(BUILD_INIT_RIGID_TOOLS   "Compiles tools for writing RBD data" ON)
option(BUILD_RIGID_SNDGEN       "Compiles the rigid-body sound generator" ON)
option(BUILD_ESTIMATE_TIMESCALE "Compiles the Hertz timescale estimator" ON)

# Try enabling diffs here...
#   This is hacked in as a way to consolidate different codebases
#   We should eventually go through and resolve all of these diffs
add_definitions(-DDIFF_DEFINE)

#===================================================================
## Compiler
# set compiler flags for debug/release
if ( USE_DEBUG )
    MESSAGE(STATUS "Using debug mode")
    add_definitions(-DDEBUG)
    set(CMAKE_BUILD_TYPE Debug)
else ( USE_DEBUG )
    MESSAGE(STATUS "Using release mode")
    set(CMAKE_BUILD_TYPE Release)
endif ( USE_DEBUG )

if ( USE_64BIT_ARCH )
    add_definitions(-m64)
    set(LINK_FLAGS -m64)
endif ( USE_64BIT_ARCH )

#===================================================================
## Libraries
# check boost libraries
find_package(Boost 1.33 REQUIRED)

option(USE_MKL_BLAS             "Get BLAS/LAPACK from MKL" ON)

if ( USE_MKL_BLAS )
  # check for Intel MKL
  find_package(MKL REQUIRED)
  set(BLAS_LAPACK_LIBS mkl_intel_lp64 mkl_intel_thread mkl_core iomp5)
  add_definitions(-DUSE_MKL)
  add_definitions(-DUSE_LAPACKE)
else ()
  # Try using default libraries
  set(BLAS_LAPACK_LIBS blas lapack)
endif ()

# ARPACK for eigensolvers
find_package(ARPACK REQUIRED)

# SuperLU is also required for the eigensolver
#
# The prepackaged version that can be installed with Ubuntu doesn't
# seem to work.  Try downloading and compiling the new version (v4.3).
# 
# When setting up the cmake project, use the flag
#   -DSUPERLU_LIBRARIES=<path to version 4.3 SuperLU library>
#     NOTE: for now, includes aren't needed
#   -DSUPERLU_INCLUDES=<SuperLU include directory>
find_package(SuperLU REQUIRED)

if ( BUILD_ISOSTUFFER_GUI OR BUILD_TETVIEWER )
    # package for opengl and glut
    find_package(OpenGL REQUIRED)
    find_package(GLUT REQUIRED)

    find_package(Qt4 REQUIRED)

    # check the existence of libQGLViewer
    find_package(QGLViewer 1.0 REQUIRED)
endif ( BUILD_ISOSTUFFER_GUI OR BUILD_TETVIEWER )

# LibConfig needed for configuration files in certain projects
find_package(LibConfig 1.0 REQUIRED)

# Whether or not we need the GNU triangulated surface (GTS) package
if ( BUILD_PAN_PRECOMPUTE_GUI
    OR BUILD_RUN_PULSE_SOLVE
    OR BUILD_MULTITERM_PAN 
    OR BUILD_ESTIMATE_TIMESCALE )
    option(USE_GTS    "Compile with GTS"  ON)
    option(USE_CGAL   "Compile with CGAL" ON) 
    option(USE_GSL    "Compile with GSL"  ON)
else ( BUILD_PAN_PRECOMPUTE_GUI
       OR BUILD_RUN_PULSE_SOLVE
       OR BUILD_MULTITERM_PAN 
       OR BUILD_ESTIMATE_TIMESCALE )
    option(USE_GTS    "Compile with GTS"  OFF)
    option(USE_CGAL   "Compile with CGAL" OFF) 
    option(USE_GSL    "Compile with GSL"  OFF)
endif ( BUILD_PAN_PRECOMPUTE_GUI
        OR BUILD_RUN_PULSE_SOLVE
        OR BUILD_MULTITERM_PAN 
        OR BUILD_ESTIMATE_TIMESCALE )

if ( USE_GTS )
    find_package(GTS REQUIRED)

    # GTS requires GLIB
    find_package(GLIB 2.4 REQUIRED)
endif ( USE_GTS )

if ( USE_CGAL )
    find_package(CGAL REQUIRED)

    add_definitions(-frounding-math)
    add_definitions(-DUSE_CGAL)
endif ( USE_CGAL )

if ( USE_GSL )
    find_package(GSL REQUIRED)
endif ( USE_GSL )


if ( USE_DEBUG )
  add_definitions(-Wall)
endif ( USE_DEBUG )

if ( USE_OPENMP )
    add_definitions(-DUSE_OPENMP)
endif ( USE_OPENMP )

set(LINK_FLAGS " ")

# check if we're using Intel's compiler
if ( CMAKE_CXX_COMPILER MATCHES ".*icpc$" )
    add_definitions(-wd981 -wd383 -wd444 -wd1224 -wd1572)
    if ( NOT USE_DEBUG )
        set(CMAKE_CXX_FLAGS_RELEASE "-ipo -O3 -no-prec-div -xHost -DNDEBUG")
    endif ( NOT USE_DEBUG )

    if ( USE_OPENMP )
        find_package(Threads 1.0 REQUIRED)
        set(ICC_LIBS iomp5 ${CMAKE_THREAD_LIBS_INIT})
        add_definitions(-openmp)
    endif ( USE_OPENMP )
else ( CMAKE_CXX_COMPILER MATCHES ".*icpc$" )
    if ( USE_OPENMP )
        add_definitions(-fopenmp)
    endif ( USE_OPENMP )

    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif ( CMAKE_CXX_COMPILER MATCHES ".*icpc$" )

add_subdirectory(src)

