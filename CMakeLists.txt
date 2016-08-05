# CMake Project file for isosurface stuffing  
# Changxi Zheng (cxzheng@cs.cornell.edu)
# Copyright @ Jan 2009
project(isostuffer)

cmake_minimum_required(VERSION 3.0.2)

# set(CMAKE_C_COMPILER "/usr/bin/gcc-4.6")
# set(CMAKE_CXX_COMPILER "/usr/bin/g++-4.6")

# set(CMAKE_C_COMPILER "/usr/bin/gcc")
# set(CMAKE_CXX_COMPILER "/usr/bin/g++")

set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY true)
SET(CMAKE_CXX_FLAGS "-std=c++11")

# General Configuration: find all the required libraries.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" "${CMAKE_CURRENT_SOURCE_DIR}/cmake/bk")

# compile option
option(USE_QT5 "Turn on for QT5. Turn off for QT4" ON)
option(USE_DEBUG "Turn on the debug mode" OFF)
option(USE_OPENMP "Turn on the OpenMP feature when compiling" ON)
option(USE_64BIT_ARCH "Compile the 64bit executable" OFF)
option(USE_CUDA "Use CUDA" OFF)
option(USE_VLFEAT "Use VLFeat for KNN" ON)

# Options for building individual tools
option(BUILD_ISOSTUFFER         "Compile the command-line isostuffer" ON)
option(BUILD_ISOSTUFFER_GUI     "Compile the isostuffer GUI" ON)
option(BUILD_TETVIEWER          "Compile tetmesh viewer" ON)
option(BUILD_ELASTICITY_SOLVER  "Compiles elasticity solver" ON)
option(BUILD_EIGENSOLVER        "Compile modal eigensolver" ON)
option(BUILD_SINGLE_TET         "Compile single tet builder" ON)
option(BUILD_FBEM_INPUT_GEN     "Compile FastBEM input file builder" ON)
option(BUILD_RIGIDSIM           "Compile the command-line RBD solver" ON)
option(BUILD_RIGIDSIM_GUI       "Compile rigid-body solver GUI" ON)
option(BUILD_RIGIDSIM_PLAYER    "Compile the RBD solver replayer" ON)
option(BUILD_PAT_PRECOMPUTE_GUI "Compile GUI for PAT precomputation" ON)
option(BUILD_PAN_PRECOMPUTE_GUI "Compile GUI for PAN precomputation" ON)
option(BUILD_RUN_PULSE_SOLVE    "Compile command-line PAN precomputation" ON)
option(BUILD_MULTITERM_PAN      "Compile multi-term PAN system builder" ON)
option(BUILD_COMPRESSED_PAN     "Compile compressed PAN system builder" ON)
option(BUILD_INIT_RIGID_TOOLS   "Compiles tools for writing RBD data" ON)
option(BUILD_WAVESOLVER_TOOLS   "Compiles tools for wave solver data" ON)
option(BUILD_RIGID_SNDGEN       "Compiles the rigid-body sound generator" ON)
option(BUILD_ESTIMATE_TIMESCALE "Compiles the Hertz timescale estimator" ON)
option(BUILD_PYTHON_EXTENSIONS  "Compiles the Python extensions"    OFF)

option(BUILD_DOWNSAMPLE_IMPULSE_RESPONSE     "Compiles the tool for downsampling the impulse response" ON)
option(BUILD_PRECOMPUTE_IMPULSE_RESPONSE     "Compiles the impulse response wavesolver" ON)
option(BUILD_PRECOMPUTE_IMPULSE_RESPONSE_GUI "Compiles the gui of the impulse response wave solver" ON)
option(BUILD_UNIT_TESTING                    "Compiles unit testing for some modules in the repository" ON) 
option(BUILD_DIFFERENTIATE_INTERPOLATE_HRIR  "Compiles HRIR differentiation tools" ON) 
option(BUILD_FDTD_ACOUSTIC_SIMULATOR_VIEWER  "Compile FDTD Acoustic Simulator viewer" ON)

# Try enabling diffs here...
#   This is hacked in as a way to consolidate different codebases
#   We should eventually go through and resolve all of these diffs
add_definitions(-DDIFF_DEFINE)
add_definitions(-DCGAL_INTERSECTION_VERSION=1)

#===================================================================
## Compiler
# set compiler flags for debug/release
if ( USE_DEBUG )
    MESSAGE(STATUS "Using debug mode")
    add_definitions(-DDEBUG)
    add_definitions(-Wall)
    set(CMAKE_BUILD_TYPE Debug)
else ( USE_DEBUG )
    MESSAGE(STATUS "Using release mode")
    set(CMAKE_BUILD_TYPE Release)
    add_definitions(-Wall)
endif ( USE_DEBUG )

if ( USE_64BIT_ARCH )
    add_definitions(-m64)
    set(LINK_FLAGS -m64)
endif ( USE_64BIT_ARCH )

#===================================================================
## Libraries
# check boost libraries
find_package(Boost 1.33 REQUIRED)
add_definitions(-DUSE_BOOST)

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
find_package(Eigen REQUIRED)

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
find_package(LibConfig++ 1.0 REQUIRED)

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


if ( USE_OPENMP )
    find_package(OpenMP)
    add_definitions(-DUSE_OPENMP)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ( USE_OPENMP )

set(LINK_FLAGS " ")

if (USE_CUDA)
    find_package(CUDA)
    add_definitions(-DUSE_CUDA)
endif (USE_CUDA)

if (USE_VLFEAT)
    find_package(VLFeat REQUIRED)
    add_definitions(-DUSE_VLFeat)
    include_directories(${VLFEAT_INCLUDE_DIR})
endif (USE_VLFEAT)

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


SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

add_subdirectory(src)

