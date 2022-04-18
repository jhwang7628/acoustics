# CMake Project file for isosurface stuffing
# Changxi Zheng (cxzheng@cs.cornell.edu)
# Copyright @ Jan 2009
project(isostuffer)

cmake_minimum_required(VERSION 3.16.0)

set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY true)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# General Configuration: find all the required libraries.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" "${CMAKE_CURRENT_SOURCE_DIR}/cmake/bk")
list(PREPEND CMAKE_MODULE_PATH
    "${PROJECT_SOURCE_DIR}/cmake/recipes/external/"
    "${PROJECT_SOURCE_DIR}/cmake/recipes/internal/"
)
list(PREPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/acoustics/")

# compile option
option(USE_QT5                                  "Turn on for QT5. Turn off for QT4"                             ON)
option(USE_DEBUG                                "Turn on the debug mode"                                        OFF)
option(USE_DEBUG_PRINT                          "Turn on debug print (can still be compiled with release mode"  OFF)
option(USE_OPENMP                               "Turn on the OpenMP feature when compiling"                     OFF)
option(USE_64BIT_ARCH                           "Compile the 64bit executable"                                  OFF)
option(USE_CUDA                                 "Use CUDA"                                                      OFF)
option(USE_VLFEAT                               "Use VLFeat for KNN"                                            OFF)
option(USE_BOOST                                "Use Boost"                                                     OFF)
option(USE_IGL                                  "Use IGL"                                                       OFF)
option(USE_MESHFIX                              "Use meshfix"                                                   OFF)

# Options for building individual tools
option(BUILD_ISOSTUFFER                         "Compile the command-line isostuffer"                           OFF)
option(BUILD_ISOSTUFFER_GUI                     "Compile the isostuffer GUI"                                    OFF)
option(BUILD_TETVIEWER                          "Compile tetmesh viewer"                                        OFF)
option(BUILD_ELASTICITY_SOLVER                  "Compiles elasticity solver"                                    OFF)
option(BUILD_EIGENSOLVER                        "Compile modal eigensolver"                                     OFF)
option(BUILD_SINGLE_TET                         "Compile single tet builder"                                    OFF)
option(BUILD_FBEM_INPUT_GEN                     "Compile FastBEM input file builder"                            OFF)
option(BUILD_RIGIDSIM                           "Compile the command-line RBD solver"                           OFF)
option(BUILD_RIGIDSIM_GUI                       "Compile rigid-body solver GUI"                                 OFF)
option(BUILD_RIGIDSIM_PLAYER                    "Compile the RBD solver replayer"                               OFF)
option(BUILD_PAT_PRECOMPUTE_GUI                 "Compile GUI for PAT precomputation"                            OFF)
option(BUILD_PAN_PRECOMPUTE_GUI                 "Compile GUI for PAN precomputation"                            OFF)
option(BUILD_RUN_PULSE_SOLVE                    "Compile command-line PAN precomputation"                       OFF)
option(BUILD_MULTITERM_PAN                      "Compile multi-term PAN system builder"                         OFF)
option(BUILD_COMPRESSED_PAN                     "Compile compressed PAN system builder"                         OFF)
option(BUILD_INIT_RIGID_TOOLS                   "Compiles tools for writing RBD data"                           OFF)
option(BUILD_WAVESOLVER_TOOLS                   "Compiles tools for wave solver data"                           OFF)
option(BUILD_RIGID_SNDGEN                       "Compiles the rigid-body sound generator"                       OFF)
option(BUILD_ESTIMATE_TIMESCALE                 "Compiles the Hertz timescale estimator"                        OFF)
option(BUILD_PYTHOFF_EXTENSIOFFS                "Compiles the Python extensions"                                OFF)
option(USE_MKL_BLAS                             "Get BLAS/LAPACK from MKL"                                      OFF)

option(BUILD_DOWNSAMPLE_IMPULSE_RESPOFFSE       "Compiles the tool for downsampling the impulse response"       OFF)
option(BUILD_PRECOMPUTE_IMPULSE_RESPOFFSE       "Compiles the impulse response wavesolver"                      OFF)
option(BUILD_PRECOMPUTE_IMPULSE_RESPOFFSE_GUI   "Compiles the gui of the impulse response wave solver"          OFF)
option(BUILD_UNIT_TESTING                       "Compiles unit testing for some modules in the repository"      OFF)
option(BUILD_DIFFERENTIATE_INTERPOLATE_HRIR     "Compiles HRIR differentiation tools"                           OFF)
option(BUILD_FDTD_ACOUSTIC_SIMULATOR_VIEWER     "Compiles FDTD Acoustic Simulator viewer"                       OFF)
option(BUILD_KIRCHHOFF_INTEGRAL_SOLVER          "Compiles Kirchhoff integral solver for acoustic transfer"      OFF)

# Options for building Kleinpat support tools
option(BUILD_KLEINPAT_TOOLS                     "Compile the KleinPAT tools" OFF)
option(BUILD_WAVESOLVER                         "Compile the wavesolver"                                        ON)

if (BUILD_KLEINPAT_TOOLS)
    message(STATUS "Building all supporting tools for KleinPAT")
    set(USE_MKL_BLAS ON)
    set(BUILD_ELASTICITY_SOLVER ON)
    set(BUILD_EIGENSOLVER ON)
    # hack before we completely refactor find_package: use FetchContent to download and use find_package
    # to include/link
    set(FETCHCONTENT_BASE_DIR "${CMAKE_SOURCE_DIR}/external" CACHE PATH "Directory under which to collect all populated content")
    include(eigen)
    include(libigl)
endif ()

if (BUILD_WAVESOLVER)
    message(STATUS "Building all tools for the wavesolver")
    set(BUILD_FDTD_ACOUSTIC_SIMULATOR_VIEWER ON)
    set(USE_GTS ON)
    set(USE_VLFEAT ON)
    set(USE_IGL ON)
    set(USE_BOOST ON)
    set(USE_GSL ON)
    set(USE_MKL_BLAS ON)
endif ()




# Try enabling diffs here...
#   This is hacked in as a way to consolidate different codebases
#   We should eventually go through and resolve all of these diffs
add_definitions(-DDIFF_DEFINE)
add_definitions(-DCGAL_INTERSECTION_VERSION=1)
## FIXME hack
# add_definitions(-DNO_DRAW)

#===================================================================
## Compiler
# set compiler flags for debug/release
if ( USE_DEBUG )
    MESSAGE(STATUS "Using debug mode")
    add_definitions(-DDEBUG)
    add_definitions(-Wall)
    add_definitions(-O0 -g -Wno-sign-compare)
    #add_definitions(-Og -g -DCGAL_DISABLE_ROUNDING_MATH_CHECK=ON -Wno-sign-compare)
    set(CMAKE_BUILD_TYPE Debug)
else ( USE_DEBUG )
    MESSAGE(STATUS "Using release mode")
    set(CMAKE_BUILD_TYPE Release)
    #add_definitions(-Wall)
    add_definitions(-g -Wno-unused-result -Wno-sign-compare -Wno-unused-variable -Wno-unused-local-typedefs)
endif ( USE_DEBUG )

if (USE_DEBUG_PRINT)
    add_definitions(-DDEBUG_PRINT)
endif (USE_DEBUG_PRINT)

if ( USE_64BIT_ARCH )
    add_definitions(-m64)
    set(LINK_FLAGS -m64)
endif ( USE_64BIT_ARCH )

#===================================================================
## Libraries
# check boost libraries
if (USE_BOOST)
    find_package(Boost 1.33 REQUIRED)
    add_definitions(-DUSE_BOOST)
endif (USE_BOOST)

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

    if (USE_QT5)
        find_package(Qt5Widgets REQUIRED)
    else (USE_QT5)
        find_package(Qt4 REQUIRED)
    endif (USE_QT5)

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
    add_definitions(-DUSE_GTS)
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

if (USE_IGL)
    add_definitions(-DUSE_IGL)
    include_directories(external/libigl/include)
endif (USE_IGL)
if (USE_MESHFIX)
    include_directories(external/meshfix)
    include_directories(external/meshfix/JMeshExt-1.0alpha_src/include)
    include_directories(external/meshfix/JMeshExt-1.0alpha_src/JMeshLib-1.2/include)
    add_definitions(-DMESHFIX_WITH_EIGEN)
    link_directories(external/meshfix/build)
    link_directories(external/meshfix/build/JMeshExt-1.0alpha_src)
    link_directories(external/meshfix/build/JMeshExt-1.0alpha_src/OpenNL3.2.1/binaries/lib)
    link_directories(external/meshfix/build/JMeshExt-1.0alpha_src/JMeshLib-1.2)
endif (USE_MESHFIX)

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

