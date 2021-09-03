find_path(Eigen_INCLUDE_DIR Eigen/Core HINTS
    $ENV{EIGEN_ROOT}
    /usr/local/include/eigen3
    ${CMAKE_SOURCE_DIR}/external/eigen-src
  DOC "The directory where Eigen/Core resides"
)

mark_as_advanced(Eigen_INCLUDE_DIR)
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Eigen DEFAULT_MSG Eigen_INCLUDE_DIR)
