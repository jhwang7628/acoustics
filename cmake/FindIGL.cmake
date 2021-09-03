find_path(IGL_INCLUDE_ROOT libigl/include HINTS
    $ENV{LIBIGL_ROOT}
    ${CMAKE_SOURCE_DIR}/external/libigl-src
  DOC "The directory where libigl resides"
)
set(IGL_INCLUDE_DIR "${IGL_INCLUDE_ROOT}/libigl/include")
mark_as_advanced(IGL_INCLUDE_DIR)
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(IGL DEFAULT_MSG IGL_INCLUDE_DIR) 
