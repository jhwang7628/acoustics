# - Try to find VLFeat
# Once done this will define
#  VLFEAT_FOUND - System has VLFeat
#  VLFEAT_INCLUDE_DIR - The VLFeat include directories
#  VLFEAT_LIBRARY - The libraries needed to use VLFeat
#  VLFEAT_DEFINITIONS - Compiler switches required for using VLFeat

find_path(VLFEAT_INCLUDE_DIR NAMES vl/generic.h
          PATHS $ENV{VLFEAT_ROOT}
          PATH_SUFFIXES vl )

find_library(VLFEAT_LIBRARY NAMES vl libvl
          PATHS $ENV{VLFEAT_ROOT}
          PATH_SUFFIXES bin/glnxa64 bin/maci64)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(VLFeat DEFAULT_MSG VLFEAT_LIBRARY VLFEAT_INCLUDE_DIR)
mark_as_advanced(VLFEAT_INCLUDE_DIR VLFEAT_LIBRARY )
