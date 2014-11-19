# Install script for directory: /home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/arpack++/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/isostuffer/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/utils/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/linearalgebra/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/generic/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/logging/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/deformable/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/io/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/rigid/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/field/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/wavesolver/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/superlu-interface/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/eigensolver/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/distancefield/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/math/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/geometry/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/transfer/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/sndgen/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/multipole/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/scripts/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/tools/cmake_install.cmake")

endif()

