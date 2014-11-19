# Install script for directory: /home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/src/tools

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
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/tools/isostuffer/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/tools/tetviewer/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/tools/elasticity/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/tools/eigensolver/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/tools/single_tet/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/tools/fbem_input_gen/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/tools/rigidsim/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/tools/pulse/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/tools/init-rigid-tools/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/tools/soundgen/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/tools/acceleration_noise/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/tools/acoustic_transfer/cmake_install.cmake")
  include("/home/jui-hsien/Research/gpu-wavesolver-repo/acoustics/gcc-build/src/tools/impulse_response/cmake_install.cmake")

endif()

