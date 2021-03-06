# Finite-Difference Time-Domain Acoustic Simulator

This tool is built using CMake. Building has been tested on Linux with Gnu gcc compiler only. 

## Requirements

This project uses a number of open source projects to work properly. The following list might not be comprehensive as some older part of the code has complex dependencies. Please add the libraries you have problem with to the list!

### Required packages: 
* [CMake](https://cmake.org/) (v3.0.2+)
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) (v3.2.8+) - For some linear algebra operations such as svd. It is header only library, therefore it only needs to be downloaded, and placed in an appropriate location. Specify the location of `Eigen/Core` on your machine by setting the environmental variable `EIGEN_ROOT`. For example, on Linux, do `export EIGEN_ROOT=/home/<name>/opt/eigen`.

### Optional packages: 
* [Paraview](http://www.paraview.org/) - For postprocessing purpose. 
* [VLFeat](http://www.vlfeat.org/index.html) - For faster NN lookup in the TriangleMeshKDTree class. This class is only available if the flag `USE_VLFEAT` is on. As a result, this is an optional package. If CMake cannot locate this library and an error is thrown, set environmental variable `VLFEAT_ROOT` to the install directory. For example, on Linux, do `export VLFEAT_ROOT=/home/<name>/opt/src/vlfeat-0.9.20`.
* [libigl](https://github.com/libigl/libigl) - Libigl is used to replace GTS on computing mesh curvature, since it was found that GTS might fail when computing mean curvatures on tet extracted mesh (used for acceleration noise). A [fork](https://github.com/jhwang7628/libigl) from this library is used with some minor changes in namespace prefix to fix name collision with our own Matrix library and to keep it as lightweight as possible. This is currently integrated into the acoustic repo as submodules. If you are cloning the repo for the first time, you will have the directory `src/libigl` but with no contents in it. To initialize, 
``` 
cd external/libigl
git submodule init
git submodule update
```
Alternatively, if recursive option was on when cloning then no need for the above step
```
git clone --recursive git@bitbucket.org:jw969/acoustics.git
```
More on git submodules: [link](https://git-scm.com/book/en/v2/Git-Tools-Submodules). 

## Build Directions:
Create build directory within source directory and configure cmake
```
mkdir build_release
cd build_release
cmake ..
```
At the last step, you can specify build option such as `cmake -DUSE_OPENMP=ON`. To view all the option, its easiest to use package such as ccmake: `ccmake ..`.

To make a debug build, do the following:
```
mkdir build_debug
cd build_debug
cmake -DUSE_DEBUG=ON ..
```

### Simulator module
To build the simulator test module in parallel, do
```
make -j40 unit_testing_FDTD_AcousticSimulator
```

### Postprocessing module
The result of the run can be visualized by paraview. The simulator can be configured to dump raw data in `*.dat` files, and we have a converter to turn these files into `*.vtk` files, which can be recognized by paraview. To build the converter, do 
```
make -j40 convert-dat-vtk
```

## Run the code
The simulator is mainly configured by a XML file. A test template is given in `src/tools/unit_testing/test_FDTD_RigidObject.xml`. We can configure solve domain, object initial position, pressure sources etc. To run the simulator with the configuration, do 
```
./bin/unit_testing_FDTD_AcousticSimulator [config_file]
```
