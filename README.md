# Finite-Difference Time-Domain Acoustic Simulator

This tool is built using CMake. Building has been tested on Linux with Gnu gcc compiler only. 

## Requirements

This project uses a number of open source projects to work properly. The following list might not be comprehensive as some older part of the code has complex dependencies. Please add the libraries you have problem with to the list!

* [CMake](https://cmake.org/) (v3.0.2+)
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) (v3.2.8+) - For some linear algebra operations such as svd. It is header only library, therefore it only needs to be downloaded, and placed in an appropriate location. Specify the location of `Eigen/Core` on your machine by setting the environmental variable `EIGEN_ROOT`. For example, on Linux, do `export EIGEN_ROOT=/home/<name>/opt/eigen`.
* [Paraview](http://www.paraview.org/) - For postprocessing purpose. 

## Build Directions:
Create build directory within source directory and configure cmake
```
mkdir build
cd build
cmake ..
```
At the last step, you can specify build option such as `cmake -DUSE_OPENMP=ON`. To view all the option, its easiest to use package such as ccmake: `ccmake ..`.

### Simulator module
To build the simulator test module in parallel, do
```
make -j8 unit_testing_FDTD_AcousticSimulator
```

### Postprocessing module
The result of the run can be visualized by paraview. The simulator can be configured to dump raw data in `*.dat` files, and we have a converter to turn these files into `*.vtk` files, which can be recognized by paraview. To build the converter, do 
```
make -j8 convert-dat-vtk
```

## Run the code
The simulator is mainly configured by a XML file. A test template is given in `src/tools/unit_testing/test_FDTD_RigidObject.xml`. We can configure solve domain, object initial position, pressure sources etc. To run the simulator with the configuration, do 
```
./bin/unit_testing_FDTD_AcousticSimulator <config_file>
```

## TODOS
1. Improve ghost-cell boundary motion
