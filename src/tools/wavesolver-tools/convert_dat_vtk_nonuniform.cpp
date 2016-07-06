#include <sstream> 
#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense> 
#include <unistd.h>
#include <omp.h>
#include "utils/IO/IO.h"
#include "utils/IO/StringHelper.h" 
#include "utils/STL_Wrapper.h"
#include "vtkConverter/vtkConverter.h"
#include "signalprocessing/resample.h" 

#ifndef SDUMP
#define SDUMP(x) " [" << #x << "] = " << x
#endif
#ifndef COUT_SDUMP
#define COUT_SDUMP(x) std::cout << SDUMP(x) << std::endl;
#endif

int main(int argc, char ** argv) {

    if (argc < 8)
    {
        std::cerr << "**Usage: " << argv[0] << " <Nx> <Ny> <Nz> <vertex_position_dat_file> <data_dir> <p/vx/vy/vz (p for pressure, vx for x-component of velocity)> <stop_file_index (-1 will disable)> [start_file_index]\n"; 
        exit(1);
    }

    const int Nx            = atoi(argv[1]); 
    const int Ny            = atoi(argv[2]); 
    const int Nz            = atoi(argv[3]); 
    const Eigen::Vector3i dimension(Nx, Ny, Nz); 
    const char *positionFile= argv[4]; 
    const char *datadir     = argv[5]; 
    const int indexStop     = atoi(argv[7]); 
    const int indexStart    = (argc==9) ? atoi(argv[8]) : 0; 

    // get all filenames
    std::vector<std::string> filenames, filenamesAll; 
    const char *pvFlag = argv[6];
    std::string scalarName; 
    if (strcmp(pvFlag, "p")==0)
    {
        IO::listDirectoryMatch(datadir, ".*pressure_[[:digit:]]*\\.dat", filenamesAll); 
        scalarName = std::string("pressure"); 
    }
    else if (strcmp(pvFlag, "vx")==0) 
    {
        IO::listDirectoryMatch(datadir, ".*velocity_0_[[:digit:]]*\\.dat", filenamesAll); 
        scalarName = std::string("velocity_x"); 
    }
    else if (strcmp(pvFlag, "vy")==0) 
    {
        IO::listDirectoryMatch(datadir, ".*velocity_1_[[:digit:]]*\\.dat", filenamesAll); 
        scalarName = std::string("velocity_y"); 
    }
    else if (strcmp(pvFlag, "vz")==0) 
    {
        IO::listDirectoryMatch(datadir, ".*velocity_2_[[:digit:]]*\\.dat", filenamesAll); 
        scalarName = std::string("velocity_z"); 
    }
    else 
        throw std::runtime_error("**ERROR** wrong p/v flag for the input argument. Valid argument: p / vx / vy / vz"); 

    // read position for the data. 
    Eigen::MatrixXd position; 
    IO::readMatrixX<double>(position, positionFile, IO::BINARY, 0); 

    for (size_t ii=0; ii<filenamesAll.size(); ii++) 
    {
        if (((int)ii>indexStop || (int)ii<indexStart) && indexStop!=-1)  // early termination
            continue; 
        if (!IO::ExistFile(filenamesAll[ii]+".vtk"))
            filenames.push_back(filenamesAll[ii]); 
    }

    const int N_files = filenames.size(); 
    STL_Wrapper::PrintVectorContent(std::cout, filenames, 6);
    int count=0;
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(static) default(shared)
    #endif
    for (int ii=0; ii<N_files; ii++) 
    {
        // read the data 
        Eigen::MatrixXd data; 
        const std::string absFilePath = IO::AssembleFilePath(std::string(datadir), filenames[ii]); 
        IO::readMatrixX<double>(data, absFilePath.c_str(), IO::BINARY, 0); 

        // convert the data into vtk and write it
        const std::string filename = filenames[ii] + std::string(".vtk"); 
        VTKConverter::VTKStructureGridWithScalarFromEigen(position, data, filename, scalarName, VTKConverter::BINARY, dimension); 
        #pragma omp critical
        {
            std::cout << count << " " << std::flush;
            ++count;
        }
    }
    std::cout << std::endl;

    return 0;
}
