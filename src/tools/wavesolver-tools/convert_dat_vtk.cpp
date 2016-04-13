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
#include "grid/GridWithObject.h" 
#include "vtkConverter/vtkConverter.h"
#include "signalprocessing/resample.h" 

#ifndef SDUMP
#define SDUMP(x) " [" << #x << "] = " << x
#endif
#ifndef COUT_SDUMP
#define COUT_SDUMP(x) std::cout << SDUMP(x) << std::endl;
#endif

int main(int argc, char ** argv) {

    if (argc < 4)
    {
        std::cerr << "**Usage: " << argv[0] << " <config_file> <data_dir> <stop_file_index (-1 will disable)> [start_file_index]\n"; 
        exit(1);
    }

    const char *configFile  = argv[1]; 
    const char *datadir     = argv[2]; 
    const int indexStop     = atoi(argv[3]); 
    const int indexStart    = (argc==5) ? atoi(argv[4]) : 0; 

    UniformGridWithObject grid; 

    grid.Reinitialize(std::string(configFile)); 

    Eigen::MatrixXd centroids; 

    grid.GetAllCellCenterPosition(centroids); 

    std::vector<std::string> filenames, filenamesAll; 
    IO::listDirectoryMatch(datadir, ".*pressure_[[:digit:]]*\\.dat", filenamesAll); 

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
        UniformGridWithObject gridThread(grid);
        std::shared_ptr<Eigen::MatrixXd> dataBuffer(new Eigen::MatrixXd(grid.N_cells(), 3)); 
        gridThread.InsertCellCenteredData("dataBuffer", dataBuffer); 
        const std::string absFilePath = IO::AssembleFilePath(std::string(datadir), filenames[ii]); 
        IO::readMatrixX<double>(*dataBuffer, absFilePath.c_str(), IO::BINARY, 0); 
        //std::cout << " max: " << dataBuffer->maxCoeff() << std::endl;
        gridThread.WriteVTKCellCentered(absFilePath, "dataBuffer", "pressure");
        #pragma omp critical
        {
            std::cout << count << " " << std::flush;
            ++count;
        }
    }
    std::cout << std::endl;

    return 0;
}
