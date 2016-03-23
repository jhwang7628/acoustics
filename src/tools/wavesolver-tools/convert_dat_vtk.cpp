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

    //COUT_SDUMP(grid.GetCellCenterPosition(50,21,50)); 
    //COUT_SDUMP(grid.FlattenIndicies(50,21,50)); // 502150
    //COUT_SDUMP(grid.GetCellCenterPosition(100,42,100)); 
    //COUT_SDUMP(grid.FlattenIndicies(100,42,100)); // 502150

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
    std::shared_ptr<Eigen::MatrixXd> dataBuffer(new Eigen::MatrixXd(grid.N_cells(), 3)); 
    grid.InsertCellCenteredData("dataBuffer", dataBuffer); 

    STL_Wrapper::PrintVectorContent(std::cout, filenames, 6);

    std::string absFilePath; 
    for (int ii=0; ii<N_files; ii++) 
    {
        absFilePath = IO::AssembleFilePath(std::string(datadir), filenames[ii]); 
        IO::readMatrixX<double>(*dataBuffer, absFilePath.c_str(), IO::BINARY, 0); 
        std::cout << " max: " << dataBuffer->maxCoeff() << std::endl;
        grid.WriteVTKCellCentered(absFilePath, "dataBuffer", "pressure");
    }

    return 0;
}
