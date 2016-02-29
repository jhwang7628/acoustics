#include <sstream> 
#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense> 

#include <unistd.h>

#include <omp.h>

#include "utils/IO/IO.h"
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

    if (argc != 3)
    {
        std::cerr << "**Usage: " << argv[0] << " <config_file> <data_dir> \n"; 
        exit(1);
    }

    const char *configFile = argv[1]; 
    const char *datadir    = argv[2]; 

    UniformGridWithObject grid; 
    grid.Reinitialize(std::string(configFile)); 

    Eigen::MatrixXd centroids; 

    grid.GetAllCellCenterPosition(centroids); 

    std::vector<std::string> filenames; 
    IO::listDirectoryMatch(datadir, ".*pressure_[[:digit:]]*\\.dat", filenames); 

    const int N_files = filenames.size(); 
    std::shared_ptr<Eigen::MatrixXd> dataBuffer(new Eigen::MatrixXd(grid.N_cells(), 3)); 
    grid.InsertCellCenteredData("dataBuffer", dataBuffer); 

    STL_Wrapper::PrintVectorContent(std::cout, filenames);

    std::string absFilePath; 
    for (int ii=0; ii<N_files; ii++) 
    {
        absFilePath = IO::AssembleFilePath(std::string(datadir), filenames[ii]); 
        IO::readMatrixX<double>(*dataBuffer, absFilePath.c_str(), IO::BINARY); 
        grid.WriteVTKCellCentered(absFilePath, "dataBuffer", "pressure");
    }

    return 0;
}
