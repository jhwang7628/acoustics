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

    if (argc != 4)
    {
        std::cerr << "**Usage: " << argv[0] << " <config_file> <data_dir> <max_number_files (-1 will disable)>\n"; 
        exit(1);
    }

    const char *configFile = argv[1]; 
    const char *datadir    = argv[2]; 
    const int N            = atoi(argv[3]); 

    UniformGridWithObject grid; 

    grid.Reinitialize(std::string(configFile)); 

    std::cout << grid << std::endl;

    Eigen::MatrixXd centroids; 

    grid.GetAllCellCenterPosition(centroids); 

    std::vector<std::string> filenames, filenamesAll; 
    IO::listDirectoryMatch(datadir, ".*pressure_[[:digit:]]*\\.dat", filenamesAll); 

    for (size_t ii=0; ii<filenamesAll.size(); ii++) 
    {
        if (!IO::ExistFile(filenamesAll[ii]+".vtk"))
            filenames.push_back(filenamesAll[ii]); 

        if (filenames.size()>N && N!=-1)  // early termination
            break; 
    }

    const int N_files = filenames.size(); 
    std::shared_ptr<Eigen::MatrixXd> dataBuffer(new Eigen::MatrixXd(grid.N_cells(), 3)); 
    grid.InsertCellCenteredData("dataBuffer", dataBuffer); 

    STL_Wrapper::PrintVectorContent(std::cout, filenames, 6);

    std::string absFilePath; 
    for (int ii=0; ii<N_files; ii++) 
    {
        absFilePath = IO::AssembleFilePath(std::string(datadir), filenames[ii]); 
        IO::readMatrixX<double>(*dataBuffer, absFilePath.c_str(), IO::BINARY); 
        grid.WriteVTKCellCentered(absFilePath, "dataBuffer", "pressure");
    }

    return 0;
}
