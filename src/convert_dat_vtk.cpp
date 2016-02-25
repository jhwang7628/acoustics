#include <sstream> 
#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense> 

#include <unistd.h>

#include <omp.h>

#include "IO/IO.h"
#include "grid/GridWithObject.h" 
#include "vtkConverter/vtkConverter.h"
#include "signalprocessing/resample.h" 

/* 
 * simple tool interpolate the hrir computed with the wavesolver given
 * interpolating positions using the grid structure 
 */
int main(int argc, char ** argv) {

    //if (argc != 14) 
    //{
    //    std::cerr << "**Usage: " << argv[0] << " <data_dir> <data_name> <Nx> <Ny> <Nz> <minBoundx> <minBoundy> <minBoundz> <maxBoundx> <maxBoundy> <maxBoundz> <listening_points> <out_file>" << std::endl;
    //    exit(1);
    //}
    //

    ///// read parameters 
    //const char * datadir  = argv[1]; 
    //const char * dataname = argv[2]; 
    //const int Nx = atoi( argv[3] ); 
    //const int Ny = atoi( argv[4] ); 
    //const int Nz = atoi( argv[5] );

    //const double minBoundx = atof( argv[6] ); 
    //const double minBoundy = atof( argv[7] ); 
    //const double minBoundz = atof( argv[8] ); 

    //const double maxBoundx = atof( argv[9] ); 
    //const double maxBoundy = atof( argv[10] );
    //const double maxBoundz = atof( argv[11]);

    //const char * listenfile = argv[12];
    //const char * outfile = argv[13];


    //Eigen::Vector3i cellCount; 
    //cellCount << Nx, Ny, Nz; 

    //Eigen::Vector3d minBound; 
    //minBound << minBoundx, minBoundy, minBoundz; 

    //Eigen::Vector3d maxBound; 
    //maxBound << maxBoundx, maxBoundy, maxBoundz; 

    ///// create grid for interpolation
    //UniformGrid<double> grid( minBound, maxBound, cellCount ); 
    //std::cout << grid << std::endl;

    ///// read IR data
    //std::vector<std::string> filenames; 
    //IO::listDirectory(datadir, dataname, filenames); 
    //std::sort(filenames.begin(), filenames.end());

    ///// read listening position 
    //Eigen::MatrixXd listeningPositions; 
    //IO::readMatrixXd( listeningPositions, listenfile, IO::BINARY ); 


    ////const int NCell=Nx*Ny*Nz; 
    //const int NCellListened=listeningPositions.rows(); 
    //const int Nts=filenames.size(); 

    ///// store the listened values 
    //Eigen::MatrixXd listenedValues = Eigen::MatrixXd::Zero( NCellListened, Nts ); 
    //int globalCount = 0;

    //#pragma omp parallel for num_threads(46) 
    //for ( int ii=0; ii<Nts; ii++ ) 
    //{

    //    UniformGrid<double> grid_thread(grid); // cannot put it as openmp clause, don't know why

    //    /// for progress tracking
    //    #pragma omp critical
    //    {
    //        globalCount += 1;
    //        std::stringstream buf; 
    //        buf << globalCount << " ";
    //        std::cout << buf.rdbuf() << std::flush;
    //    }
    //    std::string fileAndPath = IO::AssembleFilePath( datadir, filenames[ii] );

    //    std::shared_ptr<Eigen::MatrixXd> IR_step(new Eigen::MatrixXd()); // its pressure defined cell centered
    //    IO::readMatrixXd( *IR_step, fileAndPath.c_str(), IO::BINARY, 0 );

    //    // put data in grid, interpolate, and remove it right away
    //    const std::string& key = filenames[ii]; 
    //    grid_thread.InsertCellCenteredData( key, IR_step );

    //    Eigen::VectorXd valueTimestep = Eigen::VectorXd::Zero(NCellListened); 

    //    for ( int jj=0; jj<NCellListened; jj++ ) 
    //    {
    //        Eigen::VectorXd listeningRows = listeningPositions.row(jj);
    //        Eigen::VectorXd val = grid_thread.InterpolateCellCenteredData( key, listeningRows ); 
    //        if ( val.size() != 1 ) throw runtime_error("**ERROR** wrong size of interpolated data for IR" );
    //        valueTimestep(jj) = val(0);
    //    }
    //    grid_thread.DeleteCellCenteredData( filenames[ii] ); 

    //    /// synchronize threads to write to listenedvalues
    //    #pragma omp critical 
    //    listenedValues.col(ii) = valueTimestep;

    //}
    //std::cout << std::endl << std::flush;

    //const int skipEvery = 8; 
    //SIGNAL_PROCESSING::NaiveDownsample(listenedValues, skipEvery); 

    ////IO::writeMatrixXd( listenedValues, outfile, IO::BINARY );
    //IO::writeMatrixXd( listenedValues, outfile, IO::BINARY );


    return 0;
}
