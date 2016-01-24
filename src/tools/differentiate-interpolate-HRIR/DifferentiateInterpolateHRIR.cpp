#include <sstream> 
#include <iostream>
#include <fstream>
#include <string>
#define EIGEN_RUNTIME_NO_MALLOC
#include <Eigen/Dense> 

#include <unistd.h>

#include <omp.h>

#include <boost/timer/timer.hpp>

#include "IO/IO.h"
#include "Grid.h" 

#define OMP_THREADS 36

inline double nanosecond_to_double( const boost::timer::nanosecond_type & t ) 
{
    return static_cast<double>(t)/1.0E9; 
}
    

/* 
 * compute d^2 G / dx dx on uniform grid and interpolation
 */
int main(int argc, char ** argv) {

    if (argc != 14) 
    {
        std::cerr << "**Usage: " << argv[0] << " <data_dir> <data_name> <Nx> <Ny> <Nz> <minBoundx> <minBoundy> <minBoundz> <maxBoundx> <maxBoundy> <maxBoundz> <listening_points> <out_file_prefix> (need 14 you have " << std::to_string(argc) << ")" << std::endl;
        exit(1);
    }


    /// read parameters 
    const char * datadir  = argv[1]; 
    const char * dataname = argv[2]; 
    const int Nx = atoi( argv[3] ); 
    const int Ny = atoi( argv[4] ); 
    const int Nz = atoi( argv[5] );

    const double minBoundx = atof( argv[6] ); 
    const double minBoundy = atof( argv[7] ); 
    const double minBoundz = atof( argv[8] ); 

    const double maxBoundx = atof( argv[9] ); 
    const double maxBoundy = atof( argv[10] );
    const double maxBoundz = atof( argv[11]);

    const char * listenfile = argv[12];
    std::string outPrefix = std::string(argv[13]);

    const double dx = ( maxBoundx - minBoundx ) / Nx; 
    const double dy = ( maxBoundy - minBoundy ) / Ny; 
    const double dz = ( maxBoundz - minBoundz ) / Nz; 


    Eigen::Vector3i cellCount; 
    cellCount << Nx, Ny, Nz; 

    Eigen::Vector3d minBound; 
    minBound << minBoundx, minBoundy, minBoundz; 

    Eigen::Vector3d maxBound; 
    maxBound << maxBoundx, maxBoundy, maxBoundz; 

    /// create grid for interpolation
    UniformGrid grid( minBound, maxBound, cellCount ); 
    std::cout << grid << std::endl;

    /// read IR data
    std::vector<std::string> filenames; 
    IO::listDirectory(datadir, dataname, filenames); 
    std::sort(filenames.begin(), filenames.end());

    /// read listening position 
    std::cout << "reading listening matrix : " << listenfile << std::endl;
    Eigen::MatrixXd listeningPositions; 
    IO::readMatrixXd( listeningPositions, listenfile, IO::BINARY, 0 ); 

    const int NCell=Nx*Ny*Nz; 
    const int NCellListened=listeningPositions.rows(); 
    const int Nts=filenames.size(); 

    omp_set_num_threads(OMP_THREADS); 
    std::cout << "start application with " << OMP_THREADS << " threads " << std::endl;


    std::cout << "initializing buffer matrices. might take a while..." << std::endl;
    /// store the listened values 
    std::vector<std::shared_ptr<Eigen::MatrixXd>> pressureBuffer(OMP_THREADS); // pressure read from disk
    std::vector<std::vector<std::shared_ptr<Eigen::MatrixXd>>> gradientBuffer(OMP_THREADS); 
    std::vector<std::vector<std::shared_ptr<Eigen::MatrixXd>>> hessianBuffer(OMP_THREADS); 


    std::cout << " check memory now" << std::endl; 
    usleep(3E6); 
    std::cout << " allocating listen values " << std::endl;
    Eigen::MatrixXd listenedValues_xx = Eigen::MatrixXd::Zero(NCellListened, Nts); 
    Eigen::MatrixXd listenedValues_xy = Eigen::MatrixXd::Zero(NCellListened, Nts); 
    Eigen::MatrixXd listenedValues_xz = Eigen::MatrixXd::Zero(NCellListened, Nts); 
    Eigen::MatrixXd listenedValues_yy = Eigen::MatrixXd::Zero(NCellListened, Nts); 
    Eigen::MatrixXd listenedValues_yz = Eigen::MatrixXd::Zero(NCellListened, Nts); 
    Eigen::MatrixXd listenedValues_zz = Eigen::MatrixXd::Zero(NCellListened, Nts); 

    std::cout << " allocating buffers " << std::endl;

    for (int jj=0; jj<OMP_THREADS; jj++)
    {
        std::cout << "  - for thread : " << jj << std::endl;
        pressureBuffer[jj].reset(new Eigen::MatrixXd(NCell,1)); 
        pressureBuffer[jj]->setZero(); 
        gradientBuffer[jj].resize(3);
        hessianBuffer[jj].resize(6);
        for (int ii=0; ii<3; ii++) 
        {
            gradientBuffer[jj][ii].reset(new Eigen::MatrixXd(NCell,1)); 
            gradientBuffer[jj][ii]->setZero(); 
        }
        for (int ii=0; ii<6; ii++) 
        {
            hessianBuffer[jj][ii].reset(new Eigen::MatrixXd(NCell,1)); 
            hessianBuffer[jj][ii]->setZero(); 
        }
    }
    std::cout << " check memory again" << std::endl; 
    usleep(3E6); 

    std::cout << "start computing hessian" << std::endl;

    int global_count=0;

    #pragma omp parallel for
    for ( int pp=0; pp<Nts; pp++ ) 
    {

        boost::timer::cpu_timer timer; 
        boost::timer::cpu_times dataReadCPUTime; 
        boost::timer::cpu_times hessianCPUTime; 
        boost::timer::cpu_times interpolationCPUTime; 
        boost::timer::nanosecond_type lastTime(0); 
        boost::timer::nanosecond_type dataReadTime; 
        boost::timer::nanosecond_type hessianTime; 
        boost::timer::nanosecond_type interpolationTime; 

        const int thread_id = omp_get_thread_num(); 
        UniformGrid<double> grid_thread(grid); 

        if (thread_id==0)
            timer.start(); 

        #pragma omp critical
        {
            global_count ++; 
            std::cout << global_count << " " << std::flush; 
        }

        //std::cout << "Computing file number : " << pp << std::endl; 
        //boost::timer::auto_cpu_timer timerHessian(" time elapsed : %w seconds\n"); 

        /// for progress tracking
        std::string fileAndPath = IO::AssembleFilePath( datadir, filenames[pp] );

        IO::readMatrixX<double>( *(pressureBuffer[thread_id]), fileAndPath.c_str(), IO::BINARY, 0 );

        // put data in grid, interpolate, and remove it right away
        const std::string& key = filenames[pp]; 
        grid_thread.InsertCellCenteredData( key, pressureBuffer[thread_id] );

        ///// differentiate /////
        if (global_count%OMP_THREADS==0)
        {
            dataReadCPUTime = timer.elapsed(); 
            dataReadTime = dataReadCPUTime.wall - lastTime; 
            lastTime = dataReadCPUTime.wall; 
        }

        grid_thread.CellCenteredScalarHessian( key, gradientBuffer[thread_id], hessianBuffer[thread_id] ); 

        if (global_count%OMP_THREADS==0)
        {
            hessianCPUTime = timer.elapsed(); 
            hessianTime = hessianCPUTime.wall - lastTime; 
            lastTime = hessianCPUTime.wall; 
        }



        // debug
        //std::ofstream ofxx("data/differentiated_xx_" + std::to_string(pp) + ".txt"); 
        //std::ofstream ofxy("data/differentiated_xy_" + std::to_string(pp) + ".txt"); 
        //std::ofstream ofxz("data/differentiated_xz_" + std::to_string(pp) + ".txt"); 
        //std::ofstream ofyy("data/differentiated_yy_" + std::to_string(pp) + ".txt"); 
        //std::ofstream ofyz("data/differentiated_yz_" + std::to_string(pp) + ".txt"); 
        //std::ofstream ofzz("data/differentiated_zz_" + std::to_string(pp) + ".txt"); 
        //for (int ii=0; ii<Nx; ii++) 
        //{
        //    for (int jj=0; jj<Ny; jj++) 
        //    {
        //        for (int kk=0; kk<Nz; kk++) 
        //        {
        //            double x = minBound[0] + ii*(maxBoundx-minBoundx)/Nx + dx/2.; 
        //            double y = minBound[1] + jj*(maxBoundy-minBoundy)/Ny + dy/2.; 
        //            double z = minBound[2] + kk*(maxBoundz-minBoundz)/Nz + dz/2.; 
        //            ofxx << x << "," << y << "," << z << "," << (*hessianBuffer[0])(ii+jj*Nx+kk*Ny*Nz) << std::endl;
        //            ofxy << x << "," << y << "," << z << "," << (*hessianBuffer[1])(ii+jj*Nx+kk*Ny*Nz) << std::endl;
        //            ofxz << x << "," << y << "," << z << "," << (*hessianBuffer[2])(ii+jj*Nx+kk*Ny*Nz) << std::endl;
        //            ofyy << x << "," << y << "," << z << "," << (*hessianBuffer[3])(ii+jj*Nx+kk*Ny*Nz) << std::endl;
        //            ofyz << x << "," << y << "," << z << "," << (*hessianBuffer[4])(ii+jj*Nx+kk*Ny*Nz) << std::endl;
        //            ofzz << x << "," << y << "," << z << "," << (*hessianBuffer[5])(ii+jj*Nx+kk*Ny*Nz) << std::endl;
        //        }
        //    }
        //}
        //ofxx.close(); 
        //ofxy.close(); 
        //ofxz.close(); 
        //ofyy.close(); 
        //ofyz.close(); 
        //ofzz.close(); 

        ///// interpolate /////
        // add hessian for interpolation
        grid_thread.InsertCellCenteredData( key+"_xx", hessianBuffer[thread_id][0] ); 
        grid_thread.InsertCellCenteredData( key+"_xy", hessianBuffer[thread_id][1] ); 
        grid_thread.InsertCellCenteredData( key+"_xz", hessianBuffer[thread_id][2] ); 
        grid_thread.InsertCellCenteredData( key+"_yy", hessianBuffer[thread_id][3] ); 
        grid_thread.InsertCellCenteredData( key+"_yz", hessianBuffer[thread_id][4] ); 
        grid_thread.InsertCellCenteredData( key+"_zz", hessianBuffer[thread_id][5] ); 


        for ( int jj=0; jj<NCellListened; jj++ ) 
        {
            const double buffer_xx = grid_thread.InterpolateCellCenteredData( key+"_xx", listeningPositions.row(jj) )(0); 
            const double buffer_xy = grid_thread.InterpolateCellCenteredData( key+"_xy", listeningPositions.row(jj) )(0); 
            const double buffer_xz = grid_thread.InterpolateCellCenteredData( key+"_xz", listeningPositions.row(jj) )(0); 
            const double buffer_yy = grid_thread.InterpolateCellCenteredData( key+"_yy", listeningPositions.row(jj) )(0); 
            const double buffer_yz = grid_thread.InterpolateCellCenteredData( key+"_yz", listeningPositions.row(jj) )(0); 
            const double buffer_zz = grid_thread.InterpolateCellCenteredData( key+"_zz", listeningPositions.row(jj) )(0); 
            #pragma omp critical
            {
                listenedValues_xx(jj,pp) = buffer_xx;
                listenedValues_xy(jj,pp) = buffer_xy;
                listenedValues_xz(jj,pp) = buffer_xz;
                listenedValues_yy(jj,pp) = buffer_yy;
                listenedValues_yz(jj,pp) = buffer_yz;
                listenedValues_zz(jj,pp) = buffer_zz;
            }
        }

        if (global_count%OMP_THREADS==0)
        {
            interpolationCPUTime = timer.elapsed(); 
            interpolationTime = interpolationCPUTime.wall - lastTime; 
            lastTime = interpolationCPUTime.wall; 
        }

        if (pp%20==0) 
        {
            grid_thread.WriteVTKCellCentered( "dGdxx_"+std::to_string(pp), key+"_xx", "dGdxx");
            grid_thread.WriteVTKCellCentered( "dGdxy_"+std::to_string(pp), key+"_xy", "dGdxy");
            grid_thread.WriteVTKCellCentered( "dGdxz_"+std::to_string(pp), key+"_xz", "dGdxz");
            grid_thread.WriteVTKCellCentered( "dGdyy_"+std::to_string(pp), key+"_yy", "dGdyy");
            grid_thread.WriteVTKCellCentered( "dGdyz_"+std::to_string(pp), key+"_yz", "dGdyz");
            grid_thread.WriteVTKCellCentered( "dGdzz_"+std::to_string(pp), key+"_zz", "dGdzz");
        }

        grid_thread.DeleteCellCenteredData(key+"_xx"); 
        grid_thread.DeleteCellCenteredData(key+"_xy"); 
        grid_thread.DeleteCellCenteredData(key+"_xz"); 
        grid_thread.DeleteCellCenteredData(key+"_yy"); 
        grid_thread.DeleteCellCenteredData(key+"_yz"); 
        grid_thread.DeleteCellCenteredData(key+"_zz"); 

        #pragma omp critical
        {
            if (thread_id == 0) 
            {
                std::cout << std::endl; 
                std::cout << "Timing: " << std::endl;
                std::cout << " data read takes "           << nanosecond_to_double(dataReadTime     ) << std::endl;
                std::cout << " hessian computation takes " << nanosecond_to_double(hessianTime      ) << std::endl;
                std::cout << " interpolation takes "       << nanosecond_to_double(interpolationTime) << std::endl;
            }
        }

    }
    std::cout << std::endl;

    IO::writeMatrixX<double>(listenedValues_xx, (outPrefix+"_xx.dat").c_str(), IO::BINARY); 
    IO::writeMatrixX<double>(listenedValues_xy, (outPrefix+"_xy.dat").c_str(), IO::BINARY); 
    IO::writeMatrixX<double>(listenedValues_xz, (outPrefix+"_xz.dat").c_str(), IO::BINARY); 
    IO::writeMatrixX<double>(listenedValues_yy, (outPrefix+"_yy.dat").c_str(), IO::BINARY); 
    IO::writeMatrixX<double>(listenedValues_yz, (outPrefix+"_yz.dat").c_str(), IO::BINARY); 
    IO::writeMatrixX<double>(listenedValues_zz, (outPrefix+"_zz.dat").c_str(), IO::BINARY); 
    //IO::writeMatrixXd( listenedValues, outfile, IO::BINARY );
    //IO::writeMatrixXd( listenedValues, outfile, IO::ASCII );

    return 0;
}
