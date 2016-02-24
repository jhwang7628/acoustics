#include "OpenfoamMesh.h" 
#include "Grid.h" 
#include "IO/IO.h"
#include "signalprocessing/FilterDesign.h"
#include "IFFT_Synthesis.h"
#include <boost/timer/timer.hpp>

#include <unistd.h>



/// read and analyze the frequency response of IR
void SpectralAnalysis_ImpulseResponse()
{

    const int N_FFT = 256; 
    const std::string dirHRIR("/home/jui-hsien/code/acoustics/work/impulse-response/data_Gaussian0p0025_res250_offset2stddev_simulation"); 
    const std::string searchString("head_right_pressure_"); 

    const std::string dataDirectory = IO::AssembleFilePath(dirHRIR,"out"); 
    IO::CreateDirectoryForce(dataDirectory); 

    std::cout << "analyzing spectral component of HRIR at directory: " << dirHRIR << " with start string " << searchString << std::endl;

    std::vector<std::string> filenames; 
    IO::listDirectory(dirHRIR.c_str(), searchString.c_str(), filenames); 
    for (std::string &s : filenames) s = IO::AssembleFilePath(dirHRIR,s); 

    const int N_timesteps = filenames.size(); 



    std::cout << " number of files found : " << N_timesteps << std::endl;

    // define grid for projection 
    const int gridDimension = 250;
    const int gridDimension3= gridDimension*gridDimension*gridDimension;
    Eigen::Vector3d minBound(-0.2178337591671500,-0.2017749603171500,-0.2192579190471500); 
    Eigen::Vector3d maxBound(0.2204221281671500,0.2364809270171500,0.2189979682871500); 
    Eigen::Vector3i cellCount(gridDimension,gridDimension,gridDimension); 
    UniformGrid<double> grid(minBound, maxBound, cellCount);

    IFFT_Synthesis fftHelper(N_FFT); 

    Eigen::MatrixXd spectrum(gridDimension3, N_FFT/2+1); 
    spectrum.setZero(); 

    Eigen::MatrixXd dataBuffer(gridDimension3, N_FFT); 
    bool initialized_dataBuffer = false; 

    int index = 0; 
    int count_spectrumComputed = 0; 
    Eigen::MatrixXd dataTimestep(gridDimension3,1); 
    Eigen::MatrixXd spectrumTimestep(gridDimension3,1); 
    while (index+N_FFT < N_timesteps)
    {

        // fill the buffer
        for (int ii=0; ii<N_FFT; ii++) 
        {
            std::cout << "reading time step : " << filenames[index+ii] << std::endl;
            {
                boost::timer::auto_cpu_timer processingTimer(" elapsed wall clock time : %w seconds\n"); 
                IO::readMatrixX<double>(dataTimestep, filenames[index+ii].c_str(), IO::BINARY, 0); 
                dataBuffer.col(ii) = dataTimestep.col(0); 
            }
        }
        index += N_FFT/2; // backup so there is overlap 


        // compute frequency content
        std::cout << "computing spectrum" << std::endl;
        {
            boost::timer::auto_cpu_timer timer_FFT(" elaspsed wall clock time : %w seconds"); 
            for (int ii=0; ii<dataBuffer.rows(); ii++) 
            {
                if (ii%100==0) std::cout << " progress : " << ii << "\r" << std::flush; 
                const Eigen::MatrixXd fftData = dataBuffer.row(ii).transpose(); 
                fftHelper.FFTspectrumAvgCol(fftData,spectrumTimestep); 
                spectrumTimestep.transposeInPlace(); 
                spectrum.row(ii) += spectrumTimestep; 
            }
            std::cout << endl; 
            count_spectrumComputed ++; 
        }

    }

    spectrum /= static_cast<double>(count_spectrumComputed); 

    for (int ii=0; ii<N_FFT/2+1; ii++) 
    {
        std::string filePrefix("HRIR_"+std::to_string(gridDimension)+"_spectrum_"+std::to_string(ii)); 

        std::shared_ptr<Eigen::MatrixXd> spectrumComponent = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd(spectrum.col(ii))); 
        grid.InsertCellCenteredData(filePrefix, spectrumComponent); 
        grid.WriteVTKCellCentered(IO::AssembleFilePath(dataDirectory,filePrefix+std::string(".vtk")), filePrefix, "spectrum"); 
    }

    IO::writeMatrixX<double>(spectrum, IO::AssembleFilePath(dataDirectory,"HRIR_"+std::to_string(gridDimension)+"_spectrum.dat").c_str(), IO::BINARY);

}


void SpectralAnalysis_OpenfoamQuadrupole()
{
    const std::string root("/hdd1/research_data/turbsound_data/OpenFOAM-run/velocity_dump/smooth_BL2_1m_refineconcave_y90"); 
    const std::string zone("FLUID_SMALL"); 

    const std::string dataDirectory = IO::AssembleFilePath(root,"out"); 
    IO::CreateDirectoryForce(dataDirectory); 

    OpenfoamCase foamCase(root,zone); 
    foamCase.ReinitializeCase(); 
    const double startTime = 1.80; 
    const double stopTime = -1;
    const int N_timesteps = foamCase.PrepareDataRead("Ug", 9, startTime, stopTime); 
    const int N_FFT = 256; 
    const int N_hops = N_timesteps/(N_FFT/2); 

    if (N_timesteps < N_FFT) 
        throw std::runtime_error("**WARNING** number of time steps data : "+std::to_string(N_timesteps)+"; requested NFFT : "+std::to_string(N_FFT)); 

    IFFT_Synthesis fftHelper(N_FFT); 

    Eigen::MatrixXd centroids; 
    foamCase.GetCentroids(centroids); 

    Eigen::MatrixXd dataBuffer; 
    Eigen::MatrixXd spectrumQuadrupole; 

    int N_cells = 0;
    int N_processors = 0; 
    // initialize
    {
        DataTimestep dataTimestep; 
        foamCase.ReadTimestep(1,dataTimestep); 
        N_processors = dataTimestep.processorIndex.size(); 
        for (int ii=0; ii<N_processors; ii++)
        {
            const int N_cellsProcessor = dataTimestep.fieldData[ii].rows()-foamCase.GetCellIndexOffsetProcessor(ii); 
            N_cells+=N_cellsProcessor; 
        }

        dataBuffer.resize(N_cells, N_FFT); 
        spectrumQuadrupole.resize(N_cells, N_FFT/2+1);
        spectrumQuadrupole.setZero(); 

    }
    std::cout << "total number of cells: " << N_cells << std::endl;

    int index = 0; 
    int count_spectrumComputed = 0; 
    while (index+N_FFT < N_timesteps)
    {
        // fill the buffer
        for (int ii=0; ii<N_FFT; ii++) 
        {
            boost::timer::auto_cpu_timer timer_read(" elapsed wall clock time : %w seconds\n"); 
            DataTimestep dataTimestep;

            foamCase.ReadTimestep(index+ii,dataTimestep); 

            std::string timestep = dataTimestep.timestep; 
            std::cout << "processing time step : " << timestep << std::endl;


            int fillIndex = 0; 
            for (int jj=0; jj<N_processors; jj++) 
            {
                const int cellIndexOffset = foamCase.GetCellIndexOffsetProcessor(jj); 
                const int dataRows = dataTimestep.fieldData[jj].rows(); 
                const int fillRows = dataRows-cellIndexOffset;
                dataBuffer.block(fillIndex,ii,fillRows,1) = 
                    dataTimestep.fieldData[jj].block(cellIndexOffset,0,fillRows,1).cwiseProduct(dataTimestep.fieldData[jj].block(cellIndexOffset,0,fillRows,1)) + 
                    dataTimestep.fieldData[jj].block(cellIndexOffset,4,fillRows,1).cwiseProduct(dataTimestep.fieldData[jj].block(cellIndexOffset,4,fillRows,1)) + 
                    dataTimestep.fieldData[jj].block(cellIndexOffset,8,fillRows,1).cwiseProduct(dataTimestep.fieldData[jj].block(cellIndexOffset,8,fillRows,1)) + 
                  ( dataTimestep.fieldData[jj].block(cellIndexOffset,1,fillRows,1).cwiseProduct(dataTimestep.fieldData[jj].block(cellIndexOffset,3,fillRows,1)) + 
                    dataTimestep.fieldData[jj].block(cellIndexOffset,2,fillRows,1).cwiseProduct(dataTimestep.fieldData[jj].block(cellIndexOffset,6,fillRows,1)) + 
                    dataTimestep.fieldData[jj].block(cellIndexOffset,5,fillRows,1).cwiseProduct(dataTimestep.fieldData[jj].block(cellIndexOffset,7,fillRows,1)) ) * 2.0; 

                fillIndex += fillRows; 
            }

        }
        index += N_FFT/2; // backup so there is overlap 


        // compute frequency content
        //


        int count_omp = 0; 
        std::cout << "computing spectrum" << std::endl;
        Eigen::MatrixXd fftBuffer(N_FFT,1); 
        //#pragma omp parallel for
        {
            boost::timer::auto_cpu_timer timer_FFT(" elaspsed wall clock time : %w seconds\n"); 

            for (int ii=0; ii<N_cells; ii++) 
            {
                const Eigen::MatrixXd fftData = dataBuffer.row(ii).transpose(); 
                Eigen::MatrixXd spectrumTimestep; 
                fftHelper.FFTspectrumAvgCol(fftData,spectrumTimestep); 
                spectrumTimestep.transposeInPlace(); 
                //#pragma omp critical
                {
                    count_omp++; 
                    //if ( count_omp%100 ==0 ) 
                    //    std::cout << " progress : " << (double)count_omp/dataBuffer.rows() << "\r" << std::flush; 
                    spectrumQuadrupole.row(ii) += spectrumTimestep; 
                }
            }
            std::cout << endl; 
            count_spectrumComputed ++; 
        }
    }

    spectrumQuadrupole /= static_cast<double>(count_spectrumComputed); 

    Eigen::VectorXd frequency = fftHelper.FFTFreq(10000.0); 
    
    for (int ii=0; ii<N_FFT/2+1; ii++) 
    {
        std::string filename("Tet_UgUg_spectrum_"+std::to_string(ii)+".csv"); 
        std::ofstream of(IO::AssembleFilePath(dataDirectory,filename).c_str()); 
        for (int jj=0; jj<N_cells; jj++) 
            of << centroids(jj,0) << "," << centroids(jj,1) << "," << centroids(jj,2) << "," << spectrumQuadrupole(jj,ii) << std::endl;

        of.close(); 
    }

    IO::writeMatrixX<double>(frequency, IO::AssembleFilePath(dataDirectory,"Tet_UgUg_spectrum_frequency.txt").c_str(), IO::ASCII);
}





/// perform out-of-core fft spectrum analysis 
void SpectralAnalysis_Quadrupole()
{
    const std::string root("/hdd1/research_data/turbsound_data/OpenFOAM-run/velocity_dump/smooth_BL2_1m_refineconcave_y90"); 
    const std::string zone("FLUID_SMALL"); 

    const std::string dataDirectory = IO::AssembleFilePath(root,"out"); 
    IO::CreateDirectoryForce(dataDirectory); 

    OpenfoamCase foamCase(root,zone); 
    foamCase.ReinitializeCase(); 
    const double startTime = 1.9; 
    const double stopTime = -1;
    const int N_timesteps = foamCase.PrepareDataRead("U", 3, startTime, stopTime); 
    const int N_FFT = 128; 
    const int N_hops = N_timesteps/(N_FFT/2); 

    if (N_timesteps < N_FFT) 
        throw std::runtime_error("**WARNING** number of time steps data : "+std::to_string(N_timesteps)+"; requested NFFT : "+std::to_string(N_FFT)); 

    // define grid for projection 
    const int gridDimension = 150;
    Eigen::Vector3d minBound(-0.2178337591671500,-0.2017749603171500,-0.2192579190471500); 
    Eigen::Vector3d maxBound(0.2204221281671500,0.2364809270171500,0.2189979682871500); 
    Eigen::Vector3i cellCount(gridDimension,gridDimension,gridDimension); 
    UniformGrid<double> grid(minBound, maxBound, cellCount);

    // compute projection given the openfoam mesh
    Eigen::MatrixXi meshGridProjectionTable; 
    std::string file_cachedTable("meshGridProjectionTabale_"+std::to_string(gridDimension)+".indicies");
    foamCase.GetMeshGridProjectionTable(grid, IO::AssembleFilePath(dataDirectory,file_cachedTable), meshGridProjectionTable, true); 
    IFFT_Synthesis fftHelper(N_FFT); 

    Eigen::MatrixXd spectrumQuadrupole; 
    bool initialized_spectrumQuadrupole = false; 

    Eigen::MatrixXd dataBuffer; 
    bool initialized_dataBuffer = false; 

    int index = 0; 
    int count_spectrumComputed = 0; 
    while (index+N_FFT < N_timesteps)
    {

        dataBuffer.setZero(); 

        // fill the buffer
        for (int ii=0; ii<N_FFT; ii++) 
        {
            std::shared_ptr<DataTimestep> dataTimestep(new DataTimestep());
            std::shared_ptr<Eigen::MatrixXd> projectedData(new Eigen::MatrixXd());; 

            foamCase.ReadTimestep(index+ii,*dataTimestep); 
            foamCase.ProjectDataTimestep(*dataTimestep, meshGridProjectionTable, grid, *projectedData); 

            const int N_cells          = projectedData->rows(); 
            const int N_dataDimensions = projectedData->cols();

            if (!initialized_spectrumQuadrupole) 
            {
                spectrumQuadrupole.resize(N_cells, N_FFT/2+1); 
                spectrumQuadrupole.setZero(); 
                initialized_spectrumQuadrupole = true; 
            }

            if (!initialized_dataBuffer) 
            {
                dataBuffer.resize(N_cells, N_FFT); 
                dataBuffer.setZero(); 
                initialized_dataBuffer = true; 
            }

            std::string timestep = dataTimestep->timestep; 
            std::cout << "reading time step : " << timestep << std::endl;
            {
                boost::timer::auto_cpu_timer processingTimer(" elapsed wall clock time : %w seconds\n"); 
                dataTimestep.reset();  // delete unused foam data

                std::shared_ptr<Eigen::MatrixXd> Ux = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd(projectedData->col(0))); 
                std::shared_ptr<Eigen::MatrixXd> Uy = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd(projectedData->col(1))); 
                std::shared_ptr<Eigen::MatrixXd> Uz = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd(projectedData->col(2))); 
                grid.InsertCellCenteredData(timestep+"_Ux", Ux); 
                grid.InsertCellCenteredData(timestep+"_Uy", Uy); 
                grid.InsertCellCenteredData(timestep+"_Uz", Uz); 

                //Eigen::VectorXd filter = SIGNAL_PROCESSING::DiscreteGaussian1D(5,2);
                Eigen::VectorXd filter = SIGNAL_PROCESSING::DiscreteGaussian1D(3,1);
                Ux.reset(new Eigen::MatrixXd(grid.CellCenteredSmoothing(timestep+"_Ux", filter)));
                Uy.reset(new Eigen::MatrixXd(grid.CellCenteredSmoothing(timestep+"_Uy", filter)));
                Uz.reset(new Eigen::MatrixXd(grid.CellCenteredSmoothing(timestep+"_Uz", filter)));

                //// reinserting the smoothed data for gradient
                grid.UpdateCellCenteredData(timestep+"_Ux", Ux);
                grid.UpdateCellCenteredData(timestep+"_Uy", Uy);
                grid.UpdateCellCenteredData(timestep+"_Uz", Uz);


                //// compute gradient on the grid
                Eigen::MatrixXd dUx = grid.CellCenteredDataGradient(timestep+"_Ux", UniformGrid<double>::ALL); 
                Eigen::MatrixXd dUy = grid.CellCenteredDataGradient(timestep+"_Uy", UniformGrid<double>::ALL); 
                Eigen::MatrixXd dUz = grid.CellCenteredDataGradient(timestep+"_Uz", UniformGrid<double>::ALL); 
                Eigen::MatrixXd Uxx = dUx.col(0); 
                Eigen::MatrixXd Uxy = dUx.col(1); 
                Eigen::MatrixXd Uxz = dUx.col(2); 
                Eigen::MatrixXd Uyx = dUy.col(0); 
                Eigen::MatrixXd Uyy = dUy.col(1); 
                Eigen::MatrixXd Uyz = dUy.col(2); 
                Eigen::MatrixXd Uzx = dUz.col(0); 
                Eigen::MatrixXd Uzy = dUz.col(1); 
                Eigen::MatrixXd Uzz = dUz.col(2); 

                grid.DeleteCellCenteredData(timestep+"_Ux"); 
                grid.DeleteCellCenteredData(timestep+"_Uy"); 
                grid.DeleteCellCenteredData(timestep+"_Uz"); 

                Eigen::MatrixXd quadrupoleTimestep = Uxx.cwiseProduct(Uxx) + Uyy.cwiseProduct(Uyy) + Uzz.cwiseProduct(Uzz) + 2.0*(Uxy.cwiseProduct(Uyx) + Uxz.cwiseProduct(Uzx) + Uyz.cwiseProduct(Uzy)); 

                //foamCase.WriteDataTimestepVTK(IO::AssembleFilePath(dataDirectory,"Uxx.vtk"), "Uxx", Uxx, grid);
                //foamCase.WriteDataTimestepVTK(IO::AssembleFilePath(dataDirectory,"Uxy_smoothed_5points.vtk"), "Uxy", d_Ux.col(1), grid);
                //foamCase.WriteDataTimestepVTK(IO::AssembleFilePath(dataDirectory,"Uxz_smoothed_5points.vtk"), "Uxz", d_Ux.col(2), grid);

                dataBuffer.col(ii) = quadrupoleTimestep.col(0); 
            }
        }
        index += N_FFT/2; // backup so there is overlap 


        // compute frequency content

        int count_omp = 0; 
        std::cout << "computing spectrum" << std::endl;
        //#pragma omp parallel for
        {
            boost::timer::auto_cpu_timer timer_FFT(" elaspsed wall clock time : %w seconds"); 
            for (int ii=0; ii<dataBuffer.rows(); ii++) 
            {
                const Eigen::MatrixXd fftData = dataBuffer.row(ii).transpose(); 
                Eigen::MatrixXd spectrumTimestep; 
                fftHelper.FFTspectrumAvgCol(fftData,spectrumTimestep); 
                spectrumTimestep.transposeInPlace(); 
                //#pragma omp critical
                {
                    count_omp++; 
                    //if ( count_omp%100 ==0 ) 
                    //    std::cout << " progress : " << (double)count_omp/dataBuffer.rows() << "\r" << std::flush; 
                    spectrumQuadrupole.row(ii) += spectrumTimestep; 
                }
            }
            std::cout << endl; 
            count_spectrumComputed ++; 
        }

    }

    spectrumQuadrupole /= static_cast<double>(count_spectrumComputed); 

    for (int ii=0; ii<N_FFT/2+1; ii++) 
    {
        std::string filename(std::to_string(gridDimension)+"_spectrum_"+std::to_string(ii)+".vtk"); 

        std::shared_ptr<Eigen::MatrixXd> spectrumComponent = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd(spectrumQuadrupole.col(ii))); 
        grid.InsertCellCenteredData(filename, spectrumComponent); 
        grid.WriteVTKCellCentered(IO::AssembleFilePath(dataDirectory,filename), filename, "spectrum"); 
    }

    std::string writeSpectrumFilename("spectrum_"+std::to_string(startTime)+"_"+std::to_string(stopTime)); 
    IO::writeMatrixX<double>(spectrumQuadrupole, IO::AssembleFilePath(dataDirectory,writeSpectrumFilename).c_str(), IO::BINARY);


}





int main()
{

    SpectralAnalysis_OpenfoamQuadrupole();
    //SpectralAnalysis_ImpulseResponse(); 
    //SpectralAnalysis_Quadrupole();
    return 0; 
}
