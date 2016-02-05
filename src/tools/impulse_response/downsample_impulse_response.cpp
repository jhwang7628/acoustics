#include <iostream> 
#include <Eigen/Dense> 

#include "signalprocessing/resample.h" 
#include "utils/IO/IO.h" 


int main(int argv, char ** argc) 
{

    if (argv != 4) 
    {
        std::cerr << "**Usage " << argc[0] << " <in_HRTF_file> <out_HRTF_file> <skip>" << std::endl;
        exit(1); 
    }

    Eigen::MatrixXd HRTF; 
    std::cout << "reading " << std::endl;
    IO::readMatrixX<double>(HRTF, argc[1], IO::BINARY, 2);

    const int skipEvery = atoi(argc[3]);

    std::cout << "downsampling " << std::endl;
    SIGNAL_PROCESSING::NaiveDownsample(HRTF, skipEvery);

    std::cout << "writing " << std::endl;
    IO::writeMatrixX<double>(HRTF, argc[2], IO::BINARY); 

    return 0; 
}
