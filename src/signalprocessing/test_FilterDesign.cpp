#include "FilterDesign.h" 
#include <Eigen/Dense>
#include <iostream> 
#include <fstream> 

int main()
{
    Eigen::VectorXd tmp = SIGNAL_PROCESSING::DiscreteGaussian1D(3,1);

    std::ofstream of("out.txt"); 
    of << tmp << std::endl;

    return 0; 
}
