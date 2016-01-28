#include "FilterDesign.h" 
#include <Eigen/Dense>
#include <iostream> 

int main()
{
    Eigen::VectorXd tmp = SIGNAL_PROCESSING::DiscreteGaussian1D(11,4);
    std::cout << tmp.transpose() << std::endl;

    return 0; 
}
