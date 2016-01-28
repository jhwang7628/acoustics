#ifndef FILTERDESIGN_H
#define FILTERDESIGN_H

#include <iostream> 
#include <Eigen/Dense>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288419716939937510582
#endif

namespace SIGNAL_PROCESSING
{

Eigen::VectorXd DiscreteGaussian1D(const int &N, const double &stddev) 
{
  
    double center = static_cast<double>(N-1)/2.0;
    Eigen::VectorXd gaussian = Eigen::VectorXd::LinSpaced(N,0,N-1);
    gaussian.array() -= center;
    gaussian /= stddev; 
    gaussian = gaussian.array().square();
    gaussian /= -2.0;
    gaussian = gaussian.array().exp();
    gaussian /= stddev*sqrt(2.0*M_PI);

    return gaussian; 



}



};

#endif
