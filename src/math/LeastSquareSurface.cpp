#include <iostream> 
#include "LeastSquareSurface.h" 

void LeastSquareSurfaceLinear3D::ComputeCoefficients(const Eigen::MatrixXd &samplePoints, const Eigen::VectorXd &sampleValues)
{
    std::cout << "computing coefficients\n";

    if (samplePoints.cols()!=3) throw std::runtime_error("**ERROR** sample points is not three dimensional. input matrix number of columns : " + std::to_string(samplePoints.cols())); 
    if (samplePoints.rows()!=sampleValues.size()) throw std::runtime_error("**ERROR** sample size mismatch. sample points number : " + std::to_string(samplePoints.rows()) + "; sample values number : " + std::to_string(sampleValues.size())); 

    const int N_points = samplePoints.rows(); 

    Eigen::MatrixXd A(N_points, 4);  // linear surface 

    A.col(0) = Eigen::VectorXd::Ones(N_points); 
    A.col(1) = samplePoints.col(0); 
    A.col(2) = samplePoints.col(1); 
    A.col(3) = samplePoints.col(2); 

    coefficients = A.householderQr().solve(sampleValues); 

}

double LeastSquareSurfaceLinear3D::Interpolate(const Eigen::Vector3d &evaluatePosition)
{
    return coefficients[0] + coefficients[1]*evaluatePosition[0] + coefficients[2]*evaluatePosition[1] + coefficients[3]*evaluatePosition[2]; 
}/
