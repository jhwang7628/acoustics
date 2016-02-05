#include <iostream> 
#include <Eigen/Dense> 
#include "math/LeastSquareSurface.h" 

int main() 
{
    LeastSquareSurfaceBase *surface = new LeastSquareSurfaceLinear3D(); 

    const int N = 4; 

    const Eigen::MatrixXd samplePoints = Eigen::MatrixXd::Random(N,3); 
    const Eigen::VectorXd sampleValues = Eigen::VectorXd::Random(N); 

    std::cout << "test points = \n" << samplePoints << std::endl; 
    std::cout << "test values = \n" << sampleValues << std::endl; 
    surface->ComputeCoefficients(samplePoints, sampleValues); 

    Eigen::VectorXd interpolatedResults(N); 

    std::cout << "evaluating surface on test points ... " << std::flush; 
    bool passed = true; 
    int failedPoint = -1; 
    for (int ii=0; ii<N; ii++) 
    {
        interpolatedResults(ii) = surface->Interpolate(samplePoints.row(ii)); 
        std::cout << sampleValues[ii] << " -> " << interpolatedResults(ii) << std::endl;
        if (fabs(sampleValues[ii] - interpolatedResults(ii))>1E-12)
        {
            passed = false; 
            failedPoint = ii; 
            break; 
        }
    }

    if (!passed) 
    {
        std::cerr << "\n**WARNING** self testing on interpolation does not pass. the error is too high. check the matrix solve\n"; 
        std::cerr << "  failed point: input = " << sampleValues(failedPoint) << "; interpolatedResult = " << interpolatedResults(failedPoint) << std::endl;
    }
    else 
    {
        std::cout << "OK\n"; 
    }






}
