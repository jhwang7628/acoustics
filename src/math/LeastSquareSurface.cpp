#include <iostream> 
#include "LeastSquareSurface.h" 


void LeastSquareSurfaceLinear3D::ComputeCoefficients(const Eigen::MatrixXd &samplePoints, const Eigen::VectorXd &sampleValues)
{

    if (samplePoints.cols()!=3) throw std::runtime_error("**ERROR** sample points is not three dimensional. input matrix number of columns : " + std::to_string(samplePoints.cols())); 
    if (samplePoints.rows()!=sampleValues.size()) throw std::runtime_error("**ERROR** sample size mismatch. sample points number : " + std::to_string(samplePoints.rows()) + "; sample values number : " + std::to_string(sampleValues.size())); 

    const int N_points = samplePoints.rows(); 

    Eigen::MatrixXd A(N_points, 4);  // linear surface 

    A.col(0) = Eigen::VectorXd::Ones(N_points); 
    A.col(1) = samplePoints.col(0); 
    A.col(2) = samplePoints.col(1); 
    A.col(3) = samplePoints.col(2); 

    // For solving the least square system on A=mxn matrix, QR is a little bit faster
    // 2mn^2-2/3n^3 flops, for SVD its 2mn^2+11n^3 (pg 84 Trefethen Bau)
    // to 4mn^2-4/3n^3 (lecture 31 TB). However, SVD has the advantage of
    // finding the minimum norm solution to the least square problems (for
    // rank-deficient matrices the solution is not unique) and is
    // more stable in that it perserves the 2-norm. 
    //
    // quote from ref [2]: The power of the SVD lies in the fact that it always exists and can be
    // computed stably. The computed SVD will be well-conditioned because orthogonal
    // matrices preserve the 2-norm. Any perturbation in A will not be amplified by the SVD. 
    //
    // [1] Trefethen and Bau, Numerilcal Linear Algebra
    // [2] D.Q. Lee, numerically efficient methods for solving least squares problems
    // [3] good reference on column-pivoted QR:
    //      http://www.math.usm.edu/lambers/mat610/sum10/lecture11.pdf
    //
    // Note that using the svd solve, the solution is implicitly a least square
    // solution, thin SVD is enough. 
    // (http://eigen.tuxfamily.org/dox/classEigen_1_1JacobiSVD.html)
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    _coefficients = svd.solve(sampleValues); 
    const double conditionNumber = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
    if (conditionNumber > 1E5) 
    {
        std::cout << "**WARNING** condition number for the least square solve is = " << conditionNumber << "\n"; 
        std::cout << "            the solution can be inaccurate.\n";
        std::cout << "            the problematic matrix: \n"; 
        std::cout << samplePoints << std::endl;
    }

    //coefficients = A.householderQr().solve(sampleValues); 
    



}

double LeastSquareSurfaceLinear3D::Evaluate(const Eigen::Vector3d &evaluatePosition)
{
    return _coefficients[0] + _coefficients[1]*evaluatePosition[0] + _coefficients[2]*evaluatePosition[1] + _coefficients[3]*evaluatePosition[2]; 
}



double WeightedLeastSquareSurfaceLinear3D::Evaluate(const Eigen::MatrixXd &samplePoints, const Eigen::VectorXd &sampleValues, const Eigen::Vector3d &evaluatePosition)
{
    if (samplePoints.cols()!=3) throw std::runtime_error("**ERROR** sample points is not three dimensional. input matrix number of columns : " + std::to_string(samplePoints.cols())); 
    if (samplePoints.rows()!=sampleValues.size()) throw std::runtime_error("**ERROR** sample size mismatch. sample points number : " + std::to_string(samplePoints.rows()) + "; sample values number : " + std::to_string(sampleValues.size())); 

    const int N_points = samplePoints.rows(); 

    Eigen::MatrixXd A(N_points, 4);  // linear surface 

    A.col(0) = Eigen::VectorXd::Ones(N_points); 
    A.col(1) = samplePoints.col(0); 
    A.col(2) = samplePoints.col(1); 
    A.col(3) = samplePoints.col(2); 


    // compute weights to rows of A as well as RHS of the equation
    if (!_setWeights)
        ComputeWeights(samplePoints, evaluatePosition); 

    //std::cout << "A before: " << std::endl; 
    //std::cout << A << std::endl;
    A.array().colwise() *= _weights.array();
    //std::cout << "A after: " << std::endl; 
    //std::cout << A << std::endl;

    Eigen::VectorXd weightedSampleValues = sampleValues.array()*_weights.array(); 

    // for choice of linear solver please see leastsquaresurface
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    _coefficients = svd.solve(weightedSampleValues); 
    const double conditionNumber = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
    if (conditionNumber > 1E5) 
    {
        std::cout << "**WARNING** condition number for the least square solve is = " << conditionNumber << "\n"; 
        std::cout << "            the solution can be inaccurate.\n";
        std::cout << "problematic points: \n"; 
        std::cout << samplePoints << std::endl;
        std::cout << "weights: \n"; 
        std::cout << _weights.transpose() << std::endl;
    }

    //coefficients = A.householderQr().solve(sampleValues); 

    return _coefficients[0] + _coefficients[1]*evaluatePosition[0] + _coefficients[2]*evaluatePosition[1] + _coefficients[3]*evaluatePosition[2]; 
}

void WeightedLeastSquareSurfaceLinear3D::ComputeWeights(const Eigen::MatrixXd &samplePoints, const Eigen::Vector3d &evaluatePosition)
{
    _weights.resize(samplePoints.rows()); 
    double distanceBuffer; 

    switch (_weightingType)
    {
        case GAUSSIAN: 
            for (int ii=0; ii<samplePoints.rows(); ii++)
            {
                distanceBuffer  = pow(samplePoints(ii,0)-evaluatePosition(0),2) 
                                + pow(samplePoints(ii,1)-evaluatePosition(1),2) 
                                + pow(samplePoints(ii,2)-evaluatePosition(2),2); 
                _weights(ii) = exp(-distanceBuffer/_gaussianSpacingSquared); 
            }
            break; 

        case INVERSE_QUADRATIC: 
            for (int ii=0; ii<samplePoints.rows(); ii++)
            {
                distanceBuffer  = pow(samplePoints(ii,0)-evaluatePosition(0),2) 
                                + pow(samplePoints(ii,1)-evaluatePosition(1),2) 
                                + pow(samplePoints(ii,2)-evaluatePosition(2),2); 
                _weights(ii) = 1./(distanceBuffer + _inverseQuadraticEpsSquared); 
            }
            break; 
    }

    _weights.normalize();
}



