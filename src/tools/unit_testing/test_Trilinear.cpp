#include <iostream> 

#include <distancefield/trilinearInterpolation.h>
#include <tools/unit_testing/testing.h> 

#include <Eigen/Dense> 

void FillVandermondeMatrix(const int &row, const Eigen::Vector3d &point, Eigen::Vector3d offset, Eigen::MatrixXd &V) 
{

    const double x = point[0] + offset[0]; 
    const double y = point[1] + offset[1]; 
    const double z = point[2] + offset[2]; 
    V(row, 0) = x*y*z; 
    V(row, 1) = x*y; 
    V(row, 2) = x*z; 
    V(row, 3) = y*z; 
    V(row, 4) = x; 
    V(row, 5) = y; 
    V(row, 6) = z; 
    V(row, 7) = 1; 

}


void FillVandermondeBoundary(const int &row, const Eigen::Vector3d &boundaryPosition, const Eigen::Vector3d &boundaryNormal, Eigen::MatrixXd &V)
{
    const double &x = boundaryPosition[0];    
    const double &y = boundaryPosition[1]; 
    const double &z = boundaryPosition[2]; 
    
    const double &nx= boundaryNormal[0];    
    const double &ny= boundaryNormal[1];
    const double &nz= boundaryNormal[2];
    
    V(row,0) = nx*y*z + ny*x*z + nz*x*y;                                                                                                                      
    V(row,1) = nx*y + ny*x; 
    V(row,2) = nx*z + nz*x;                                                                                                                                   
    V(row,3) = ny*z + nz*y;                                                                                                                                   
    V(row,4) = nx;  
    V(row,5) = ny;                                                                                                                                            
    V(row,6) = nz;                                                                                                                                            
    V(row,7) = 0;
}


int main()
{

    const double v000 = 1; 
    const double v001 = 2; 
    const double v010 = 3; 
    const double v011 = 4; 
    const double v100 = 5; 
    const double v101 = 6; 
    const double v110 = 7; 
    const double v111 = 8; 




    const Eigen::Vector3d low(0.85, 0.85, 0.85); 
    const double h = 0.01;
    const Eigen::Vector3d high = low.array() + h; 

    const Eigen::Vector3d x = low.array() + h*0.9; 
    const Eigen::Vector3d boundaryPoint = low.array() + h*0.3; 
    const Eigen::Vector3d boundaryNormal(1,2,3); 

    const Eigen::Vector3d wx = (x - low).array()/(high-low).array(); 

    const double answer = TRILINEAR_INTERPOLATION(wx[0],wx[1],wx[2],v000,v100,v110,v010,v001,v101,v111,v011); 


    Eigen::MatrixXd V(8,8); 
    Eigen::VectorXd phiSample(8);
    FillVandermondeMatrix(0, low, Eigen::Vector3d(0,0,0), V); 
    FillVandermondeMatrix(1, low, Eigen::Vector3d(0,0,h), V); 
    FillVandermondeMatrix(2, low, Eigen::Vector3d(0,h,0), V); 
    FillVandermondeMatrix(3, low, Eigen::Vector3d(0,h,h), V); 
    FillVandermondeMatrix(4, low, Eigen::Vector3d(h,0,0), V); 
    FillVandermondeMatrix(5, low, Eigen::Vector3d(h,0,h), V); 
    FillVandermondeMatrix(6, low, Eigen::Vector3d(h,h,0), V); 
    //FillVandermondeMatrix(7, low, Eigen::Vector3d(h,h,h), V); 
    FillVandermondeBoundary(7, boundaryPoint, boundaryNormal, V); 

    phiSample(0) = v000; 
    phiSample(1) = v001; 
    phiSample(2) = v010; 
    phiSample(3) = v011; 
    phiSample(4) = v100; 
    phiSample(5) = v101; 
    phiSample(6) = v110; 
    //phiSample(7) = v111; 
    phiSample(7) = 1000;  // neumann

    Eigen::MatrixXd b(1,8); 
    FillVandermondeMatrix(0, x, Eigen::Vector3d(0,0,0), b); 

    Eigen::JacobiSVD<Eigen::MatrixXd> svd_V(V, Eigen::ComputeThinU | Eigen::ComputeThinV); 
    Eigen::VectorXd C = svd_V.solve(phiSample);

    //b.transposeInPlace();
    b.transposeInPlace(); 
    V.transposeInPlace();

    Eigen::JacobiSVD<Eigen::MatrixXd> svd_VT(V, Eigen::ComputeThinU | Eigen::ComputeThinV);

    Eigen::VectorXd beta = svd_VT.solve(b);

    COUT_SDUMP(beta);
    const double conditionNumber = svd_VT.singularValues()(0) / svd_VT.singularValues()(svd_VT.singularValues().size()-1);
    //if (conditionNumber > 1E5)
    {
        std::cout << "**WARNING** condition number for the least square solve is = " << conditionNumber << "\n"
            << "            the solution can be inaccurate.\n"
            << "largest  singular value = " << svd_VT.singularValues()(0) << "\n"
            << "smallest singular value = " << svd_VT.singularValues()(svd_VT.singularValues().size()-1) << "\n"
            << "problematic matrix: \n"
            << V << std::endl;
    }

    COUT_SDUMP(beta);
    COUT_SDUMP(beta.dot(phiSample));

    Eigen::VectorXd B = b; 
    std::cout << "C^T B : " << C.dot(B) << std::endl;
    std::cout << "beta^T phi_sample : " << beta.dot(phiSample) << std::endl;


    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > testPoints(8); 

    testPoints[0] = low.array()+ Eigen::Vector3d(0,0,0).array(); 
    testPoints[1] = low.array()+ Eigen::Vector3d(0,0,h).array();
    testPoints[2] = low.array()+ Eigen::Vector3d(0,h,0).array();
    testPoints[3] = low.array()+ Eigen::Vector3d(0,h,h).array();
    testPoints[4] = low.array()+ Eigen::Vector3d(h,0,0).array();
    testPoints[5] = low.array()+ Eigen::Vector3d(h,0,h).array();
    testPoints[6] = low.array()+ Eigen::Vector3d(h,h,0).array();


    for (int ii=0; ii<8; ii++) 
    {
        if (ii<7)
        {
            std::cout << "testPoint " << ii << " : " << testPoints[ii].transpose() << std::endl; 
            Eigen::MatrixXd b_testPoint(1,8); 
            Eigen::Vector3d tmp(0,0,0);
            FillVandermondeMatrix(0, testPoints[ii], tmp, b_testPoint); 

            b_testPoint.transposeInPlace();

            Eigen::VectorXd B_testPoint = b_testPoint; 

            std::cout << " C^T B_testPoint = " << C.dot(B_testPoint) << std::endl;
        }
        else 
        {
            std::cout << "testPoint " << ii << " : " << boundaryPoint.transpose() << std::endl; 
            Eigen::MatrixXd b_testPoint(1,8); 
            FillVandermondeBoundary(0, boundaryPoint, boundaryNormal, b_testPoint); 

            b_testPoint.transposeInPlace();

            Eigen::VectorXd B_testPoint = b_testPoint; 

            std::cout << " B_testPoint = " << B_testPoint << std::endl;
            std::cout << " C^T B_testPoint = " << C.dot(B_testPoint) << std::endl;
        }
    }

    


    return 0;
}
