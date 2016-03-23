#ifndef LEASTSQUARESURFACE_H 
#define LEASTSQUARESURFACE_H

#include <Eigen/Dense> 

// base class for least square interpolation surface 
class LeastSquareSurfaceBase
{
    public: 
        virtual void ComputeCoefficients(const Eigen::MatrixXd &samplePoints, const Eigen::VectorXd &sampleValues) = 0; 
        virtual double Evaluate(const Eigen::Vector3d &evaluatePosition) = 0;
};

// compute a least square given n points. The least square has four unknowns,
// a_0, a_1, a_2, a_3. After coefficients being computed, can use this surface
// to evaluate at any given point in R3. 
//
// phi_1 = a_0 + a_1 x_1 + a_2 y_1 + a_3 z_1 
// phi_2 = a_0 + a_1 x_2 + a_2 y_2 + a_3 z_2 
// phi_3 = a_0 + a_1 x_3 + a_2 y_3 + a_3 z_3 
//                  .
//                  .
//                  .
// phi_n = a_0 + a_1 x_n + a_2 y_n + a_3 z_n 
//
// object can be reused, there is no critical initialization involved. 
//
class LeastSquareSurfaceLinear3D : public LeastSquareSurfaceBase
{

    public: 
        enum WeightingType{ UNIFORM, GAUSSIAN, FLEISHMAN }; 

    private: 
        Eigen::Vector4d coefficients; // the coefficients a_0, a_1, a_2, a_3
        Eigen::VectorXd weights; 

    public: 
        virtual void ComputeCoefficients(const Eigen::MatrixXd &samplePoints, const Eigen::VectorXd &sampleValues); 
        virtual double Evaluate(const Eigen::Vector3d &evaluatePosition); 
        
};




#endif
