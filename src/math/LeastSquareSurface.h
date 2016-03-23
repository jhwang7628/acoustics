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

    private: 
        Eigen::Vector4d _coefficients; // the coefficients a_0, a_1, a_2, a_3

    public: 
        virtual void ComputeCoefficients(const Eigen::MatrixXd &samplePoints, const Eigen::VectorXd &sampleValues); 
        virtual double Evaluate(const Eigen::Vector3d &evaluatePosition); 

}; 


// base class for least square interpolation surface 
class WeightedLeastSquareSurfaceBase
{
    public: 
        virtual double Evaluate(const Eigen::MatrixXd &samplePoints, const Eigen::VectorXd &sampleValues, const Eigen::Vector3d &evaluatePosition) = 0;
};
        
// similar to LeastSquareSurfaceLinear3D but all points are weighted by
// a weighting function. coefficients depends on weights, which in turn
// depends on the evauluating position.
//
// See ref. ASAP Intro to LS, WLS, MLS by Andrew Nealen for more info
class WeightedLeastSquareSurfaceLinear3D : public WeightedLeastSquareSurfaceBase
{
    public: 
        enum WeightingType{ GAUSSIAN, INVERSE_QUADRATIC }; 

    private: 
        WeightingType   _weightingType; 
        Eigen::Vector4d _coefficients; // the coefficients a_0, a_1, a_2, a_3
        Eigen::VectorXd _weights; 

        // for these parameters please see eq(9) and eq(11) in Nealen's paper
        // essentially, 
        //  gaussian spacing is for smoothing out small features in the function
        //  inverse quadratic epsilon can be viewed as how much interpolating we want our weighted least square to be
        //
        double          _gaussianSpacing; 
        double          _gaussianSpacingSquared; 
        double          _inverseQuadraticEps; 
        double          _inverseQuadraticEpsSquared; 

        bool            _setWeights; 

    public: 
        virtual double Evaluate(const Eigen::MatrixXd &samplePoints, const Eigen::VectorXd &sampleValues, const Eigen::Vector3d &evaluatePosition);

        // set the weights
        void ComputeWeights(const Eigen::MatrixXd &samplePoints, const Eigen::Vector3d &evaluatePosition); 

        void SetWeights(const Eigen::VectorXd &weights) { _weights = weights; _setWeights = true; } 

        inline Eigen::VectorXd GetWeights(){ return _weights; } 
        inline void SetGaussianSpacing(const double &spacing)
        { 
            _gaussianSpacing = spacing; 
            _gaussianSpacingSquared = _gaussianSpacing*_gaussianSpacing; 
        } 
        inline void SetInverseQuadraticEps(const double &e)
        {
            _inverseQuadraticEps = e; 
            _inverseQuadraticEps = _inverseQuadraticEps*_inverseQuadraticEps; 
        }



        WeightedLeastSquareSurfaceLinear3D(): 
            _weightingType(GAUSSIAN), 
            _gaussianSpacing(0.005), 
            _gaussianSpacingSquared(_gaussianSpacing*_gaussianSpacing),
            _inverseQuadraticEps(1E-6), 
            _inverseQuadraticEpsSquared(_inverseQuadraticEps*_inverseQuadraticEps), 
            _setWeights(false)
        {
        }
};




#endif
