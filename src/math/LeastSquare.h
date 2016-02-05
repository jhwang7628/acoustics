#ifndef LEASTSQUARE_H 
#define LEASTSQUARE_H


class LeastSquareSurfaceBase
{
    public: 
        virtual void ComputeCoefficients() = 0; 
        virtual void Interpolate(const Eigen::Vector3d &evaluatePosition) = 0;
}


#endif
