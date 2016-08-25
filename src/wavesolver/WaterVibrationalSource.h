#ifndef WATER_VIBRATIONAL_SOURCE_H 
#define WATER_VIBRATIONAL_SOURCE_H 

#include <TYPES.h>
#include <interp/CSpline.hpp>
#include <wavesolver/VibrationalSource.h> 
#include <wavesolver/FDTD_RigidObject.h> 
#include <wavesolver/Wavesolver_ConstantsAndTypes.h> 

//##############################################################################
// This class handles source evaluation for water surface objects.
//##############################################################################
class WaterVibrationalSource : public VibrationalSource
{
    public: 
        typedef CSpline<REAL, false> T_CSpline; 
    private: 
        std::shared_ptr<TriangleMesh<REAL> >    _surfaceMesh; 
        REAL                                    _sampleRate; 
        REAL                                    _startTime = 0.0; 
        FloatArray                              _oscillatorTime; 
        FloatArray                              _oscillatorDisplacement;  
        FloatArray                              _oscillatorVelocity; 
        FloatArray                              _oscillatorAcceleration;
        T_CSpline                               _interpolatorVelocity; 
        T_CSpline                               _interpolatorAcceleration; 

    public:
        WaterVibrationalSource(RigidObjectPtr owner, const std::string &wavFile);

        virtual REAL Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time); 
        virtual REAL EvaluateVelocity(const Vector3d &position, const Vector3d &normal, const REAL &time); 
        virtual REAL EvaluateDisplacement(const Vector3d &position, const Vector3d &normal, const REAL &time); 
        void Initialize(const std::string &wavFile); 
        void ReadOscillatorFromWav(const std::string &wavFile); 
        void ComputeVelocityAndAcceleration(); 
        void PrecomputeInterpolation(); 

        ///// debug functions /////
        void TestSpline(); 
};

#endif
