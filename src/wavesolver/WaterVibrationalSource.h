#ifndef WATER_VIBRATIONAL_SOURCE_H 
#define WATER_VIBRATIONAL_SOURCE_H 

#include <TYPES.h>
#include <wavesolver/VibrationalSource.h> 
#include <wavesolver/FDTD_RigidObject.h> 
#include <wavesolver/Wavesolver_ConstantsAndTypes.h> 
//#include <geometry/TriangleMesh.hpp>

//##############################################################################
// This class handles source evaluation for water surface objects.
//##############################################################################
class WaterVibrationalSource : public VibrationalSource
{
    private: 
        std::shared_ptr<TriangleMesh<REAL> >    _surfaceMesh; 
        FloatArray                              _oscillatorDisplacement;  
        FloatArray                              _oscillatorVelocity; 
        FloatArray                              _oscillatorAcceleration;

    public:
        WaterVibrationalSource(RigidObjectPtr owner);

        virtual REAL Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time); 
        virtual REAL EvaluateVelocity(const Vector3d &position, const Vector3d &normal, const REAL &time); 
        virtual REAL EvaluateDisplacement(const Vector3d &position, const Vector3d &normal, const REAL &time); 
        void Initialize(const std::string &wavFile); 
        void ReadOscillatorFromWav(const std::string &wavFile); 
        void ComputeVelocityAndAcceleration(); 
};

#endif
