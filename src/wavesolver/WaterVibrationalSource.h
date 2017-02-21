#ifndef WATER_VIBRATIONAL_SOURCE_H 
#define WATER_VIBRATIONAL_SOURCE_H 

#include <TYPES.h>
#include <interp/BasicInterp.hpp>
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
        typedef std::shared_ptr<TriangleMesh<REAL> > TriangleMeshPtr;
        typedef CSpline<REAL, false, FixAcc<REAL> > T_CSpline; 
        struct GaussianDecayModel
        {
            Vector3d    center; 
            REAL        stddev; 
        }; 

    private: 
        Vector3d            _wantedNormal = Vector3d(0, 1, 0);
        REAL                _validAngleThreshold = 0.5; // use to determine if the vertex has source. See Evaluate() for usage. default: 0.5
        TriangleMeshPtr     _surfaceMesh; 
        REAL                _sampleRate; 
        REAL                _startTime = 0.0; 
        FloatArray          _oscillatorTime; 
        FloatArray          _oscillatorDisplacement;  
        FloatArray          _oscillatorVelocity; 
        FloatArray          _oscillatorAcceleration;
        T_CSpline           _interpolatorVelocity; 
        T_CSpline           _interpolatorAcceleration; 
        GaussianDecayModel  _decayModel; 
        REAL                _decayRadius;

    public:
        WaterVibrationalSource(RigidObjectPtr owner, const std::string &wavFile, const REAL &decayRadius);

        virtual REAL Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time, const int &hintTriangle); 
        virtual REAL Evaluate(const int &vertexID, const Vector3d &vertexNormal, const REAL &time); 
        virtual REAL EvaluateVelocity(const Vector3d &position, const Vector3d &normal, const REAL &time); 
        virtual REAL EvaluateDisplacement(const Vector3d &position, const Vector3d &normal, const REAL &time); 
        inline REAL Decay(const Vector3d &samplePoint){
            return exp(-(_owner->WorldToObjectPoint(samplePoint) - _decayModel.center).lengthSqr() / (2.0 * pow(_decayModel.stddev,2)));} 
        void Initialize(const std::string &wavFile); 
        void InitializeDecayModel(); 
        void ReadOscillatorFromWav(const std::string &wavFile); 
        void ComputeVelocityAndAcceleration(); 
        void PrecomputeInterpolation(); 

        ///// debug functions /////
        void TestSpline(); 
};

#endif
