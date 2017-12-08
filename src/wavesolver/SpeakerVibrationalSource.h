#ifndef SPEAKER_VIBRATIONAL_SOURCE_H 
#define SPEAKER_VIBRATIONAL_SOURCE_H 

#include <TYPES.h> 
#include <wavesolver/VibrationalSource.h> 
#include <wavesolver/FDTD_RigidObject.h> 

//##############################################################################
// Speaker vibrations. 
//##############################################################################
class SpeakerVibrationalSource : public VibrationalSource
{
    private: 

    public:
        SpeakerVibrationalSource(RigidObjectPtr owner) 
            : VibrationalSource(owner) 
        {}

        virtual REAL Evaluate(const Vector3d &position, 
                              const Vector3d &normal, 
                              const REAL &time, 
                              const int &hintTriangle=-1); 
        virtual REAL Evaluate(const int &vertexID, 
                              const Vector3d &vertexNormal, 
                              const REAL &time); 
        virtual Vector3d Evaluate(const int &vertexID, 
                                  const REAL &time);
        virtual REAL EvaluateVelocity(const Vector3d &position, 
                                      const Vector3d &normal, 
                                      const REAL &time); 
        virtual REAL EvaluateDisplacement(const Vector3d &position, 
                                          const Vector3d &normal, 
                                          const REAL &time); 
};

#endif
