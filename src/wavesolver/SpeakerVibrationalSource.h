#ifndef SPEAKER_VIBRATIONAL_SOURCE_H 
#define SPEAKER_VIBRATIONAL_SOURCE_H 

#include <unordered_set>
#include <TYPES.h> 
#include <wavesolver/VibrationalSource.h> 
#include <wavesolver/FDTD_RigidObject.h> 

//##############################################################################
// Speaker vibrations. 
//##############################################################################
class SpeakerVibrationalSource : public VibrationalSource
{
    private: 
        std::vector<REAL> _speakerData;
        REAL              _speakerDataSampleRate; 
        std::unordered_set<int> _handles; 

    public:
        SpeakerVibrationalSource(RigidObjectPtr owner) 
            : VibrationalSource(owner) 
        {}

        void Initialize(const std::string &speakerFile, 
                        const std::vector<int> &handleVIds);
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
