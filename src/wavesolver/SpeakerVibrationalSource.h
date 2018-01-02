#ifndef SPEAKER_VIBRATIONAL_SOURCE_H 
#define SPEAKER_VIBRATIONAL_SOURCE_H 
#include <queue>
#include <unordered_set>
#include <TYPES.h> 
#include <wavesolver/VibrationalSource.h> 
#include <wavesolver/FDTD_RigidObject.h> 

//##############################################################################
// Speaker vibrations. 
//##############################################################################
class SpeakerVibrationalSource : public VibrationalSource
{
    public: 
        struct DataStep
        {
            int frame; 
            std::string objFilePrefix; 
            std::vector<int> handles; 
        };
    private: 
        std::string       _speakerWavFile; 
        std::vector<REAL> _speakerData;
        std::vector<char> _speakerDataMask;
        REAL              _speakerDataSampleRate; 
        REAL              _speakerDataStartTime = 0.0; 
        std::unordered_set<int> _handles; 
        std::queue<DataStep> _objSeqData;
        REAL                 _objSeqSampleRate; 
        // keep this low, user should provide proper mask
        const REAL        _speakerDataThreshold = 1E-15;

    public:
        SpeakerVibrationalSource() = default;
        SpeakerVibrationalSource(RigidObjectPtr owner) 
            : VibrationalSource(owner) 
        {}

        DataStep ReadObjSeqMetaData(const std::string &dir, const std::string &objPrefix,
                                const std::string &speakerVIdsDir, 
                                const std::string &speakerVIdsSuf);
        void SetSeqSampleRate(const REAL &sampleRate){_objSeqSampleRate = sampleRate;}
        void Initialize(const std::string &speakerFile, 
                        const std::vector<int> &handleVIds);
        virtual bool UpdateTime(const REAL time); 
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
        virtual REAL EarliestEventTime(const REAL startTime) const; 
        ///// debug /////
        void PrintHandlesTotalArea(); 
};

#endif
