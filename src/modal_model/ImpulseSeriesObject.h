#ifndef IMPULSE_SERIES_OBJECT_H 
#define IMPULSE_SERIES_OBJECT_H 
#include <string> 
#include <vector> 
#include <memory>
#include <config.h>
#include <linearalgebra/Vector3.hpp>
#include <geometry/TriangleMesh.hpp>
#include <rigid/RigidsimConfigData.h>

//##############################################################################
// Class that stores the impulse information for an object.
//
// Can be created by the IO class io/ImpulseSeriesReader
//##############################################################################
class ImpulseSeriesObject
{
    public: 
        typedef std::shared_ptr<TriangleMesh<REAL> > TriangleMeshPtr; 
        typedef std::shared_ptr<RigidsimConfigData> RigidsimConfigDataPtr; 
        // impactVector could have been impulse or force. force will be computed by
        // member function
        struct ImpactRecord
        {
            Vector3d impactVector; 
            Vector3d impactPosition; // in object space
            REAL timestamp; 
            REAL supportLength = 0.0;  // tau
            REAL contactSpeed = 0.0; 
            REAL gamma;                // gamma = pi * norm(J) / (2 tau)
            int appliedVertex;
        }; 

    protected: 
        // object that owns this impulse series
        TriangleMeshPtr     _objectMesh; 

        // impulse data can be read from modalImpulse.txt
        REAL                    _firstImpulseTime; 
        REAL                    _lastImpulseTime; 
        std::vector<ImpactRecord> _impulses; 

        // rigidbody configurations can be read from config file
        RigidsimConfigDataPtr   _rigidsimConfigData; 

    public: 
        ImpulseSeriesObject(); 
        ImpulseSeriesObject(const TriangleMeshPtr &meshPtr); 

        inline int N_Impulses(){return _impulses.size();}
        inline REAL GetFirstImpulseTime(){return _firstImpulseTime;}
        inline void SetMesh(const TriangleMeshPtr &meshPtr){_objectMesh = meshPtr;}
        inline bool Initialized(){return _objectMesh!=nullptr && N_Impulses()!=0;}
        inline RigidsimConfigDataPtr GetRigidsimConfigData(){assert(_rigidsimConfigData); return _rigidsimConfigData;}
        inline void SetRigidsimConfigData(RigidsimConfigDataPtr configData){_rigidsimConfigData = configData;}
        inline REAL GetRigidsimTimeStepSize(){assert(_rigidsimConfigData); return _rigidsimConfigData->simulation_step_size;}
        void Initialize(); 
        void AddImpulse(const ImpactRecord &record); 
        void GetImpulse(const int &index, REAL &timestamp, int &vertex, Vector3d &impulse); 
        void GetImpulse(const REAL &timeStart, const REAL &timeStop, std::vector<ImpactRecord> &records); 
        void GetImpulseWithinSupport(const REAL &timeStart, std::vector<ImpactRecord> &records); 
        void GetForces(const REAL &time, std::vector<ImpactRecord> &records); 
        void GetRangeOfImpulses(REAL &firstImpulseTime, REAL &lastImpulseTime);
        void Filter(); 

        //// debug methods ////
        void PrintAllImpulses(); 

    friend std::ostream &operator <<(std::ostream &os, const ImpactRecord &record);
};

#endif 
