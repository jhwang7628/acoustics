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
        struct ImpactRecord{Vector3d impactVector; REAL timestamp; int appliedVertex;}; 

    protected: 
        // object that owns this impulse series
        TriangleMeshPtr     _objectMesh; 

        // impulse data can be read from modalImpulse.txt
        int                     _lengthImpulses; 
        REAL                    _firstImpulseTime; 
        REAL                    _lastImpulseTime; 
        std::vector<Vector3d>   _impulses; 
        std::vector<REAL>       _impulseTimestamps; 
        std::vector<int>        _impulseAppliedVertex; 

        // rigidbody configurations can be read from config file
        RigidsimConfigDataPtr   _rigidsimConfigData; 

    public: 
        ImpulseSeriesObject(); 
        ImpulseSeriesObject(const TriangleMeshPtr &meshPtr); 

        inline Vector3d ConvertImpulseToForce(const Vector3d &impulse){return impulse / _rigidsimConfigData->simulation_step_size; }
        inline int Size(){return _lengthImpulses;}
        inline int N_Impulses(){return _lengthImpulses;}
        inline REAL GetFirstImpulseTime(){return _firstImpulseTime;}
        inline void SetMesh(const TriangleMeshPtr &meshPtr){_objectMesh = meshPtr;}
        inline bool Initialized(){return _objectMesh!=nullptr && _lengthImpulses!=0;}
        inline RigidsimConfigDataPtr GetRigidsimConfigData(){return _rigidsimConfigData;}
        inline void SetRigidsimConfigData(RigidsimConfigDataPtr configData){_rigidsimConfigData = configData;}
        inline REAL GetRigidsimTimeStepSize(){return _rigidsimConfigData->simulation_step_size;}
        void Initialize(); 
        // Add impulse for this object. The vertex index should be for surface
        // triangle mesh (not volumetric tetrahedron). 
        void AddImpulse(const REAL &timestamp, const int &appliedVertex, const Vector3d &impulse); 
        // Get single impulse frame by index
        void GetImpulse(const int &index, REAL &timestamp, int &vertex, Vector3d &impulse); 
        // Get all impulses that have timestamps within [timeStart, timeStop]
        void GetImpulse(const REAL &timeStart, const REAL &timeStop, std::vector<ImpactRecord> &records); 
        // Get all forces (scaled impulses) that have timestamps [timeStart, timeStop]
        void GetForces(const REAL &timeStart, const REAL &timeStop, std::vector<ImpactRecord> &records); 
        void GetImpulseRange(REAL &firstImpulseTime, REAL &lastImpulseTime);
};

#endif 
