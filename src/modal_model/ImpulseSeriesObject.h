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

        inline int Size(){return _lengthImpulses;}
        inline void SetMesh(const TriangleMeshPtr &meshPtr){_objectMesh = meshPtr;}
        inline bool Initialized(){return _objectMesh!=nullptr && _lengthImpulses!=0;}
        inline RigidsimConfigDataPtr GetRigidsimConfigData(){return _rigidsimConfigData;}
        inline void SetRigidsimConfigData(RigidsimConfigDataPtr configData){_rigidsimConfigData = configData;}
        void Initialize(); 
        // Add impulse for this object. The vertex index should be for surface
        // triangle mesh (not volumetric tetrahedron). 
        void AddImpulse(const REAL &timestamp, const int &appliedVertex, const Vector3d &impulse); 
        void GetImpulse(const int &index, REAL &timestamp, int &vertex, Vector3d &impulse); 
        void GetImpulse(const REAL &timeStart, const REAL &timeStop, Vector3d &impulse); 
        void GetImpulseRange(REAL &firstImpulseTime, REAL &lastImpulseTime);
};

#endif 
