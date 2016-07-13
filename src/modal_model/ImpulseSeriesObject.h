#ifndef IMPULSE_SERIES_OBJECT_H 
#define IMPULSE_SERIES_OBJECT_H 
#include <string> 
#include <vector> 
#include <memory>
#include <config.h>
#include <linearalgebra/Vector3.hpp>
#include <geometry/TriangleMesh.hpp>

//##############################################################################
// Class that handles impulse series of an object.
//
// TODO move comments
// the impulses created by rigidsim tool in the repo
// The dumped file that this class deals with usually has the hard-coded name 
// "modalImpulses.txt" from the simulator
//##############################################################################
class ImpulseSeriesObject
{
    public: 
        typedef std::shared_ptr<TriangleMesh<REAL> > TriangleMeshPtr; 

    private: 
        // object that is managing this impulse series
        int                 _objectID; 
        TriangleMeshPtr     _objectMesh; 

        // impulse data
        std::string             _impulseFile; 
        int                     _lengthImpulses; 
        REAL                    _firstImpulseTime; 
        REAL                    _lastImpulseTime; 
        std::vector<Vector3d>   _impulses; 
        std::vector<REAL>       _impulseTimestamps; 
        std::vector<int>        _impulseAppliedVertex; 

    public: 
        ImpulseSeriesObject(const int &objectID, const TriangleMeshPtr &meshPtr, const std::string &impulseFile); 

        inline int Size(){return _lengthImpulses;}
        void Initialize(); 
        // Add impulse for this object. The vertex index should be for surface
        // triangle mesh (not volumetric tetrahedron). 
        void AddImpulse(const REAL &timestamp, const int &appliedVertex, const Vector3d &impulse); 
        void GetImpulse(const int &index, REAL &timestamp, int &vertex, Vector3d &impulse); 
        void GetImpulse(const REAL &time); 
        void GetImpulseRange(REAL &firstImpulseTime, REAL &lastImpulseTime);
};

#endif 
