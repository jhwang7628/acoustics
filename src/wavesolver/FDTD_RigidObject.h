#ifndef FDTD_RIGID_OBJECT_H 
#define FDTD_RIGID_OBJECT_H 

#include <TYPES.h>
#include <config.h>
#include <distancefield/closestPointField.h>
#include <distancefield/AdaptiveDistanceField.h>
#include <distancefield/FieldBuilder.h>
#include <parser/Parser.h>
#include <geometry/Point3.hpp>
#include <linearalgebra/Vector3.hpp>
#include <wavesolver/VibrationalSource.h> 
#include <wavesolver/Wavesolver_ConstantsAndTypes.h>
#include <wavesolver/FDTD_MovableObject.h> 
#include <geometry/TriangleMeshKDTree.hpp>
#include <geometry/TetMeshIndexToSurfaceMesh.h> 
#include <Eigen/Geometry> 
#include <Eigen/Dense> 
#include <utils/SimpleTimer.h>
#include <utils/IO/IO.h>

//##############################################################################
// Stores FDTD simulation Object that contains mesh, sdf, transformation matrix
// and can be used to query for scattering, collision etc.
//##############################################################################
class FDTD_RigidObject : public FDTD_MovableObject
{
    protected: 
        std::string                         _workingDirectory;
        std::string                         _objectPrefix; 

        int                                 _meshID; 
        REAL                                _meshScale; 
        std::string                         _meshName;
        std::shared_ptr<TriangleMesh<REAL>> _mesh; 
        Vector3d                            _meshObjectCentroid;  // in object space
        REAL                                _volume = -1.0; // mesh volume that can be estimated if built from tet mesh
        std::shared_ptr<TetMeshIndexToSurfaceMesh> _tetMeshIndexToSurfaceMesh; 

        std::vector<VibrationalSourcePtr>   _vibrationalSources; 

        int                                 _signedDistanceFieldResolution;
        std::string                         _signedDistanceFieldFilePrefix; 
#ifndef USE_ADF
        std::shared_ptr<ClosestPointField>  _signedDistanceField; 
#else
        std::shared_ptr<AdaptiveDistanceField>  _signedDistanceField; 
#endif

        bool                                _parsed; 

        std::vector<Vector3d>               _debugArrowStart; 
        std::vector<Vector3d>               _debugArrowNormal; 

    public: 
        FDTD_RigidObject()
            : _workingDirectory("NOT_IDENTIFIED"),
              _objectPrefix("NOT_IDENTIFIED"),
              _meshID(-1), 
              _meshScale(1.0), 
              _meshName("NOT_IDENTIFIED"),
              _parsed(false) 
        {
        }

        FDTD_RigidObject(const std::string &workingDirectory, const int &resolution, const std::string &objectPrefix, const bool &buildFromTetMesh, const std::string &meshName="NOT_IDENTIFIED", const int &scale=1.0)
            : _workingDirectory(workingDirectory), 
              _objectPrefix(objectPrefix),
              _meshID(-1), 
              _meshScale(scale), 
              _meshName(meshName),
              _signedDistanceFieldResolution(resolution)
        {
            _parsed = true; 
            Initialize(buildFromTetMesh); 
        }

        inline bool Exist(){return (_mesh!=0);}
        inline std::shared_ptr<TriangleMesh<REAL>> GetMeshPtr(){return _mesh;} 
        inline const std::shared_ptr<TriangleMesh<REAL>> &GetMeshPtr() const {return _mesh;} 
#ifndef USE_ADF
        inline std::shared_ptr<ClosestPointField> GetSignedDistanceFieldPtr(){return _signedDistanceField;} 
#else
        inline std::shared_ptr<AdaptiveDistanceField> GetSignedDistanceFieldPtr(){return _signedDistanceField;} 
#endif
        inline std::string &GetMeshName(){return _meshName;} 
        inline void AddVibrationalSource(VibrationalSourcePtr &sourcePtr){_vibrationalSources.push_back(std::move(sourcePtr));}
        inline void SetMeshID(const int &ID){_meshID = ID;}
        inline int GetMeshID(){return _meshID;}
        void Initialize(const bool &buildFromTetMesh); 
        virtual void UpdateBoundingBox(); 
        virtual void ResetUnionBox(); 
        // in-place query for object sdf distance from world x,y,z
        REAL DistanceToMesh(const double &x, const double &y, const double &z); 
        REAL DistanceToMesh(const Vector3d &position); 
        // in-place query for object sdf normal from world x,y,z
        // return success or not (could be that query point is outside of bbox,
        // then normal is not defined
        bool NormalToMesh(const double &x, const double &y, const double &z, Vector3d &queriedNormal); 
        bool NormalToMesh(const Vector3d &position, Vector3d &queriedNormal); 
        REAL EvaluateBoundaryAcceleration(const Vector3d &boundaryPoint, const Vector3d &boundaryNormal, const REAL &time); 
        REAL EvaluateBoundaryVelocity(const Vector3d &boundaryPoint, const Vector3d &boundaryNormal, const REAL &time); 
        bool ReflectAgainstBoundary(const Vector3d &originalPoint, Vector3d &reflectedPoint, Vector3d &boundaryPoint, Vector3d &erectedNormal, REAL &distanceTravelled);
        bool FindImageFreshCell(const Vector3d &originalPoint, Vector3d &imagePoint, Vector3d &boundaryPoint, Vector3d &erectedNormal, REAL &distanceTravelled);
        Vector3d MeshCentroid(); 

        //// debug methods //// 
        void TestQueryDistance(); 
        void TestObjectBoundaryCondition();
        void WriteDebugArrow(const std::string &file);
        void ClearDebugArrow();
        REAL EvaluateAccelerationNoiseAnalytical(const Vector3d &listeningPoint, const REAL &time, const REAL &density, const REAL &soundSpeed, const REAL &sphereRadius); 
};
//##############################################################################


#endif 
