#ifndef FDTD_RIGID_OBJECT_H 
#define FDTD_RIGID_OBJECT_H 

#include <TYPES.h>
#include <config.h>
#include <distancefield/closestPointField.h>
#include <distancefield/AdaptiveDistanceField.h>
#include <distancefield/FieldBuilder.h>
#include <parser/Parser.h>
#include <geometry/Point3.hpp>
#include <linearalgebra/Matrix3.hpp>
#include <linearalgebra/Vector3.hpp>
#include <wavesolver/VibrationalSource.h> 
#include <wavesolver/Wavesolver_ConstantsAndTypes.h>
#include <wavesolver/PML_WaveSolver_Settings.h>
#include <wavesolver/FDTD_MovableObject.h> 
#include <geometry/TriangleMeshKDTree.hpp>
#include <geometry/TriangleMeshGraph.hpp>
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
    public: 
        struct OptionalAttributes
        {
            bool isFixed=false; 
        };

    protected: 
        std::string                         _workingDirectory;
        std::string                         _objectPrefix; 

        REAL                                _meshScale; 
        std::string                         _meshName;
        std::shared_ptr<TriangleMesh<REAL>> _mesh; 
        std::shared_ptr<TriangleMeshGraph<REAL>> _meshGraph; 
        Vector3d                            _meshObjectCentroid;  // in object space
        std::shared_ptr<TetMeshIndexToSurfaceMesh> _tetMeshIndexToSurfaceMesh; 
        std::vector<VibrationalSourcePtr>   _vibrationalSources; 
        std::shared_ptr<PML_WaveSolver_Settings> _solverSettings; 

        int                                 _signedDistanceFieldResolution;
        std::string                         _signedDistanceFieldFilePrefix; 
#ifndef USE_ADF
        std::shared_ptr<ClosestPointField>  _signedDistanceField; 
#else
        std::shared_ptr<AdaptiveDistanceField>  _signedDistanceField; 
#endif

        bool                                _parsed; 
        OptionalAttributes                  _optionalAttributes; 

        // volumetric field 
        bool                                _hasVolume = false; 
        REAL                                _volume; // mesh volume that can be estimated if built from tet mesh
        Matrix3d                            _volumeInertiaTensor; // I = _volumeInertiaTensor * density
        Point3d                             _volumeCenter; 

        std::vector<Vector3d>               _debugArrowStart; 
        std::vector<Vector3d>               _debugArrowNormal; 

        bool _animated = false; // if an object is being animated by rbd sim, then this should be set to true

    public: 
        FDTD_RigidObject(const ObjectType &type)
            : FDTD_MovableObject(type),
              _workingDirectory("NOT_IDENTIFIED"),
              _objectPrefix("NOT_IDENTIFIED"),
              _meshScale(1.0), 
              _meshName("NOT_IDENTIFIED"),
              _solverSettings(nullptr),
              _parsed(false) 
        {
        }

        FDTD_RigidObject(const ObjectType &type,
                         const std::string &workingDirectory, 
                         const int &resolution, 
                         const std::string &objectPrefix, 
                         const bool &buildFromTetMesh, 
                         const std::shared_ptr<PML_WaveSolver_Settings> &solverSettings, 
                         const std::string &meshName="NOT_IDENTIFIED", 
                         const int &scale=1.0)
            : FDTD_MovableObject(type),
              _workingDirectory(workingDirectory), 
              _objectPrefix(objectPrefix),
              _meshScale(scale), 
              _meshName(meshName),
              _solverSettings(solverSettings),
              _signedDistanceFieldResolution(resolution)
        {
            _parsed = true; 
            Initialize(buildFromTetMesh); 
        }

        inline bool Animated(){return _animated;}
        inline void SetAnimated(const bool &is){_animated = is;}
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
        inline void ClearVibrationalSources(){_vibrationalSources.clear();}
        inline void SetOptionalAttributes(const OptionalAttributes &attr){_optionalAttributes = attr;}
        inline const OptionalAttributes &GetOptionalAttributes(){return _optionalAttributes;}
        void Initialize(const bool &buildFromTetMesh); 
        virtual void UpdateBoundingBox(); 
        virtual void ResetUnionBox(); 
        virtual void ApplyScale(const REAL scale); 
        // in-place query for object sdf distance from world x,y,z
        virtual REAL DistanceToMesh(const double &x, const double &y, const double &z); 
        virtual REAL DistanceToMesh(const Vector3d &position); 
        // in-place query for object sdf normal from world x,y,z
        // return success or not (could be that query point is outside of bbox,
        // then normal is not defined
        virtual bool NormalToMesh(const double &x, const double &y, const double &z, Vector3d &queriedNormal); 
        virtual bool NormalToMesh(const Vector3d &position, Vector3d &queriedNormal); 
        REAL EvaluateBoundaryAcceleration(const Vector3d &boundaryPoint, const Vector3d &boundaryNormal, const REAL &time, const int &hintTriangle=-1); 
        REAL EvaluateBoundaryAcceleration(const int &vertexID, const Vector3d &vertexNormal, const REAL &time); 
        Vector3d EvaluateBoundaryAcceleration(const int &vertexID, const REAL &time); 
        REAL EvaluateBoundaryVelocity(const Vector3d &boundaryPoint, const Vector3d &boundaryNormal, const REAL &time); 
        int ReflectAgainstBoundary(const Vector3d &originalPoint, Vector3d &reflectedPoint, Vector3d &boundaryPoint, Vector3d &erectedNormal, REAL &distanceTravelled, const int &startFromTriangle=-1);
        bool FindImageFreshCell(const Vector3d &originalPoint, Vector3d &imagePoint, Vector3d &boundaryPoint, Vector3d &erectedNormal, REAL &distanceTravelled);
        void SetRigidBodyTransform(const Point3d &newCOM, const Quaternion<REAL> &quaternion); 
        Vector3d MeshCentroid(); 

        //// debug methods //// 
        void TestQueryDistance(); 
        void TestObjectBoundaryCondition();
        void WriteDebugArrow(const std::string &file);
        void ClearDebugArrow();
        REAL EvaluateAccelerationNoiseAnalytical(const Vector3d &listeningPoint, const REAL &time, const REAL &density, const REAL &soundSpeed, const REAL &sphereRadius); 
        int FindLowestVertex(const int &dimension, Vector3d &position); 
};
//##############################################################################

#endif 
