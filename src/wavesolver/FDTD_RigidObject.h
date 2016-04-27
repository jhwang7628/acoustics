#ifndef FDTD_RIGID_OBJECT_H 
#define FDTD_RIGID_OBJECT_H 

#include <TYPES.h>
#include <config.h>
#include <distancefield/closestPointField.h>
#include <distancefield/FieldBuilder.h>
#include <parser/Parser.h>
#include <geometry/Point3.hpp>
#include <linearalgebra/Vector3.hpp>
#include <wavesolver/VibrationalSource.h> 
#include <wavesolver/Wavesolver_ConstantsAndTypes.h>
#include <wavesolver/FDTD_MovableObject.h> 
#include <Eigen/Geometry> 
#include <Eigen/Dense> 

//##############################################################################
// Stores FDTD simulation Object that contains mesh, sdf, transformation matrix
// and can be used to query for scattering, collision etc.
//##############################################################################
class FDTD_RigidObject : public FDTD_MovableObject
{
    private: 
        int                                 _meshID; 
        REAL                                _meshScale; 
        std::string                         _meshName;
        std::string                         _meshFileName;
        std::shared_ptr<TriangleMesh<REAL>> _mesh; 
        std::vector<VibrationalSourcePtr>   _vibrationalSources; 

        int                                 _signedDistanceFieldResolution;
        std::string                         _signedDistanceFieldFilePrefix; 
        std::shared_ptr<ClosestPointField>  _signedDistanceField; 

        bool                                _parsed; 

    public: 
        FDTD_RigidObject()
            : _meshScale(1.0), 
              _meshName("NOT_IDENTIFIED"),
              _meshFileName("NOT_IDENTIFIED"), 
              _parsed(false) 
        {
        }

        FDTD_RigidObject(const std::string &fileName, const int &resolution, const std::string &sdfFilePrefix, const std::string &meshName="NOT_IDENTIFIED", const int &scale=1.0)
            : _meshScale(scale), 
              _meshName(meshName),
              _meshFileName(fileName), 
              _signedDistanceFieldResolution(resolution), 
              _signedDistanceFieldFilePrefix(sdfFilePrefix)
        {
            _parsed = true; 
        }

        inline bool Exist(){return (_mesh!=0);}
        inline std::string &GetMeshName(){return _meshName;} 
        inline void AddVibrationalSource(VibrationalSourcePtr &sourcePtr){_vibrationalSources.push_back(std::move(sourcePtr));}
        void Initialize(); 
        virtual void UpdateBoundingBox(); 
        // in-place query for object sdf distance from world x,y,z
        REAL DistanceToMesh(const double &x, const double &y, const double &z); 
        // in-place query for object sdf normal from world x,y,z
        // return success or not (could be that query point is outside of bbox,
        // then normal is not defined
        bool NormalToMesh(const double &x, const double &y, const double &z, Vector3d &queriedNormal); 
        REAL EvaluateBoundaryCondition(const Vector3d &boundaryPoint, const Vector3d &boundaryNormal, const REAL &time); 
        bool ReflectAgainstBoundary(const Vector3d &originalPoint, Vector3d &reflectedPoint, Vector3d &boundaryPoint, Vector3d &erectedNormal, REAL &distanceTravelled);

        //// debug methods //// 
        void TestQueryDistance(); 
        void TestObjectBoundaryCondition();
};
//##############################################################################


#endif 
