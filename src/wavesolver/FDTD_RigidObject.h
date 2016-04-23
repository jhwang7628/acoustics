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
#include <Eigen/Geometry> 
#include <Eigen/Dense> 

//##############################################################################
// Stores FDTD simulation Object that contains mesh, sdf, transformation matrix
// and can be used to query for scattering, collision etc.
//##############################################################################
class FDTD_RigidObject
{
    public: 
        typedef Eigen::Transform<double,3,Eigen::Affine> Affine3; 
        // bounding box to check 
        struct BoundingBox
        {
            Point3d     minBound; 
            Point3d     maxBound; 
            Vector3d    dimension; 
            Vector3d    centroid; 

            BoundingBox() 
                : minBound(Point3d(F_LOW,F_LOW,F_LOW)), 
                  maxBound(Point3d(F_MAX,F_MAX,F_MAX)), 
                  dimension(Vector3d(F_MAX,F_MAX,F_MAX)), 
                  centroid(Vector3d(0,0,0))
            {
            }
            inline bool Inside(const double &x, const double &y, const double &z)
            {
                return (x<maxBound[0]&&x>minBound[0] &&
                        y<maxBound[1]&&y>minBound[1] &&
                        z<maxBound[2]&&z>minBound[2]); 
            }
            inline bool Inside(const double &x, const double &y, const double &z, const double &scale)
            {
                const Vector3d scaledDimension = dimension*scale; 
                return (x<(centroid[0]+scaledDimension[0])&&x>(centroid[0]-scaledDimension[0]) &&
                        y<(centroid[1]+scaledDimension[1])&&y>(centroid[1]-scaledDimension[1]) &&
                        z<(centroid[2]+scaledDimension[2])&&z>(centroid[2]-scaledDimension[2])); 
            }
            inline void Update(const Point3<double> &minBound_c, const Point3<double> &maxBound_c)
            {
                minBound = minBound_c; 
                maxBound = maxBound_c; 
                dimension = maxBound - minBound; 
                centroid = (maxBound + minBound)/2.0;
            }
        };
        

    private: 
        int                                 _meshID; 
        REAL                                _meshScale; 
        std::string                         _meshName;
        std::string                         _meshFileName;
        std::shared_ptr<TriangleMesh<REAL>> _mesh; 
        std::vector<VibrationalSourcePtr>   _vibrationalSources; 

        // bounding box in the work coordinate system
        BoundingBox                         _bboxWorld; 

        int                                 _signedDistanceFieldResolution;
        std::string                         _signedDistanceFieldFilePrefix; 
        std::shared_ptr<ClosestPointField>  _signedDistanceField; 

        // premultiply this transform to point takes it from object space to
        // work space
        Affine3                             _modelingTransform; 
        Affine3                             _modelingTransformInverse; 

        bool                                _parsed; 


    public: 
        FDTD_RigidObject()
            : _meshScale(1.0), 
              _meshName("NOT_IDENTIFIED"),
              _meshFileName("NOT_IDENTIFIED"), 
              _modelingTransform(Affine3::Identity()), 
              _modelingTransformInverse(Affine3::Identity()),
              _parsed(false) 
        {
        }

        FDTD_RigidObject(const std::string &fileName, const int &resolution, const std::string &sdfFilePrefix, const std::string &meshName="NOT_IDENTIFIED", const int &scale=1.0)
            : _meshScale(scale), 
              _meshName(meshName),
              _meshFileName(fileName), 
              _signedDistanceFieldResolution(resolution), 
              _signedDistanceFieldFilePrefix(sdfFilePrefix), 
              _modelingTransform(Affine3::Identity()), 
              _modelingTransformInverse(Affine3::Identity())
        {
            _parsed = true; 
        }

        inline bool Exist(){return (_mesh!=0);}
        inline std::string &GetMeshName(){return _meshName;} 
        inline void AddVibrationalSource(VibrationalSourcePtr &sourcePtr){_vibrationalSources.push_back(std::move(sourcePtr));}
        void Initialize(); 
        void UpdateBoundingBox(); 
        // in-place query for object sdf distance from world x,y,z
        REAL DistanceToMesh(const double &x, const double &y, const double &z); 
        // in-place query for object sdf normal from world x,y,z
        // return success or not (could be that query point is outside of bbox,
        // then normal is not defined
        bool NormalToMesh(const double &x, const double &y, const double &z, Vector3d &queriedNormal); 
        REAL EvaluateBoundaryCondition(const Vector3d &boundaryPoint, const Vector3d &boundaryNormal, const REAL &time); 

        //// debug methods //// 
        void PrintBoundingBox(); 
        void TestQueryDistance(); 
        void TestObjectBoundaryCondition();
};
//##############################################################################


#endif 
