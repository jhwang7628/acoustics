#ifndef FDTD_MOVABLE_OBJECT_H 
#define FDTD_MOVABLE_OBJECT_H 

#include <TYPES.h>
#include <config.h> 
#include <linearalgebra/Vector3.hpp>
#include <geometry/Point3.hpp> 
#include <Eigen/Geometry> 
#include <Eigen/Dense> 

//##############################################################################
// Objects that can be animated using rigid body simulator. has bounding box to
// accelerate collision detection etc. 
//##############################################################################
class FDTD_MovableObject
{
    public: 
        typedef Eigen::Transform<REAL,3,Eigen::Affine> Affine3; 
        // bounding box to check 
        struct BoundingBox
        {
            Vector3d     minBound; 
            Vector3d     maxBound; 
            Vector3d    dimension; 
            Vector3d    centroid; 

            BoundingBox() 
                : minBound(Vector3d(F_LOW,F_LOW,F_LOW)), 
                  maxBound(Vector3d(F_MAX,F_MAX,F_MAX)), 
                  dimension(Vector3d(F_MAX,F_MAX,F_MAX)), 
                  centroid(Vector3d(0,0,0))
            {
            }
            BoundingBox(const Vector3d &i_minBound, const Vector3d &i_maxBound, const Vector3d &i_dimension, const Vector3d &i_centroid) 
                : minBound(i_minBound), 
                  maxBound(i_maxBound), 
                  dimension(i_dimension), 
                  centroid(i_centroid)
            {
            }
            BoundingBox(const Vector3d &i_minBound, const Vector3d &i_maxBound) 
                : minBound(i_minBound), 
                  maxBound(i_maxBound), 
                  dimension(Vector3d(maxBound.x-minBound.x, maxBound.y-minBound.z, maxBound.z-minBound.z)), 
                  centroid(Vector3d((maxBound.x+minBound.x)/2.0, (maxBound.y+minBound.y)/2.0, (maxBound.z+minBound.z)/2.0)) 
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
            inline void Update(const Vector3<double> &minBound_c, const Vector3<double> &maxBound_c)
            {
                minBound = minBound_c; 
                maxBound = maxBound_c; 
                dimension = maxBound - minBound; 
                centroid = (maxBound + minBound)/2.0;
            }
            // union two bounding box to form a bigger box in each directions
            inline void Union(const BoundingBox &inBox) 
            {
                for (int d=0; d<3; ++d) 
                {
                    minBound[d] = std::min<double>(inBox.minBound[d], minBound[d]); 
                    maxBound[d] = std::max<double>(inBox.maxBound[d], maxBound[d]); 
                }
                dimension = maxBound - minBound; 
                centroid = (maxBound + minBound)/2.0; 
            }
        };

    protected: 
        // bounding box in the work coordinate system
        BoundingBox                         _bboxWorld; 
        BoundingBox                         _bboxWorldUnion2Steps; 

        // premultiply this transform to point takes it from object space to
        // world space
        Affine3                             _modelingTransform; 
        Affine3                             _modelingTransformInverse; 

    public: 
        FDTD_MovableObject()
            : _modelingTransform(Affine3::Identity()), 
              _modelingTransformInverse(Affine3::Identity())
        {
        }

        inline const BoundingBox &GetUnionBBox(){return _bboxWorldUnion2Steps;}
        inline bool InsideBoundingBox(const double &x, const double &y, const double &z, const double &scale)
        {
            return _bboxWorld.Inside(x, y, z, scale);
        }
        virtual void UpdateBoundingBox()=0; 
        virtual void ApplyTranslation(const double &x, const double &y, const double &z);
        

        //// debug methods //// 
        void PrintBoundingBox(); 
        void PrintTransformation();
};

#endif 
