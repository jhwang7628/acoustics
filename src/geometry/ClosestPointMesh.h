//////////////////////////////////////////////////////////////////////
// ClosestPointMesh.h: Interface
//
//////////////////////////////////////////////////////////////////////

#ifndef CLOSEST_POINT_MESH_H
#define CLOSEST_POINT_MESH_H

#if 0
#include <distancefield/distanceField.h>
#endif

#include <geometry/MeshTree.h>
#include <geometry/TriangleMesh.hpp>

#include <linearalgebra/Vector3.hpp>

#include <TYPES.h>

//////////////////////////////////////////////////////////////////////
// ClosestPointMesh class
//
// Comments
//////////////////////////////////////////////////////////////////////
class ClosestPointMesh {
    public:
        // Provide a tri mesh, as well as a SDF resolution and file path
#if 0
        ClosestPointMesh( const TriangleMesh<REAL> &mesh,
                          int sdfResolution, const std::string &sdfFilePrefix );
#endif
        ClosestPointMesh( const TriangleMesh<REAL> &mesh );

        // Destructor
        virtual ~ClosestPointMesh();

        struct MeshPoint {
            int                    _triangleID;

            Point3d                _position;
            Vector3d               _barycentricPosition;
        };

        // Returns the closest point on the mesh to the given position
        MeshPoint closestPoint( const Point3d &x ) const;

#if 0
        const DistanceField *sdf() const
        {
            return _sdf;
        }
#endif

    protected:

    private:
        const TriangleMesh<REAL>    &_mesh;

#if 0
        const DistanceField         *_sdf;
#endif

        MeshTree                    *_rayTree;

};

#endif
