//////////////////////////////////////////////////////////////////////
// MeshTree.h: Interface for the MeshTree class
//
//////////////////////////////////////////////////////////////////////

#ifndef MESH_TREE_H
#define MESH_TREE_H

#include <TYPES.h>

#include <geometry/TriangleMesh.hpp>

#include <linearalgebra/Matrix3.hpp>
#include <linearalgebra/Vector3.hpp>

#ifdef USE_CGAL
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/Simple_cartesian.h>
#endif

//////////////////////////////////////////////////////////////////////
// MeshTree class
//
// Wrapper around CGAL's AABB tree class.  Used for performing
// ray tracing queries on a mesh (or set of meshes)
//////////////////////////////////////////////////////////////////////
class MeshTree {
    public:
#ifdef USE_CGAL
        typedef CGAL::Simple_cartesian<REAL>   Space;

        // Definitions for CGAL geometric primitives
        typedef Space::Ray_3                   Ray;
        typedef Space::Line_3                  Line;
        typedef Space::Point_3                 Point;
        typedef Space::Vector_3                Vector;
        typedef Space::Triangle_3              SceneTriangle;
#endif

        // Mesh-triangle index pair
        typedef std::pair<int, int>            MeshReference;

        struct RayHit {
            MeshReference                        _hitReference;

            // Hit geometry
            Point3d                              _position;
            Point3d                              _barycentricPosition;
            Vector3d                             _normal;

        };

        MeshTree();

        // Destructor
        virtual ~MeshTree();

        // Inserts a mesh in to the tree
        void addMesh( const TriangleMesh<REAL> *mesh );

        // Assembles the tree data structure
        void buildTree();

        // Ray intersection function
        RayHit intersectRay( const Point3d &source, const Vector3d &direction ); 

        RayHit closestPoint( const Point3d &x );

    protected:

    private:
#ifdef USE_CGAL
        // Intersects the given ray with the scene tree.  Returns a RayHit
        // object
        RayHit intersectRay( Ray &r );

        RayHit closestPoint( Point &p );

        inline static Point Vector3d_to_Point( const Vector3d &x )
        {
            return Point( x[ 0 ], x[ 1 ], x[ 2 ] );
        }

        inline static Vector3d Point_to_Vector3d( const Point &p )
        {
            return Vector3d( p.x(), p.y(), p.z() );
        }

        inline static Vector Vector3d_to_Vector( const Vector3d &x )
        {
            return Vector( x[ 0 ], x[ 1 ], x[ 2 ] );
        }

        inline static Point Point3d_to_Point( const Point3d &x )
        {
            return Point( x[ 0 ], x[ 1 ], x[ 2 ] );
        }

        inline static Point3d Point_to_Point3d( const Point &p )
        {
            return Point3d( p.x(), p.y(), p.z() );
        }

        inline static Vector Point3d_to_Vector( const Point3d &x )
        {
            return Vector( x[ 0 ], x[ 1 ], x[ 2 ] );
        }

        inline static Ray BuildRay( const Point3d &source, const Vector3d &direction )
        {
            return Ray( Point3d_to_Point( source ), Vector3d_to_Vector( direction ) );
        }
#endif

    private:
#ifdef USE_CGAL
        class NamedTriangle : public SceneTriangle {
            public:
                NamedTriangle( Point p, Point q, Point r, int index )
                    : SceneTriangle( p, q, r ),
                    _index( index )
                {
                }

                int index() const
                {
                    return _index;
                }

            private:
                int                  _index;
        };

        typedef std::vector<NamedTriangle>     TriangleList;
        typedef TriangleList::iterator         TriangleIterator;

        typedef CGAL::AABB_triangle_primitive<Space, TriangleIterator> Primitive;
        typedef CGAL::AABB_traits<Space, Primitive> TriangleTraits;
        typedef CGAL::AABB_tree<TriangleTraits> TriangleTree;

        // All meshes in the scene
        std::vector<const TriangleMesh<REAL> *>           _meshes;

        // Data structures used to store a tree representation of the
        // scene (for ray intersection)
        TriangleList                           _allTriangles;

        // Use this to refer back to the appropriate mesh, given a triangle
        // index in the full list
        std::vector<MeshReference>             _triangleMap;

        // The actual tree
        TriangleTree                          *_sceneTree;
#endif

};

#endif
