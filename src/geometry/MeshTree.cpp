//////////////////////////////////////////////////////////////////////
// MeshTree.cpp: Implementation of the MeshTree class
//
//////////////////////////////////////////////////////////////////////

#include "MeshTree.h"

#include <utils/IO.h>
#include <utils/trace.h>

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
MeshTree::MeshTree()
#ifdef USE_CGAL
    : _sceneTree( NULL )
#endif
{
#ifndef USE_CGAL
    cerr << "Mesh Tree not implemented; aborting" << endl;
    abort();
#endif
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
MeshTree::~MeshTree()
{
#ifdef USE_CGAL
    delete _sceneTree;
#endif
}

//////////////////////////////////////////////////////////////////////
// Inserts a mesh in to the tree
//////////////////////////////////////////////////////////////////////
void MeshTree::addMesh( const TriangleMesh<REAL> *mesh )
{
#ifdef USE_CGAL
    int                        full_triangle_idx;
    int                        mesh_idx = _meshes.size();

    const vector<Point3d>     &vertices = mesh->vertices();
    const vector<Tuple3ui>    &triangles = mesh->triangles();

    TRACE_ASSERT( mesh->has_normals(), "No normals defined for this mesh" );

    _meshes.push_back( mesh );

    TRACE_ASSERT( !_sceneTree, "Scene tree is already built" );

    for ( size_t tri_idx = 0; tri_idx < triangles.size(); tri_idx++ )
    {
        const Tuple3ui        &triangle = triangles[ tri_idx ];

        full_triangle_idx = _allTriangles.size();

        _allTriangles.push_back(
                NamedTriangle( Point3d_to_Point( vertices[ triangle[ 0 ] ] ),
                               Point3d_to_Point( vertices[ triangle[ 1 ] ] ),
                               Point3d_to_Point( vertices[ triangle[ 2 ] ] ),
                               full_triangle_idx ) );

        // Push an ID in to the triangle map, so that we can do a reverse
        // lookup later on
        _triangleMap.push_back( MeshReference( mesh_idx, tri_idx ) );
    }
#endif
}

//////////////////////////////////////////////////////////////////////
// Assembles the tree data structure
//////////////////////////////////////////////////////////////////////
void MeshTree::buildTree()
{
#ifdef USE_CGAL
    TRACE_ASSERT( !_sceneTree, "Scene tree is already built" );

    _sceneTree = new TriangleTree( _allTriangles.begin(), _allTriangles.end() );
#endif
}

//////////////////////////////////////////////////////////////////////
// Mainly for testing
//////////////////////////////////////////////////////////////////////
MeshTree::RayHit MeshTree::intersectRay( const Point3d &source,
                                         const Vector3d &direction )
{
#ifdef USE_CGAL
    Ray                        r = BuildRay( source, direction );

    return intersectRay( r );
#endif
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
MeshTree::RayHit MeshTree::closestPoint( const Point3d &x )
{
#ifdef USE_CGAL
    Point                      p( x[ 0 ], x[ 1 ], x[ 2 ] );

    return closestPoint( p );
#endif
}

#ifdef USE_CGAL
//////////////////////////////////////////////////////////////////////
// Intersects the given ray with the scene tree.  Returns a RayHit
// object
//////////////////////////////////////////////////////////////////////
MeshTree::RayHit MeshTree::intersectRay( MeshTree::Ray &r )
{
    typedef std::list<TriangleTree::Object_and_primitive_id>::iterator ListIter;

    RayHit                     hit;

    Point                      source = r.source();
    Point                      hitPoint;
    Vector                     raySegment;
    ListIter                   iter;
    CGAL::Object               object;
    TriangleList::iterator     triangleID;

    REAL                       rayDistance;
    REAL                       closestPoint = FLT_MAX;

    hit._hitReference = MeshReference( -1, -1 );

    std::list<TriangleTree::Object_and_primitive_id> intersections;

    _sceneTree->all_intersections( r, std::back_inserter( intersections ) );

    if ( intersections.size() == 0 )
    {
        return hit;
    }

    iter = intersections.begin();

    for ( ; iter != intersections.end(); iter++ )
    {
        // Get the distance of the first point
        object = iter->first;
        triangleID = iter->second;

        // Get the position of the ray intersection
        if ( CGAL::assign( hitPoint, object ) )
        {
            raySegment = hitPoint - source;
            rayDistance = raySegment.squared_length();

            if ( rayDistance < closestPoint )
            {
                closestPoint = rayDistance;

                hit._hitReference = _triangleMap[ triangleID->index() ];
                hit._position = Point_to_Point3d( hitPoint );
            }
        }
    }

#if 0
    if ( hit._hitReference.first >= 0 )
    {
        hit._normal = _meshes[ hit._hitReference.first ]->triangles()
            [ hit._hitReference.second ]->getNormal();
    }

    // Determine barycentric coordinates for the hit position
    const Triangle            *tri;

    tri = _meshes[ hit._hitReference.first ]->triangles()
        [ hit._hitReference.second ];
#endif

    if ( hit._hitReference.first < 0 ) {
        cout << "Warning: no ray intersection found" << endl;
        return hit;
    }

    const Tuple3ui           &tri = _meshes[ hit._hitReference.first ]->triangles()
                                                                [ hit._hitReference.second ];
    const vector<Point3d>    &vertices = _meshes[ hit._hitReference.first ]->vertices();
    const vector<Vector3d>   &normals = _meshes[ hit._hitReference.first ]->normals();

#if 0
    MATRIX3                    barySystem( tri->getX( 0 ), tri->getX( 1 ), tri->getX( 2 ) );

    barySystem = barySystem.transpose();
    barySystem = barySystem.inverse();
#endif

    // FIXME: this is ugly and potentially ill-conditioned
    Matrix3d                   barySystem( vertices[ tri[ 0 ] ],
                                           vertices[ tri[ 1 ] ],
                                           vertices[ tri[ 2 ] ] );

    barySystem = barySystem.inverse();

    hit._barycentricPosition = barySystem * hit._position;

    hit._normal = 0;
    for ( size_t corner_idx = 0; corner_idx < 3; corner_idx++ ) {
        hit._normal += hit._barycentricPosition[ corner_idx ] * normals[ tri[ corner_idx ] ];
    }

    return hit;
}
#endif

#ifdef USE_CGAL
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
MeshTree::RayHit MeshTree::closestPoint( Point &p )
{
    TriangleTree::Point_and_primitive_id   hitID;
    TriangleList::iterator                 triangleID;
    RayHit                                 hit;

    hitID = _sceneTree->closest_point_and_primitive( p );

    triangleID = hitID.second;

    hit._hitReference = _triangleMap[ triangleID->index() ];
    hit._position = Point_to_Point3d( hitID.first );

#if 0
    if ( hit._hitReference.first >= 0 )
    {
        hit._normal = _meshes[ hit._hitReference.first ]->triangles()
            [ hit._hitReference.second ]->getNormal();
    }

    // Determine barycentric coordinates for the hit position
    const Triangle            *tri;

    tri = _meshes[ hit._hitReference.first ]->triangles()
        [ hit._hitReference.second ];

    MATRIX3                    barySystem( tri->getX( 0 ),
            tri->getX( 1 ),
            tri->getX( 2 ) );

    barySystem = barySystem.transpose();
    barySystem = barySystem.inverse();
#endif

    // This should always work
    TRACE_ASSERT( hit._hitReference.first >= 0, "No closest point found!" );

    // Determine barycentric coordinates for the hit position
    const Tuple3ui            &tri = _meshes[ hit._hitReference.first ]->triangles()
                                                                [ hit._hitReference.second ];
    const vector<Point3d>     &vertices = _meshes[ hit._hitReference.first ]->vertices();
    const vector<Vector3d>    &normals = _meshes[ hit._hitReference.first ]->normals();

    // FIXME: this is ugly and potentially ill-conditioned
    Matrix3d                   barySystem( vertices[ tri[ 0 ] ],
                                           vertices[ tri[ 1 ] ],
                                           vertices[ tri[ 2 ] ] );

    barySystem = barySystem.inverse();

    hit._barycentricPosition = barySystem * hit._position;

    // TODO: hit normal
    hit._normal = 0;
    for ( size_t corner_idx = 0; corner_idx < 3; corner_idx++ ) {
        hit._normal += hit._barycentricPosition[ corner_idx ] * normals[ tri[ corner_idx ] ];
    }

    return hit;
}
#endif
