//////////////////////////////////////////////////////////////////////
// ClosestPointMesh.cpp: Implementation
//
//////////////////////////////////////////////////////////////////////

#include "ClosestPointMesh.h"

#include <utils/IO.h>
#include <utils/trace.h>

#include <float.h>

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
#if 0
ClosestPointMesh::ClosestPointMesh( const TriMesh &mesh,
                                    int sdfResolution,
                                    const string &sdfFilePrefix )
    : _mesh( mesh ),
      _sdf( NULL ),
      _rayTree( NULL )
{
    _sdf = mesh.buildDistanceField( sdfResolution, sdfFilePrefix.c_str() );

    _rayTree = new MeshTree();
    _rayTree->addMesh( &mesh );
    _rayTree->buildTree();
}
#endif
ClosestPointMesh::ClosestPointMesh( const TriangleMesh<REAL> &mesh )
    : _mesh( mesh )
{
    _rayTree = new MeshTree();
    _rayTree->addMesh( &mesh );
    _rayTree->buildTree();
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
ClosestPointMesh::~ClosestPointMesh()
{
#if 0
    delete _sdf;
#endif
    delete _rayTree;
}

//////////////////////////////////////////////////////////////////////
// Returns the closest point on the mesh to the given position
//////////////////////////////////////////////////////////////////////
ClosestPointMesh::MeshPoint ClosestPointMesh::closestPoint( const Point3d &x ) const
{
#if 0
    Vector3d                   gradient;
#endif
    MeshTree::RayHit           hit1, hit2;
#if 0
    Real                       distance1 = FLT_MAX;
    Real                       distance2 = FLT_MAX;
#endif

#if 0
    bool                       hitFound = false;
#endif

    MeshPoint                  pointResult;

#if 0
    // Figure out a direction along which to search for the closest point
    gradient = _sdf->gradient( x );
    gradient.normalize();

    hit1 = _rayTree->intersectRay( x, gradient );

    if ( hit1._hitReference.first == 0 )
    {
        // Found a mesh point
        distance1 = norm( hit1._position - x );

        hitFound = true;
    }

    hit2 = _rayTree->intersectRay( x, -1.0 * gradient );

    if ( hit2._hitReference.first == 0 )
    {
        // Found a mesh point
        distance2 = norm( hit2._position - x );

        hitFound = true;
    }

    if ( hitFound )
    {
        // Get a triangle ID from the ray hit
        pointResult._triangleID = ( distance1 < distance2 ) ? 
            hit1._hitReference.second
            : hit2._hitReference.second;

        // Get the position and barycentric coordinates of the hit
        pointResult._position = ( distance1 < distance2 ) ?
            hit1._position : hit2._position;

        pointResult._barycentricPosition = ( distance1 < distance2 ) ?
            hit1._barycentricPosition
            : hit2._barycentricPosition;
    }
    else
    {
        pointResult._triangleID = -1;
    }

    if ( !hitFound )
    {
        cerr << "Warning: no closest point found on mesh" << endl;
    }
    //TRACE_ASSERT( hitFound, "Did not find a point on the mesh" );
#endif

    hit1 = _rayTree->closestPoint( x );

    if ( hit1._hitReference.first != 0
            || hit1._hitReference.second < 0
            || hit1._hitReference.second >= _mesh.triangles().size() )
    {
        cerr << "Warning: no closest point found on mesh" << endl;

        pointResult._triangleID = -1;
    }
    else
    {
        pointResult._triangleID = hit1._hitReference.second;
        pointResult._position = hit1._position;
        pointResult._barycentricPosition = hit1._barycentricPosition;
    }

    return pointResult;
}
