//////////////////////////////////////////////////////////////////////
// GTS_TriMesh.cpp: Implementation of the GTS_TriMesh class
//
//////////////////////////////////////////////////////////////////////

#include "GTS_TriMesh.h"

#include <utils/trace.h>

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
GTS_TriMesh::GTS_TriMesh( const TriangleMesh<REAL> &mesh )
    : _mesh( mesh ),
      _surface( NULL )
{
    initSurface();
#if 0
    printf( "GTS surface has %d vertices\n",
            (int)gts_surface_vertex_number( _surface ) );
    printf( "GTS surface has %d faces\n",
            (int)gts_surface_face_number( _surface ) );

    if ( gts_surface_is_manifold( _surface ) )
    {
        cout << "Surface is manifold" << endl;
    }
    else
    {
        cout << "Surface is NOT manifold" << endl;
    }

    if ( gts_surface_is_closed( _surface ) )
    {
        cout << "Surface is closed" << endl;
    }
    else
    {
        cout << "Surface is NOT closed" << endl;
    }
#endif
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
GTS_TriMesh::~GTS_TriMesh()
{
    clear();
}

//////////////////////////////////////////////////////////////////////
// Gaussian curvature at a vertex
//////////////////////////////////////////////////////////////////////
REAL GTS_TriMesh::gaussianCurvature( int vertex_idx )
{
    gboolean                   success;
    gdouble                    Kg;

    success = gts_vertex_gaussian_curvature( _vertices.at( vertex_idx ),
                                             _surface, &Kg );

    if ( success == FALSE )
    {
        cerr << "Error computing Gaussian curvature at vertex "
            << vertex_idx << " ( " << Kg << " )" << endl;
        Kg = 0.0;
    }

    return Kg;
}

//////////////////////////////////////////////////////////////////////
// Mesh curvature at a vertex
//////////////////////////////////////////////////////////////////////
Vector3d GTS_TriMesh::meanCurvature( int vertex_idx )
{
    gboolean                   success;
    GtsVector                  Kh;
    Vector3d                   Kh_retval;

    success = gts_vertex_mean_curvature_normal( _vertices.at( vertex_idx ),
                                                _surface, Kh );

    if ( success == FALSE )
    {
        cerr << "Error computing mean curvature at vertex "
             << vertex_idx << " ( "
             << Kh[ 0 ] << ", " << Kh[ 1 ] << ", " << Kh[ 2 ] << " )" << endl;

        Kh[ 0 ] = 0.0;
        Kh[ 1 ] = 0.0;
        Kh[ 2 ] = 0.0;
    }

    Kh_retval[ 0 ] = Kh[ 0 ];
    Kh_retval[ 1 ] = Kh[ 1 ];
    Kh_retval[ 2 ] = Kh[ 2 ];

    return Kh_retval;
}

//////////////////////////////////////////////////////////////////////
// Samples mean curvature on a triangle: assumes that mean
// curvatures have been precomputed
//////////////////////////////////////////////////////////////////////
REAL GTS_TriMesh::meanCurvature( int triangle_idx,
                                 const Vector3d &barycentricPosition ) const
{
    REAL                       curvature = 0.0;

    TRACE_ASSERT( _meanCurvatures.size() == _vertices.size() );

    const Tuple3ui            &triangle = _mesh.triangles()[ triangle_idx ];

    for ( int i = 0; i < 3; i++ )
    {
        curvature
            //+= barycentricPosition[ i ] * _meanCurvatures[ triangle->getIndex( 0 ) ];
            += barycentricPosition[ i ] * _meanCurvatures[ triangle[ i ] ];
    }

    return curvature;
}

//////////////////////////////////////////////////////////////////////
// Gets all curvature information at once at a vertex
//////////////////////////////////////////////////////////////////////
void GTS_TriMesh::allCurvatures( int vertex_idx,
                                 REAL &Kg, Vector3d &Kh,
                                 REAL &principalCurvature1,
                                 REAL &principalCurvature2 )
{
    gdouble                    K1, K2;

    Kg = gaussianCurvature( vertex_idx );
    Kh = meanCurvature( vertex_idx );

    //gts_vertex_principal_curvatures( Kg, 0.5 * norm( Kh ), &K1, &K2 );
    gts_vertex_principal_curvatures( Kg, 0.5 * Kh.length(), &K1, &K2 );

    principalCurvature1 = K1;
    principalCurvature2 = K2;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void GTS_TriMesh::precomputeMeanCurvatures( bool smooth )
{
    _meanCurvatures.resize( _vertices.size() );

    for ( size_t vertex_idx = 0; vertex_idx < _vertices.size(); vertex_idx++ ) {
        //_meanCurvatures[ vertex_idx ] = norm( meanCurvature( vertex_idx ) ) / 2.0;
        _meanCurvatures[ vertex_idx ] = meanCurvature( vertex_idx ).length() / 2.0;
    }

    if ( smooth ) {
        smoothCurvatures();
    }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void GTS_TriMesh::initSurface()
{
    const vector<Point3d>       &vertices = _mesh.vertices();
    const vector<Tuple3ui>      &triangles = _mesh.triangles();
    //const Vector3Array        &vertices = _mesh.vertices();
    //const vector<Triangle *>  &triangles = _mesh.triangles();

    typedef UnorderedMap<pair<int, int>, GtsEdge *>::type  EdgeMap;

    EdgeMap                    edges;

    GtsEdge                   *e1;
    GtsEdge                   *e2;
    GtsEdge                   *e3;

    _surface = gts_surface_new( gts_surface_class(),
                                gts_face_class(),
                                gts_edge_class(),
                                gts_vertex_class() );

    for ( size_t vertex_idx = 0; vertex_idx < vertices.size(); vertex_idx++ )
    {
        const Point3d           &x = vertices[ vertex_idx ];

        _vertices.push_back( gts_vertex_new( gts_vertex_class(),
                                             x[ 0 ], x[ 1 ], x[ 2 ] ) );
    }

    // Generate the faces and edges
    for ( size_t tri_idx = 0; tri_idx < triangles.size(); tri_idx++ )
    {
        const Tuple3ui          &triangle = triangles[ tri_idx ];

        pair<int, int>           edgePair1( min( triangle[ 0 ], triangle[ 1 ] ),
                                            max( triangle[ 0 ], triangle[ 1 ] ) );
        pair<int, int>           edgePair2( min( triangle[ 1 ], triangle[ 2 ] ),
                                            max( triangle[ 1 ], triangle[ 2 ] ) );
        pair<int, int>           edgePair3( min( triangle[ 2 ], triangle[ 0 ] ),
                                            max( triangle[ 2 ], triangle[ 0 ] ) );

        // Get the three face edges.  If they don't exist already, build
        // new ones
        if ( edges.find( edgePair1 ) != edges.end() ) {
            e1 = edges[ edgePair1 ];
        } else {
            e1 = gts_edge_new( gts_edge_class(),
                               _vertices[ triangle[ 0 ] ],
                               _vertices[ triangle[ 1 ] ] );
            edges[ edgePair1 ] = e1;
        }

        if ( edges.find( edgePair2 ) != edges.end() ) {
            e2 = edges[ edgePair2 ];
        } else {
            e2 = gts_edge_new( gts_edge_class(),
                               _vertices[ triangle[ 1 ] ],
                               _vertices[ triangle[ 2 ] ] );
            edges[ edgePair2 ] = e2;
        }

        if ( edges.find( edgePair3 ) != edges.end() ) {
            e3 = edges[ edgePair3 ];
        } else {
            e3 = gts_edge_new( gts_edge_class(),
                               _vertices[ triangle[ 2 ] ],
                               _vertices[ triangle[ 0 ] ] );
            edges[ edgePair3 ] = e3;
        }


        _faces.push_back( gts_face_new( gts_face_class(), e1, e2, e3 ) );

        gts_surface_add_face( _surface, _faces.back() );
    }

#if 0
    // FIXME: debugging
    for ( int tri_idx = 0; tri_idx < _faces.size(); tri_idx++ )
    {
        gboolean                 face_check;
        guint                    face_count;

        face_check = gts_face_has_parent_surface( _faces[ tri_idx ], _surface );
        face_count = gts_face_neighbor_number( _faces[ tri_idx ], _surface );

        cout << "Face " << tri_idx << " has parent surface: " << face_check << endl;
        cout << "Face " << tri_idx << " has " << face_count
            << " neighbours" << endl;

        face_check = gts_face_is_compatible( _faces[ tri_idx ], _surface );

        cout << "Face " << tri_idx << " is compatible: " << face_check << endl;
    }
#endif
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void GTS_TriMesh::clear()
{
    if ( _surface ) {
        gts_object_destroy( GTS_OBJECT( _surface ) );
    }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void GTS_TriMesh::smoothCurvatures()
{
    FloatArray                   triangleCurvatures( _mesh.triangles().size(), 0.0 );

    FloatArray                   triangleCounts( _mesh.vertices().size(), 0.0 );

    for ( size_t tri_idx = 0; tri_idx < _mesh.triangles().size(); tri_idx++ ) {
        const Tuple3ui          &tri = _mesh.triangles()[ tri_idx ];

        for ( int i = 0; i < 3; i++ ) {
            triangleCurvatures[ tri_idx ] += _meanCurvatures[ tri[ i ] ];
        }
        triangleCurvatures[ tri_idx ] /= 3.0;
    }

    _meanCurvatures.clear();
    _meanCurvatures.resize( _mesh.vertices().size(), 0.0 );

    for ( size_t tri_idx = 0; tri_idx < _mesh.triangles().size(); tri_idx++ ) {
        const Tuple3ui          &tri = _mesh.triangles()[ tri_idx ];

        for ( int i = 0; i < 3; i++ ) {
            _meanCurvatures[ tri[ i ] ] += triangleCurvatures[ tri_idx ];
            triangleCounts[ tri[ i ] ] += 1.0;
        }
    }

    for ( size_t vert_idx = 0; vert_idx < _meanCurvatures.size(); vert_idx++ ) {
        _meanCurvatures[ vert_idx ] /= triangleCounts[ vert_idx ];

        printf( "Vertex %lu smoothed curvature: %f\n",
                vert_idx, _meanCurvatures[ vert_idx ] );
        printf( "Vertex %lu smoothed radius: %f\n\n\n",
                vert_idx, 1.0 / _meanCurvatures[ vert_idx ] );
    }
}
