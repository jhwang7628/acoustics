/////////////////////////////////////////////////////////////////////////
// single_tet.cpp
//
// Builds a tet mesh with a single tet and writes it to disk
//
/////////////////////////////////////////////////////////////////////////

#include <geometry/FixVtxTetMesh.hpp>

#include <io/TetMeshWriter.hpp>

#include <iostream>

using namespace std;

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
int main( int argc, char **argv )
{
    FixVtxTetMesh<REAL>      mesh;

    if ( argc < 2 ) {
        cerr << "Usage: " << argv[ 0 ] << " <output file name>" << endl;
        return 1;
    }

    // Add four vertices along the axes
    mesh.add_vertex( Point3d( 0.0, 0.0, 0.0 ) );
    mesh.add_vertex( Point3d( 0.0, 0.0, 1.0 ) );
    mesh.add_vertex( Point3d( 1.0, 0.0, 0.0 ) );
    mesh.add_vertex( Point3d( 0.0, 1.0, 0.0 ) );

    mesh.add_tet( 0, 1, 2, 3 );

    mesh.init();

    FV_TetMeshWriter_Double::write_mesh( argv[ 1 ], mesh );

    return 0;
}
