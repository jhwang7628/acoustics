/////////////////////////////////////////////////////////////////////////
// elasticity_solver.cpp
//
// Solver which builds the linear elastic mass and stiffness matrices
// for tetrahedral finite element mesh, given the Poisson ratio and
// Young's modulus for the object.
//
/////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <boost/program_options.hpp>

#include <config.h>

#include <deformable/StVKMesh.h>

#include <geometry/FixVtxTetMesh.hpp>

#include <io/MatrixIO.hpp>
#include <io/TetMeshReader.hpp>

using namespace std;

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
    if ( argc < 4 ) {
        cerr << "Usage: " << argv[ 0 ]
             << " <tet mesh file>"
                " <Lame's first parameter>"
                " <Poisson ratio>" << endl;
        return 1;
    }

    REAL                     lambda = atof( argv[2] );
    REAL                     nu = atof( argv[3] );

    FixVtxTetMesh<double>    mesh;

    // Load the mesh
    FV_TetMeshLoader_Double::load_mesh( argv[1], mesh );

    // Initialize construct the stiffness matrix
    StVKMesh                 materialMesh( mesh, lambda, nu );

    const StVKMesh::TStiffMat &stiffness = materialMesh.update_stiffness_matrix();

    // Write the stiffness matrix to disk
    if ( argc < 5 ) {
        char                 outputFile[ 1024 ];

        sprintf( outputFile, "%s_stiffness.mat", argv[1] );

        PardisoMatrixIO::write_csc_format( stiffness, outputFile );
    } else {
        PardisoMatrixIO::write_csc_format( stiffness, argv[4] );
    }

    return 0;
}
