/////////////////////////////////////////////////////////////////////////
// elasticity_solver.cpp
//
// Solver which builds the linear elastic mass and stiffness matrices
// for tetrahedral finite element mesh, given the Poisson ratio and
// Young's modulus for the object.
//
/////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iostream>
#include <map>

#include <boost/program_options.hpp>

#include <config.h>

#include <deformable/fem.hpp>
#include <deformable/stvk.h>

#include <geometry/FixVtxTetMesh.hpp>
#include <geometry/TriangleMesh.hpp>

#include <linearalgebra/PardisoMatrix.hpp>

#include <io/MatrixIO.hpp>
#include <io/TetMeshReader.hpp>

using namespace std;

/////////////////////////////////////////////////////////////////////////
// Writes a *.geo.txt file, storing a list of tet mesh surface triangle
// IDs, as well as normal and area information for those triangles
/////////////////////////////////////////////////////////////////////////
void WriteGeoFile( const char *tetfile, const FixVtxTetMesh<REAL> &mesh );

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
    if ( argc < 4 ) {
        cerr << "Usage: " << argv[ 0 ]
             << " <tet mesh file>"
                " <Young's modulus>"
                " <Poisson ratio>" << endl;
        return 1;
    }

    REAL                     E = atof( argv[2] );
    REAL                     lambda = 0.0, G = 0.0;
    REAL                     nu = atof( argv[3] );

    FixVtxTetMesh<REAL>      mesh;

    // Load the mesh
    FV_TetMeshLoader_Double::load_mesh( argv[1], mesh );

    // Calculate Lame's first parameter based on Young's modulus and
    // Poisson ratio
    //  Source: http://en.wikipedia.org/wiki/Youngs_modulus
    lambda = E * nu / ( (1.0 + nu) * (1.0 - 2.0 * nu) );
    G = E / ( 2.0 * ( 1 + nu ) );

    // Create St. Venant-Kirchhoff material model
    StVKMaterial             material( lambda, G );

    // Mass (unscaled) and stiffness matrices
    PardisoMatrix<REAL>      mass, stiffness;

    printf( "Generating system for mesh with %d tets\n", (int)mesh.tets().size() );
    printf( "Generating system for mesh with %d tets\n", (int)mesh.num_tets() );

    // Assemble the system matrices
    DeformableFEM::mass_mat( &mesh, mass );
    DeformableFEM::stiffness_mat( &mesh, &material, stiffness );

    // Write the matrices to disk
    //
    // Choose an output file prefix
    const char              *outputPrefix = ( argc >= 5 ) ? argv[4] : argv[1];
    char                     outputFile[1024];

    // Write binary files
    sprintf( outputFile, "%s_mass.mat", outputPrefix );
    PardisoMatrixIO::write_csc_format( mass, outputFile );

    sprintf( outputFile, "%s_stiffness.mat", outputPrefix );
    PardisoMatrixIO::write_csc_format( stiffness, outputFile );

    // Write files in MATLAB format
    sprintf( outputFile, "%s_mass.bcsm", outputPrefix );
    PardisoMatrixIO::write_matlab_binary( mass, outputFile );

    sprintf( outputFile, "%s_stiffness.bcsm", outputPrefix );
    PardisoMatrixIO::write_matlab_binary( stiffness, outputFile );

    // Write out the surface <==> tet mesh mapping file (which we need for some reason)
    WriteGeoFile( argv[1], mesh );

#if 0
    const StVKMesh::TStiffMat &stiffness = materialMesh.update_stiffness_matrix();

    // Write the stiffness matrix to disk
    if ( argc < 5 ) {
        char                 outputFile[ 1024 ];

        sprintf( outputFile, "%s_stiffness.mat", argv[1] );

        PardisoMatrixIO::write_csc_format( stiffness, outputFile );
    } else {
        PardisoMatrixIO::write_csc_format( stiffness, argv[4] );
    }
#endif

    return 0;
}

/////////////////////////////////////////////////////////////////////////
// Writes a *.geo.txt file, storing a list of tet mesh surface triangle
// IDs, as well as normal and area information for those triangles
/////////////////////////////////////////////////////////////////////////
void WriteGeoFile( const char *tetfile, const FixVtxTetMesh<REAL> &mesh )
{
    char                     geoFileName[ 1024 ];
 
    sprintf( geoFileName, "%s.geo.txt", tetfile );

    ofstream                 fout( geoFileName );
    map<int, int>            surfaceIDMap;

    TriangleMesh<REAL>       surfaceMesh;

    if ( fout.fail() ) {
        cerr << "ERROR: Cannot open file: " << geoFileName << endl;
    }

    mesh.extract_surface( &surfaceMesh );
    surfaceMesh.generate_normals();
    surfaceMesh.update_vertex_areas();

    // Get a mapping from tet mesh vertex IDs to triangle mesh vertex IDs only on the surface
    mesh.surface_id_map( surfaceIDMap );

    const valarray<REAL>    &areas = surfaceMesh.vertex_areas();
    const vector<Vector3d>  &normals = surfaceMesh.normals();

    fout << surfaceIDMap.size() << endl;
    fout << setprecision(20);
    for ( map<int, int>::iterator it = surfaceIDMap.begin(); it != surfaceIDMap.end(); it++ ) {
        const Vector3d      &n = normals[it->second];
        REAL                 area = areas[it->second];

        fout << it->first << ' ' << it->second << ' '
             << n.x << ' ' << n.y << ' ' << n.z << ' '
             << area << endl;
    }
    fout.close();
}
