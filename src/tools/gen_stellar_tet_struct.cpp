/*
 * Read the given tetrahedron mesh, and use that to extract the matrices
 * such that the stiffness matrix and mass matrix
 */
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include "utils/print_msg.h"
#include "geometry/FixVtxTetMesh.hpp"
#include "geometry/TriangleMesh.hpp"
#include "io/StellarIO.hpp"
#include "io/MatrixIO.hpp"
#include "io/TglMeshWriter.hpp"
#include "deformable/stvk.h"
#include "deformable/fem.hpp"

using namespace std;

typedef FixVtxTetMesh<double>               TMesh;

static TMesh                mesh;
static TriangleMesh<double> tglmesh;
static double youngModulus;
static double poissonRatio;
static int    outType = 0;
static string tetfile = "";

static vector<Vector3d> tglNml;
static vector<Vector3d> vtxNml;

static void usage(char* cmd)
{
    printf("Usage: %s -f <tet_file> -v -y <Young's modulus value> -p <Poisson's ratio> [-m -k]\n", cmd);
    printf("    -f <tet_file>       Specify the file name tetrahedron mesh, no extension is needed.\n");
    printf("                        The extension of the tetrahedron mesh should be .tet\n");
    printf("    -v                  Specify if the tetrahedron file contain a mesh with fixed\n");
    printf("                        vertices. By default, it assumes the input mesh is\n");
    printf("                        constraint-free\n");
    printf("    -y <value>          Young's modulus value\n");
    printf("    -p <value>          Poisson's ratio\n");
    printf("    -m                  Output mass matrix, named as <tet_file>.mass.csc\n");
    printf("    -k                  Output stiffness matrix, named as <tet_file>.stiff.csc\n");
    printf("    -g                  Output geometry file, named as <tet_file>.geo.txt\n");
    printf("    -s                  Output the surface triangle mesh\n");
}

static void parse_cmd(int argc, char* argv[])
{
    int opt;
    while ( (opt = getopt(argc, argv, "hf:y:p:vmkgs")) != -1 )
    {
        switch (opt)
        {
            case 'h':
                usage(argv[0]);
                exit(0);
            case 'f':
                tetfile = optarg;
                break;
            case 'y':
                youngModulus = atof(optarg);
                break;
            case 'p':
                poissonRatio = atof(optarg);
                break;
            case 'm':
                outType |= 1;
                break;
            case 'k':
                outType |= 2;
                break;
            case 'g':
                outType |= 4;
                break;
            case 's':
                outType |= 8;
                break;
        }
    }
}

static void compute_tgl_pseudo_normals()
{
    const vector<Point3d>&  vtx = tglmesh.vertices();
    const vector<Tuple3ui>& tgl = tglmesh.surface_indices();

    tglNml.resize(tgl.size());
    for(size_t i = 0;i < tgl.size();++ i)
    {
        tglNml[i] = Triangle<double>::normal(
                vtx[tgl[i][0]], vtx[tgl[i][1]], vtx[tgl[i][2]]);
        if ( tglNml[i].lengthSqr() < 1E-24 )
        {
            fprintf(stderr, "ERROR: triangle has zero area: %.30lf\n",
                    tglNml[i].lengthSqr());
            exit(1);
        }
        tglNml[i].normalize();
    }
}

static void compute_vtx_pseudo_normals()
{
    const vector<Point3d>&  vtx = tglmesh.vertices();
    const vector<Tuple3ui>& tgl = tglmesh.surface_indices();

    vtxNml.resize(vtx.size());
    memset(&vtxNml[0], 0, sizeof(Point3d)*vtxNml.size());

    for(size_t i = 0;i < tgl.size();++ i)
    {
        const Vector3<double>& nml = tglNml[i];

        vtxNml[tgl[i][0]] += nml * Triangle<double>::angle(
              vtx[tgl[i][2]], vtx[tgl[i][0]], vtx[tgl[i][1]]);
        vtxNml[tgl[i][1]] += nml * Triangle<double>::angle(
              vtx[tgl[i][0]], vtx[tgl[i][1]], vtx[tgl[i][2]]);
        vtxNml[tgl[i][2]] += nml * Triangle<double>::angle(
              vtx[tgl[i][1]], vtx[tgl[i][2]], vtx[tgl[i][0]]);
    }
}

static void load_mesh()
{
    printf("Load mesh [%s] ... ", tetfile.c_str());
    StellarTetMeshLoader::load_mesh(tetfile.c_str(), &mesh);
    printf(" [OK]\n");

    printf("Initialize mesh ... ");
    mesh.init();
    mesh.update_surface();

    if ( outType & 4 || outType & 8 )
        mesh.extract_surface(tglmesh);
    if ( outType & 4 )
    {
        tglmesh.update_vertex_areas();
        compute_tgl_pseudo_normals();
        compute_vtx_pseudo_normals();
    }
    printf(" [OK]\n");
}

static void output_mass_matrix()
{
    char fname[128];
    sprintf(fname, "%s_mass.csc", tetfile.c_str());

    PardisoMatrix<double> M;
    DeformableFEM::mass_mat(&mesh, M);
    PardisoMatrixIO::write_csc_format(M, fname);
}

static void output_stiff_matrix()
{
    char fname[128];
    sprintf(fname, "%s_stiff.csc", tetfile.c_str());

    PardisoMatrix<REAL> K;
    StVKMaterial material(
            (youngModulus*poissonRatio)/((1.+poissonRatio)*(1.-2.*poissonRatio)),
            youngModulus / (2.*(1.+poissonRatio)));
    DeformableFEM::stiffness_mat(&mesh, &material, K);
    PardisoMatrixIO::write_csc_format(K, fname);
}

static void output_surface_mesh()
{
    char fname[128];
    sprintf(fname, "%s.obj", tetfile.c_str());
    MeshObjWriter::write(tglmesh, fname);
}

static void output_geometry()
{
    char fname[128];
    sprintf(fname, "%s_geo.txt", tetfile.c_str());
    ofstream fout(fname);
    if ( fout.fail() )
    {
        cerr << "ERROR: Cannot open file: " << fname << endl;
        return;
    }

    map<int, int> idmap;
    mesh.surface_id_map(idmap);
    const map<int, int>::iterator end = idmap.end();
    const valarray<double>& area = tglmesh.vertex_areas();

    fout << idmap.size() << endl;
    fout << setprecision(20);
    for(map<int,int>::iterator it = idmap.begin();it != end;++ it)
    {
        Vector3d n = vtxNml[it->second];
        n.normalize();
        fout << it->first << ' ' << it->second << ' '
             << n.x << ' ' << n.y << ' ' << n.z << ' '
             << area[it->second] << endl;
    }
    fout.close();
}

int main(int argc, char* argv[])
{
    parse_cmd(argc, argv);
    if ( youngModulus <= 1E-8 || poissonRatio <= 1E-8 )
    {
        PRINT_ERROR("Unreasonable Young's modulus or Poisson ratio is given");
        exit(1);
    }
    if ( tetfile.empty() ) 
    {
        PRINT_ERROR("No tetrahedron file is given");
        exit(1);
    }
    if ( !(outType & 15) )
    {
        PRINT_WARNING("No output type specified");
        exit(0);
    }

    load_mesh();
    if ( outType & 1 ) output_mass_matrix();
    if ( outType & 2 ) output_stiff_matrix();
    if ( outType & 4 ) output_geometry();
    if ( outType & 8 ) output_surface_mesh();

    return 0;
}

