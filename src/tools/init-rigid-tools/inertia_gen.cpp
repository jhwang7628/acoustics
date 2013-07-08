/*
 * compute the volumn and three principle values & correspoding eigenvectors 
 * of the inertia tensor
 */
#include <iostream>
#include <fstream>
#include <iomanip>
#include "geometry/FixVtxTetMesh.hpp"
#include "io/TetMeshReader.hpp"
#include "utils/tuple.hpp"

using namespace std;

typedef FixVtxTetMesh<double>       TTetMesh;

static string tetfileptn;
static int startid, endid;

static double inertia_tensor(const FixVtxTetMesh<double>& mesh, Matrix3<double>& I, Point3d& x0)
{
    const std::vector<double>&    ms  = mesh.masses();
    const std::vector< Point3d >& vtx = mesh.rest_positions();
    double mass = 0;
    x0.zero();

    for(size_t i = 0;i < ms.size();++ i)
    {
        mass += ms[i];
        x0.scaleAdd(ms[i], vtx[i]);
    }
    x0 /= mass; // mass center

    I.zero();
    for(size_t i = 0;i < ms.size();++ i)
    {
        Vector3d ri = vtx[i] - x0;
        I += Matrix3<double>(
                ms[i]*(M_SQR(ri.y)+M_SQR(ri.z)), -ms[i]*ri.x*ri.y, -ms[i]*ri.x*ri.z,
                -ms[i]*ri.y*ri.x, ms[i]*(M_SQR(ri.x)+M_SQR(ri.z)), -ms[i]*ri.y*ri.z,
                -ms[i]*ri.z*ri.x, -ms[i]*ri.z*ri.y, ms[i]*(M_SQR(ri.x)+M_SQR(ri.y)));
    }
    return mass;
}

/* 
 * Let the inertia matrix of the original object is T
 * then 
 *   => V'*T*V = D (D is diagonal matrix)
 *   => P*V'*T*V*P' = D_ellip (inertia matrix for the corresponding ellipsoid)
 *   P is the permutation matrix due to the eigenvalue permutation
 *   e.g. [1 2 3] => [2 3 1]
 *    P = [ 0 1 0 ;
 *          0 0 1 ;
 *          1 0 0 ]
 */
int main(int argc, char* argv[])
{
    if ( argc != 5 )
    {
        cerr << "Invalid arguments!" << endl;
        cerr << "Usage: " << argv[0] << " <tet mesh filename pattern> <start id> <end id> <output file> " << endl;
        exit(1);
    }

    tetfileptn = argv[1];
    startid = atoi(argv[2]);
    endid   = atoi(argv[3]);

    ofstream fout(argv[4]);
    if ( fout.fail() )
    {
        cerr << "Cannot open file: " << argv[4] << " for writing" << endl;
        exit(2);
    }
    fout << setprecision(12);

    int n = endid - startid + 1;
    fout << n << endl;
    char filename[128];
    Matrix3<double> inertia;
    Point3d mc;
    for(int i = startid;i <= endid;++ i)
    {
        fout << i << ' ';

        sprintf(filename, tetfileptn.c_str(), i);
        cout << "Reading tet mesh[#" << i << "]: " << filename << endl;
        TTetMesh tmesh;
        FV_TetMeshLoader_Double::load_mesh(filename, tmesh);
        tmesh.init();
        inertia_tensor(tmesh, inertia, mc);
        fout << mc.x << ' ' << mc.y << ' ' << mc.z << ' ';
        fout << inertia.cols[0].x << ' ' << inertia.cols[1].x << ' ' << inertia.cols[2].x << ' '
             << inertia.cols[0].y << ' ' << inertia.cols[1].y << ' ' << inertia.cols[2].y << ' '
             << inertia.cols[0].z << ' ' << inertia.cols[1].z << ' ' << inertia.cols[2].z << endl;
    }
    fout.close();

    return 0;
}

