#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include <geometry/FixVtxTetMesh.hpp>
#include <io/TetMeshReader.hpp>

using namespace std;

typedef FixVtxTetMesh<double>       TTetMesh;

static string tetfileptn;
static int startid, endid;

int main(int argc, char* argv[])
{
    if ( argc != 5 )
    {
        cerr << "Invalid arguments!" << endl;
        cerr << "Usage: " << argv[0] << " <tet mesh filename pattern> <start id> <end id>"
                " <output file>" << endl;
        exit(1);
    }
    tetfileptn = argv[1];
    startid = atoi(argv[2]);
    endid   = atoi(argv[3]);

    ofstream fout(argv[4], ios::app);
    if ( fout.fail() )
    {
        cerr << "Cannot open file: " << argv[4] << " for writing" << endl;
        exit(2);
    }

    char filename[128];
    Point3d x0;
    double mass;
    fout << setprecision(22);
    for(int i = startid;i <= endid;++ i)
    {
        sprintf(filename, tetfileptn.c_str(), i);
        cout << "Reading tet mesh[#" << i << "]: " << filename << endl;
        TTetMesh tmesh;
        FV_TetMeshLoader_Double::load_mesh(filename, tmesh);
        tmesh.init();

        cerr << i << "  " << tmesh.num_fixed_vertices() << endl;
        const std::vector<double>&    ms  = tmesh.masses();
        const std::vector< Point3d >& vtx = tmesh.vertices();
        cout << "Computing mass center [#" << i << ']' << endl;
        x0.zero();  mass = 0;
        for(size_t j = 0;j < ms.size();++ j)
        {
            mass += ms[j];
            x0.scaleAdd(ms[j], vtx[j]);
        }
        x0 /= mass;
        fout << i << ' ' << x0.x << ' ' << x0.y << ' ' << x0.z << endl;
    }

    return 0;
}
