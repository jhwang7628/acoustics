#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "geometry/FixVtxTetMesh.hpp"
#include "io/TetMeshReader.hpp"

using namespace std;

typedef FixVtxTetMesh<double>       TTetMesh;

static string tetfileptn;
static int startid, endid;

void writeMassCenter( Point3d x0, const char *filePrefix )
{
  char                       fileName[ 1024 ];
  double                     value;

  sprintf( fileName, "%s_centerOfMass.3vector", filePrefix );

  FILE                      *file = fopen( fileName, "wb" );
  
  if ( !file )
  {
    cerr << "Error: could not open " << fileName << " for writing" << endl;
    exit(1);
  }

  value = x0.x;
  fwrite( (void *)&value, sizeof( double ), 1, file );
  value = x0.y;
  fwrite( (void *)&value, sizeof( double ), 1, file );
  value = x0.z;
  fwrite( (void *)&value, sizeof( double ), 1, file );

  fclose( file );
}

void writeMass( double mass, const char *filePrefix )
{
  char                       fileName[ 1024 ];
  int                        size = 1;

  sprintf( fileName, "%s_mass.vector", filePrefix );

  FILE                      *file = fopen( fileName, "wb" );

  if ( !file )
  {
    cerr << "Error: could not open " << fileName << " for writing" << endl;
    exit(1);
  }

  // Write in a format that can be interpreted as a vector of
  // doubles
  fwrite( (void *)&size, sizeof( int ), 1, file );
  fwrite( (void *)&mass, sizeof( double ), 1, file );

  fclose( file );
}

int main(int argc, char* argv[])
{
#if 0
    if ( argc != 5 )
    {
        cerr << "Invalid arguments!" << endl;
        cerr << "Usage: " << argv[0] << " <tet mesh filename pattern> "
             << "<start id> <end id> <output file>" << endl;
        exit(1);
    }
#endif
    if ( argc != 3 )
    {
        cerr << "Invalid arguments!" << endl;
        cerr << "Usage: " << argv[0] << " <tet mesh filename > "
             << "<output file>" << endl;
        exit(1);
    }

    tetfileptn = argv[1];
    string outputFilePrefix = argv[2];
#if 0
    startid = atoi(argv[2]);
    endid   = atoi(argv[3]);
#endif

#if 0
    ofstream fout(argv[4]);
    if ( fout.fail() )
    {
        cerr << "Cannot open file: " << argv[4] << " for writing" << endl;
        exit(2);
    }

    char filename[128];
#endif
    Point3d x0;
    double mass;
#if 0
    fout << setprecision(22);
    for(int i = startid; i <= endid; ++ i)
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
#endif
    TTetMesh tmesh;
    FV_TetMeshLoader_Double::load_mesh( tetfileptn.c_str(), tmesh );
    tmesh.init();

    x0.zero();
    mass = 0.0;
    const std::vector<double>&    ms  = tmesh.masses();
    const std::vector< Point3d >& vtx = tmesh.vertices();

    for ( size_t j = 0; j < ms.size(); j++ )
    {
      mass += ms[ j ];
      x0.scaleAdd( ms[ j ], vtx[ j ] );
    }
    x0 /= mass;

#if 0
    // FIXME:
    mass *= 6983.0;
#endif

    std::cout << "mass center = " << x0 << std::endl;
    std::cout << "mass = " << mass << std::endl;

    writeMassCenter( x0, outputFilePrefix.c_str() );
    writeMass( mass, outputFilePrefix.c_str() );

    return 0;
}
