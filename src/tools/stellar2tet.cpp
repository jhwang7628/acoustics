#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include "utils/print_msg.h"
#include "io/StellarIO.hpp"
#include "io/TetMeshWriter.hpp"

using namespace std;

namespace fs = boost::filesystem;

static bool checkUselessVtx = false;
static FixVtxTetMesh<double>    mesh;
static FixVtxTetMesh<double>*   outMsh = NULL;
static bool* vtxUsed = NULL;

static void usage(char* cmd)
{
    cout << "Usage: " << fs::path(cmd).filename() 
         << "[-c] <input mesh(.node)> <output mesh(.tet)>" << endl;
}

static void parse_cmd(int argc, char* argv[])
{
    int opt;
    while ( (opt = getopt(argc, argv, "hc")) != -1 )
    {
        switch (opt)
        {
            case 'h':
                usage(argv[0]);
                exit(0);
            case 'c':
                checkUselessVtx = true;
                break;
        }
    }

    if ( argc - optind != 2 ) 
    {
        PRINT_ERROR("No input/output mesh is given\n");
        usage(argv[0]);
        exit(1);
    }
}

static bool check_useless_vtx()
{
    vtxUsed = new bool[ mesh.num_vertices() ];
    memset(vtxUsed, false, sizeof(bool)*mesh.num_vertices());

    const vector<TetMesh<double>::TetIdx>& idx = mesh.tet_indices();
    for(int i = 0;i < idx.size();++ i)
    {
        vtxUsed[idx[i][0]] = true;
        vtxUsed[idx[i][1]] = true;
        vtxUsed[idx[i][2]] = true;
        vtxUsed[idx[i][3]] = true;
    }

    for(int i = 0;i < mesh.num_vertices();++ i)
        if ( !vtxUsed[i] ) return true;
    return false;
}

static void regularize_mesh()
{
    outMsh = new FixVtxTetMesh<double>;
    const int NUM_VTX = mesh.num_vertices();
    const vector<Point3d>& vtx = mesh.vertices();
    int* idmap = new int[NUM_VTX];
    int usedVtxCnt = 0;

    for(int i = 0;i < NUM_VTX;++ i)
    {
        if ( vtxUsed[i] )
        {
            idmap[i] = usedVtxCnt ++;
            outMsh->add_vertex(vtx[i]);
        }
    }

    const vector<TetMesh<double>::TetIdx>& idx = mesh.tet_indices();
    for(int i = 0;i < idx.size();++ i)
    {
        outMsh->add_tet(idmap[idx[i][0]], idmap[idx[i][1]],
                        idmap[idx[i][2]], idmap[idx[i][3]]);
    }

    delete []idmap;
}

int main(int argc, char* argv[])
{
    parse_cmd(argc, argv);
    
    string infile(argv[optind]);
    if ( !boost::iends_with(infile, ".node") )
    {
        PRINT_ERROR("Input filename with .node extension is expected\n");
        return 1;
    }

    string prefix = infile.substr(0, infile.length()-5);
    // load mesh
    StellarTetMeshLoader::load_mesh(prefix.c_str(), &mesh);

    if ( checkUselessVtx && check_useless_vtx() )
        regularize_mesh();

    if ( outMsh )
    {
        FV_TetMeshWriter_Double::write_mesh(argv[optind+1], *outMsh);
        delete outMsh;
    }
    else
        FV_TetMeshWriter_Double::write_mesh(argv[optind+1], mesh);
    delete[] vtxUsed;
    return 0;
}

