/*
 * Read the mesh and merge the vertices who are on boundary and are close to each other.
 * The algorithm is as follows
 * 1. Read the mesh
 * 2. Check the vertices who are on boundary
 * 3. For each vertex on boundary, try to find another close vertex who is also 
 *    on boundary, and merge them
 */
#include <iostream>
#include <stdio.h>
#include <string>
#include <set>
#include <map>

#include "utils/print_msg.h"
#include "io/TglMeshReader.hpp"
#include "io/TglMeshWriter.hpp"

using namespace std;

static TriangleMesh<double> mesh;
static string               outfile;
static double               distThresh;
static set<int>             bdvtx;  // boundary vertex
static vector< vector<int> >   vtxtgl;   // the triangles associated with each vertex

static void usage(const char* cmd)
{
    printf("Usage: %s -d <dist threshold> -i <input mesh> -o <output mesh>\n", cmd);
}

static void parse_cmd(int argc, char* argv[])
{
    int opt;
    while ( (opt = getopt(argc, argv, "hd:i:o:")) != -1 )
    {
        switch (opt)
        {
            case 'h':
                usage(argv[0]);
                exit(0);
            case 'd':
                distThresh = atof(optarg);
                break;
            case 'i':
                MeshObjReader::read(optarg, mesh);
                break;
            case 'o':
                outfile = optarg;
                break;
        }
    }
}

static void detect_boundary_edges()
{
    bdvtx.clear();
    vtxtgl.resize(mesh.num_vertices());
    for(size_t i = 0;i < vtxtgl.size();++ i) vtxtgl[i].clear();

    const vector<Tuple3ui>& tgl = mesh.surface_indices();
    map< int, map<int, int> >  edgecnt;
    for(size_t i = 0;i < tgl.size();++ i)
    for(int j = 0;j < 3;++ j)
    {
        int v1 = tgl[i][j], v2 = tgl[i][(j+1)%3];
        if ( v1 > v2 ) swap(v1, v2);
        ++ edgecnt[v1][v2];
        vtxtgl[tgl[i][j]].push_back(i);
    }

    const map< int, map<int, int> >::iterator end = edgecnt.end();
    for(map< int, map<int, int> >::iterator it1 = edgecnt.begin();
            it1 != end;++ it1)
    {
        const map<int, int>::iterator end2 = it1->second.end();
        for(map<int,int>::iterator it2 = it1->second.begin(); 
                it2 != end2;++ it2)
            if ( it2->second == 1 )
            {
                bdvtx.insert(it1->first);
                bdvtx.insert(it2->first);
            }
            else if ( it2->second > 2 )
            {
                PRINT_ERROR("The given mesh is not a manifold\n");
                exit(1);
            }
    }
    printf("%d boundary vertices detected\n", (int)bdvtx.size());
}

static void merge_boundary_vertices()
{
    double distnow;
    vector<Point3d>& vtx  = mesh.vertices();
    vector<Tuple3ui>& tgl = mesh.triangles();
    while ( !bdvtx.empty() )
    {
        int cvtx = *(bdvtx.begin());

        //// find the closest boundary vertex
        double mindist = distThresh;
        int    minid = -1;
        const set<int>::iterator end = bdvtx.end();
        for(set<int>::iterator it = bdvtx.begin();it != end;++ it)
            if ( *it != cvtx && 
                (distnow = vtx[*it].distanceSqr(vtx[cvtx])) < mindist )
            {
                mindist = distnow;
                minid   = *it;
            }

        if ( minid < 0 )
        {
            PRINT_WARNING("Cannot find closest boundary vertex. Skip it. vtxid=%d\n", cvtx);
            bdvtx.erase(cvtx);
            continue;
        }

        //// merge it
        vtx[cvtx] = (vtx[cvtx] + vtx[minid])*0.5;
        for(size_t i = 0;i < vtxtgl[minid].size();++ i)
        {
            for(int j = 0;j < 3;++ j)
                if ( tgl[vtxtgl[minid][i]][j] == minid )
                {
                    tgl[vtxtgl[minid][i]][j] = cvtx;
                    break;
                }
        }
        bdvtx.erase(cvtx);
        bdvtx.erase(minid);
    }
    printf("Finished Merging vertices.\n");
}

int main(int argc, char* argv[])
{
    parse_cmd(argc, argv);  // parse command line
    detect_boundary_edges();
    if ( !bdvtx.empty() )
    {
        merge_boundary_vertices();
        //// output mesh
        MeshObjWriter::write(mesh, outfile.c_str());
    }
    return 0;
}
