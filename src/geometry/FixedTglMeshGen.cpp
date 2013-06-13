/*
 * =====================================================================================
 *
 *       Filename:  FixedTglMeshGen.cpp
 *
 *    Description:  Implement FixedTglMeshGen
 *
 *        Version:  1.0
 *        Created:  10/27/10 10:18:55
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#include "FixedTglMeshGen.h"
#include "utils/macros.h"
#include "utils/print_msg.h"
#include <vector>

using namespace std;

static int next_free_vtx(int freevid, int NVTX, const bool* isFixed)
{
    for(;freevid < NVTX && isFixed[freevid];++ freevid);

    if ( freevid >= NVTX )
    {
        PRINT_ERROR("Too many fixed vertices!?\n");
        SHOULD_NEVER_HAPPEN(2);
    }
    return freevid;
}

/*
 * when return, vidMap is a symmetric mapping, which maps the vertex ID from
 * the input mesh to the output mesh
 */
void gen_fixed_tgl_mesh(const TriangleMesh<double>* inMsh, const bool* isFixed,
                        TriangleMesh<double>* outMsh, int* vidMap)
{
    const int NVTX = inMsh->num_vertices();
    int freevid = 0, fixvid = NVTX - 1;
    for(int i = 0;i < NVTX;++ i) vidMap[i] = i;
    for(;;) {
        for(;fixvid >= 0 && !isFixed[fixvid];-- fixvid);
        if ( fixvid < 0 ) break;
        freevid = next_free_vtx(freevid, NVTX, isFixed);
        if ( freevid >= fixvid ) break;

        vidMap[freevid] = fixvid;
        vidMap[fixvid]  = freevid;

        ++ freevid;
        -- fixvid;
    }

    outMsh->clear();
    const vector<Point3d>& vtx = inMsh->vertices();
    const vector<Tuple3ui>& tgl = inMsh->triangles();
    for(int i = 0;i < NVTX;++ i)
        outMsh->add_vertex(vtx[vidMap[i]]);
    for(int i = 0;i < (int)tgl.size();++ i)
        outMsh->add_triangle(vidMap[tgl[i][0]], vidMap[tgl[i][1]], vidMap[tgl[i][2]]);
}
