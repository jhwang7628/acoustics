/*
 * =====================================================================================
 *
 *       Filename:  SimpleTglMeshRenderer.cpp
 *
 *        Version:  1.0
 *        Created:  08/24/10 23:33:55
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#include "SimpleTglMeshRenderer.h"
#include <GL/gl.h>
#include <vector>

using namespace std;

void SimpleTglMeshRenderer::render()
{
    if ( !pmesh_ ) return;

    const vector<Point3d>&  vtx = pmesh_->vertices();
    const vector<Tuple3ui>& tgl = pmesh_->triangles();
    const vector<Vector3d>& nml = pmesh_->normals();

    glEnable(GL_LIGHTING);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glVertexPointer(3, GL_DOUBLE, 0, (const GLvoid*)&vtx[0]);
    glNormalPointer(GL_DOUBLE, 0, (const GLvoid*)&nml[0]);
    glDrawElements(GL_TRIANGLES, tgl.size()*3, GL_UNSIGNED_INT,
                   (const GLvoid*)&tgl[0]);
}
