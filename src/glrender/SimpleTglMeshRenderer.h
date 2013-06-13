/*
 * =====================================================================================
 *
 *       Filename:  SimpleTglMeshRenderer.h
 *
 *    Description:  Simple OPENGL renderer for triangle mesh
 *
 *        Version:  1.0
 *        Created:  08/24/10 23:01:12
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#ifndef SIMPLE_TGL_MESH_RENDERER_H
#   define SIMPLE_TGL_MESH_RENDERER_H

#include "Renderer.h"
#include "geometry/TriangleMesh.hpp"

class SimpleTglMeshRenderer : public Renderer
{
    public:
        SimpleTglMeshRenderer():pmesh_(NULL) { }

        SimpleTglMeshRenderer(const TriangleMesh<double>* mesh):
                pmesh_(mesh) { }

        void set_mesh(const TriangleMesh<double>* mesh)
        {   pmesh_ = mesh; }

        void render();

    private:
        const TriangleMesh<double>* pmesh_;
};

#endif
