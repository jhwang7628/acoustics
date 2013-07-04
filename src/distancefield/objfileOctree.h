#ifndef __OBJFILEOCTREE_H__
#define __OBJFILEOCTREE_H__

#ifdef WIN32
#include <windows.h>
#endif

//==============================================================================
//  objFileOctree.h
//  Author: Jernej Barbic <barbic@cs.cmu.edu>
//
//  builds an octree on top of the geometry from a given obj file
//
//==============================================================================

#include "octree.h"
#include "triangle.h"
#include "boundingBox.h"
#include "sphere.h"
#include "objfile.h"
#include "minivector.h"

#if 0
class TriMesh;
#endif

template<class TriangleClass>
class ObjFileOctree 
{
    public:
        ObjFileOctree( ObjFile * objFile, int maxTriCount_g, int maxDepth_g, int printInfo = 1 );

        ObjFileOctree( const std::string& filename, int maxTriCount_g, int maxDepth_g,
                       int printInfo = 1 );

#if 0
        // Build an octree using a triangle mesh (which was originally generated
        // via an obj file)
        ObjFileOctree( TriMesh *mesh, int maxTriCount_g, int maxDepth_g, int printInfo = 1 );
#endif

        ~ObjFileOctree() { delete(root); }

        // range query
        // retrieves all triangles intersected by the sphere
        // return value is meaningless (always true)
        bool rangeQuery(vector< TriangleClass* > &tlist, const Sphere &ps)
        { return root->buildCollisionList(tlist,ps); }

        // line-triangles query
        // retrieves all triangles intersected by the line segment (could potentially return same triangle several times)
        // return value is meaningless (always true)
        bool lineSegmentIntersection(vector< TriangleClass* > &tlist, Vec3d segmentStart, Vec3d segmentEnd)
        { return root->buildCollisionList(tlist,segmentStart, segmentEnd); }

        void renderOctree() { root->render(); }
        void renderOctree(int level) { root->render(level); }
        void renderOctree(int level, int boxIndex) { root->render(level,boxIndex); }
        void setPrintOctreeInfo(int info) { root->setRenderInfo(info); }

        FieldBoundingBox boundingBox() { return (root->getBoundingBox()) ; }

    protected:
        vector<TriangleClass> triangles;

        Octree<TriangleClass> * root;
        int maxTriCount; // the max number of triangles in a cell
        int maxDepth;

        static const double bboxExpansionRatio;
};

#endif
