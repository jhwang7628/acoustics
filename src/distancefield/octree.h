#ifndef __OCTREE_H__
#define __OCTREE_H__

#ifdef WIN32
  #include <windows.h>
#endif

//==============================================================================
//  Octree.h
//  Author: Kurt Miller
//  Modified: Jernej Barbic <barbic@cs.cmu.edu>
//
//  Copyright (C) 2003 Gradient Studios, Inc.  All Rights Reserved.
//==============================================================================
//
//  This class is a simple octree implementation.  Its presented for
//  illustrative purposes, but feel free to use or modify it however you'd like.
//
//==============================================================================

#include <vector>
#include <queue>
#include <string>
#include <map>
#include <fstream>
#include <iostream>
#include <algorithm>
using namespace std;

#include "boundingBox.h"
#include "triangle.h"
#include "sphere.h"
#include "minivector.h"

template<class TriangleClass> 
class Octree 
{
public:

    // make empty octree
    Octree(int maxDepth_g=10, int depth_g=0);

    ~Octree() { release(); }

    bool build(std::vector<TriangleClass> &tripool, 
               FieldBoundingBox &parentCube, int maxTriCount);

    // will construct a bounding box automatically:
    bool build(std::vector<TriangleClass> &tripool, int maxTriCount);
    void setBuildPrintInfo(int info); // controls what gets printed out during construction: 0 = nothing, 1 = up to five warnings, 2 = all warnings
    void getBuildInfo(int * numMaxDepthExceededCases, int * numMaxTriInDepthExceededCases);

    // note: these two routines might return the same colliding triangle several times
    // (call <TriangleClass>::makeUniqueList if necessary to have unique triangles)
    bool buildCollisionList(vector<TriangleClass*> &tlist, const Sphere &ps);
    bool buildCollisionList(vector<TriangleClass*> &tlist, Vec3d segmentStartPoint, Vec3d segmentEndPoint);

    int getDepth(); // compute tree depth

    FieldBoundingBox getBoundingBox() { return boundingBox; }

protected:

    void release(); // free the memory

    void createChildCubes(FieldBoundingBox *nodes);

    FieldBoundingBox  boundingBox;
    Octree      *children[8];

    int depth; // how deep in the hierarchy is this octree
    
    int maxDepth; // max depth allowed

    vector<TriangleClass> triangles;

    static int buildPrintInfo;
    static int numMaxDepthExceededCases;
    static int numMaxTriInDepthExceededCases;

};

template<class TriangleClass>
inline void Octree<TriangleClass>::getBuildInfo(int * numMaxDepthExceededCases, int * numMaxTriInDepthExceededCases)
{
  *numMaxDepthExceededCases = this->numMaxDepthExceededCases;
  *numMaxTriInDepthExceededCases = this->numMaxTriInDepthExceededCases;
}

template<class TriangleClass>
inline void Octree<TriangleClass>::setBuildPrintInfo(int info)
{
  buildPrintInfo = info;
  if (info == 1)
  {
    numMaxDepthExceededCases = 0;
    numMaxTriInDepthExceededCases = 0;
  }
}


#endif
