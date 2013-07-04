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
#include <fstream>
#include <iostream>

#include "boundingBox.h"
#include "sphere.h"
#include "octree.h"

//==============================================================================
//  Octree::Octree()
//==============================================================================

template<class TriangleClass>
Octree<TriangleClass>::Octree(int maxDepth_g, int depth_g) 
{
  depth = depth_g; maxDepth = maxDepth_g;

  for(int i=0; i<8; i++)
    children[i] = NULL;
}


//==============================================================================
//  Octree::build()
//==============================================================================

template<class TriangleClass>
bool Octree<TriangleClass>::build(std::vector<TriangleClass> &tripool, int maxTriCount)
{
  FieldBoundingBox parentCube(tripool);
  parentCube.expand(1.05);
  return build(tripool,parentCube,maxTriCount);
}

template<class TriangleClass>
int Octree<TriangleClass>::buildPrintInfo = 0;

template<class TriangleClass>
int Octree<TriangleClass>::numMaxDepthExceededCases = 0;

template<class TriangleClass>
int Octree<TriangleClass>::numMaxTriInDepthExceededCases = 0;


template<class TriangleClass>
bool Octree<TriangleClass>::build(std::vector<TriangleClass> &tripool, 
                   FieldBoundingBox &parentCube, int maxTriCount)
{
        int i;
        //static int numMaxTriInDepthExceededCases;
        //static int numMaxDepthExceededCases; 
        // set the node's bounding box;
        boundingBox = parentCube;

        // total number of triangles in this octree
        int tricount = (int)tripool.size();

        // if there are fewer triangles than the threshold value, or if max depth has been reached
        // this node becomes a leaf which stores the triangles;
	// max depth checking is necessary, since otherwise vertices with high valence will always be contained in some voxel,
        //  and that voxel will be split forever
        if ((tricount <= maxTriCount) || (depth >= maxDepth))
        {
	    //cout << "L" << tricount << " " << depth << " ";
            for(int i=0; i<tricount; i++)
            {
                triangles.push_back(tripool[i]);
            }

            if ((depth >= maxDepth) && (tricount >= maxTriCount))
            {
              if (buildPrintInfo == 2)
                printf("Warning: cube max depth (%d) reached with more triangles per cell (%d) than requested (%d)\n",maxDepth,tricount,maxTriCount);

              if (buildPrintInfo == 1)
              {
                numMaxDepthExceededCases++;

                if (numMaxDepthExceededCases <= 5)
                  printf("Warning: cube max depth (%d) reached with more triangles per cell (%d) than requested (%d)\n",maxDepth,tricount,maxTriCount);

                if (numMaxDepthExceededCases == 5)
                {
                  printf("(suppressing further max depth warning messages)\n");
                }

                if (tricount > numMaxTriInDepthExceededCases)
                  numMaxTriInDepthExceededCases = tricount;
              }
            }

            return true;
        }

	//cout << "N" << tricount << " ";

        // create the child cubes;
        FieldBoundingBox childCube[8];
        createChildCubes(childCube);

        /*
        cout << "Created child bounding boxes." << endl;
        cout << "Parent: " << boundingBox.bmin() << " " << boundingBox.bmax() << " " << boundingBox.diameter() << endl;
        cout << "Children:" << endl;
        for(int j=0; j<8; j++)
          cout << childCube[j].bmin() << " " << childCube[j].bmax() << " " << childCube[j].diameter() << endl;
        */       

        std::vector<TriangleClass> childTriangles[8];

        // for each child, generate a list of triangles that
        // intersect the child cube;

        for(i=0; i<tricount; i++)
        {
            for(int j=0; j<8; j++)
            {
                if (tripool[i].doesIntersectBox(childCube[j]))
                {
                    childTriangles[j].push_back(tripool[i]);
                }
            }
        }
/* 
        for(int j=0; j<8; j++)
        {
	  cout << "Child " << j << " has " << childTriangles[j].size() << " triangles:" << endl;
          cout << "Child bounding box:" << endl;
          cout << childCube[j].bmin() << " " << childCube[j].bmax() << endl;
          for (unsigned int i=0; i<childTriangles[j].size(); i++)
          {
            cout << childTriangles[j][i].first() << " " << childTriangles[j][i].second() << " " << childTriangles[j][i].third();
            cout << endl;
          }
        }
*/
        //for (int i=0; i<8; i++)
	//  cout << "S" << i << ":" << childTriangles[i].size() << " ";


        // for any children with intersecting triangles, create and
        // recursively build the node;

        for(i=0; i<8; i++)
        {
            if (childTriangles[i].size() > 0)
            {
                children[i] = new Octree<TriangleClass>(maxDepth,depth+1);

                if(children[i] == NULL)
                {
                    return false;
                }

                if(!children[i]->build(childTriangles[i], childCube[i], maxTriCount))
                {
                    return false;
                }
            }
            else
              children[i] = NULL;
        }

        //cout << "C ";

        return true;
}

template<class TriangleClass>
int Octree<TriangleClass>::getDepth()
{
  if(triangles.size() > 0)
    return 0; // leaf
  else
  {
    int maxChildDepth=0;
    for(int i=0; i<8; i++)
    {
      if(children[i] != NULL)
      {
        int childDepth = children[i]->getDepth();
        if (childDepth > maxChildDepth)
          maxChildDepth = childDepth;
      }
    }
    return maxChildDepth + 1;
  }
}

//==============================================================================
//  Octree::buildCollisionList()
//==============================================================================

template<class TriangleClass>
bool Octree<TriangleClass>::buildCollisionList(std::vector<TriangleClass*> &tlist, const Sphere &ps)
    {
        // if the bounding sphere does not intersect this node,
        // we can disregard it and all of its children;

        if(!ps.doesBoundingBoxIntersect(boundingBox))
        {
            return true;
        }

        // if this node is a leaf node with triangles, add the children
        // to the potential collision list;  otherwise recursively check;

        if(triangles.size() > 0)
        {
            int tsize = (int)triangles.size();

            for(int i=0; i<tsize; i++)
            {
                tlist.push_back(&(triangles[i]));
            }
        }
        else
        {
            for(int i=0; i<8; i++)
            {
                if(children[i] != NULL)
                {
                    children[i]->buildCollisionList(tlist, ps);
                }
            }
        }

        return true;
    }


template<class TriangleClass>
bool Octree<TriangleClass>::buildCollisionList(std::vector<TriangleClass*> &tlist,  Vec3d segmentStartPoint, Vec3d segmentEndPoint)
    {
        // if the bounding sphere does not intersect this node,
        // we can disregard the node and all of its children;

        Vec3d intersectionPoint;
        if(!boundingBox.lineSegmentIntersection(segmentStartPoint, segmentEndPoint,&intersectionPoint))
        {
            return true;
        }

        // if this node is a leaf node with triangles, add the children
        // to the potential collision list;  otherwise recursively check;

        if(triangles.size() > 0)
        {
            int tsize = (int)triangles.size();

            for(int i=0; i<tsize; i++)
            {
                Vec3d intersectionPoint;
                int collisionStatus = triangles[i].lineSegmentIntersection(segmentStartPoint,segmentEndPoint,&intersectionPoint);
                if (collisionStatus == 1)
                  tlist.push_back(&(triangles[i]));
            }
        }
        else
        {
            for(int i=0; i<8; i++)
            {
                if(children[i] != NULL)
                {
                    children[i]->buildCollisionList(tlist, segmentStartPoint, segmentEndPoint);
                }
            }
        }

        return true;
    }

//==============================================================================
//  Octree::release()
//==============================================================================

template<class TriangleClass>
void Octree<TriangleClass>::release()
    {
        for(int i=0; i<8; i++)
        {
            if(children[i] != NULL)
            {
                children[i]->release();

                delete(children[i]);
                children[i] = NULL;
            }
        }

        triangles.clear();
    }

//==============================================================================
//  Octree::createChildCubes()
//==============================================================================

template<class TriangleClass>
void Octree<TriangleClass>::createChildCubes(FieldBoundingBox *childCubes)
    {
        Vec3d center = boundingBox.center();

        // top left, near;
        childCubes[0].setbmin(boundingBox.bmin()[0], center[1], boundingBox.bmin()[2]);
        childCubes[0].setbmax(center[0], boundingBox.bmax()[1], center[2]);

        // top right, near;
        childCubes[1].setbmin(center[0], center[1], boundingBox.bmin()[2]);
        childCubes[1].setbmax(boundingBox.bmax()[0], boundingBox.bmax()[1], center[2]);

        // bottom left, near;
        childCubes[2].setbmin(boundingBox.bmin()[0], boundingBox.bmin()[1], boundingBox.bmin()[2]);
        childCubes[2].setbmax(center[0], center[1], center[2]);

        // bottom right, near;
        childCubes[3].setbmin(center[0], boundingBox.bmin()[1], boundingBox.bmin()[2]);
        childCubes[3].setbmax(boundingBox.bmax()[0], center[1], center[2]);

        // top left, far;
        childCubes[4].setbmin(boundingBox.bmin()[0], center[1], center[2]);
        childCubes[4].setbmax(center[0], boundingBox.bmax()[1], boundingBox.bmax()[2]);

        // top right, far;
        childCubes[5].setbmin(center[0], center[1], center[2]);
        childCubes[5].setbmax(boundingBox.bmax()[0], boundingBox.bmax()[1], boundingBox.bmax()[2]);

        // bottom left, far;
        childCubes[6].setbmin(boundingBox.bmin()[0], boundingBox.bmin()[1], center[2]);
        childCubes[6].setbmax(center[0], center[1], boundingBox.bmax()[2]);

        // bottom right, far;
        childCubes[7].setbmin(center[0], boundingBox.bmin()[1], center[2]);
        childCubes[7].setbmax(boundingBox.bmax()[0], center[1], boundingBox.bmax()[2]);

        // sanity check
	for (int i=0; i<7; i++)
          childCubes[i].verifyBox();
    }

template Octree<TriangleBasic>::Octree(int maxDepth_g, int depth_g);
template bool Octree<TriangleBasic>::build(std::vector<TriangleBasic> &tripool, FieldBoundingBox &parentCube, int maxTriCount);
template bool Octree<TriangleBasic>::build(std::vector<TriangleBasic> &tripool, int maxTriCount);
template int Octree<TriangleBasic>::getDepth();
template bool Octree<TriangleBasic>::buildCollisionList(std::vector<TriangleBasic*> &tlist, const Sphere &ps);
template bool Octree<TriangleBasic>::buildCollisionList(std::vector<TriangleBasic*> &tlist, Vec3d segmentStartPoint, Vec3d segmentEndPoint);

template void Octree<TriangleBasic>::release();
template void Octree<TriangleBasic>::createChildCubes(FieldBoundingBox *childCubes);

template Octree<TriangleWithCollisionInfo>::Octree(int maxDepth_g, int depth_g);
template bool Octree<TriangleWithCollisionInfo>::build(std::vector<TriangleWithCollisionInfo> &tripool, FieldBoundingBox &parentCube, int maxTriCount);
template bool Octree<TriangleWithCollisionInfo>::build(std::vector<TriangleWithCollisionInfo> &tripool, int maxTriCount);
template bool Octree<TriangleWithCollisionInfo>::buildCollisionList(std::vector<TriangleWithCollisionInfo*> &tlist, const Sphere &ps);
template bool Octree<TriangleWithCollisionInfo>::buildCollisionList(std::vector<TriangleWithCollisionInfo*> &tlist, Vec3d segmentStartPoint, Vec3d segmentEndPoint);
template int Octree<TriangleWithCollisionInfo>::getDepth();
template void Octree<TriangleWithCollisionInfo>::release();
template void Octree<TriangleWithCollisionInfo>::createChildCubes(FieldBoundingBox *childCubes);

template Octree<TriangleWithCollisionInfoAndPseudoNormals>::Octree(int maxDepth_g, int depth_g);
template bool Octree<TriangleWithCollisionInfoAndPseudoNormals>::build(std::vector<TriangleWithCollisionInfoAndPseudoNormals> &tripool, FieldBoundingBox &parentCube, int maxTriCount);
template bool Octree<TriangleWithCollisionInfoAndPseudoNormals>::build(std::vector<TriangleWithCollisionInfoAndPseudoNormals> &tripool, int maxTriCount);
template bool Octree<TriangleWithCollisionInfoAndPseudoNormals>::buildCollisionList(std::vector<TriangleWithCollisionInfoAndPseudoNormals*> &tlist, const Sphere &ps);
template bool Octree<TriangleWithCollisionInfoAndPseudoNormals>::buildCollisionList(std::vector<TriangleWithCollisionInfoAndPseudoNormals*> &tlist, Vec3d segmentStartPoint, Vec3d segmentEndPoint);
template int Octree<TriangleWithCollisionInfoAndPseudoNormals>::getDepth();
template void Octree<TriangleWithCollisionInfoAndPseudoNormals>::release();
template void Octree<TriangleWithCollisionInfoAndPseudoNormals>::createChildCubes(FieldBoundingBox *childCubes);

//==============================--end-of-file--=================================
