#ifndef __TRIANGLE_H__
#define __TRIANGLE_H__

#ifdef WIN32
  #include <windows.h>
#endif

//==============================================================================
//  triangle.h
//  Author: Jernej Barbic <barbic@cs.cmu.edu>
//
//==============================================================================

#include <vector>
#include <queue>
#include <string>
#include <map>
#include <fstream>
#include <iostream>
#include <algorithm>
#include "math.h"
using namespace std;

#include "boundingBox.h"
#include "minivector.h"

class TriangleBasic
{
public:

  TriangleBasic(Vec3d first_g, Vec3d second_g, Vec3d third_g): 
	first_(first_g), second_(second_g), third_(third_g), index_(0) {}

  // accessors
  inline Vec3d first() {return first_ ;}
  inline Vec3d second() {return second_ ;}
  inline Vec3d third() {return third_ ;}
  inline int index() { return index_; }

	inline int firstIndex() { return firstIndex_; }
	inline int secondIndex() { return secondIndex_; }
	inline int thirdIndex() { return thirdIndex_; }

  // index can be used to keep track of what triangle this is
  inline void setIndex(int index_) { this->index_ = index_; }

  // squared 3d distance to a point
  double distanceToPoint2(Vec3d point) { cout << "Unimplemented..." << endl; return 1;} // unimplemented
  double distanceToPoint(Vec3d point) { return sqrt(distanceToPoint2(point));}
  
  bool doesIntersectBox(FieldBoundingBox & bbox);

  int lineSegmentIntersection(Vec3d segmentStart, Vec3d segmentEnd, Vec3d * intersectionPoint);
    //    Output: intersection point (when it exists)
    //    Return: -1 = triangle is degenerate (a segment or point)
    //             0 = disjoint (no intersect)
    //             1 = intersect in unique point I1
    //             2 = are in the same plane


  Vec3d getBarycentricLocation(double alpha, double beta, double gamma) { return alpha * first_ + beta * second_ + gamma * third_; }

	// Set the indices of the vertices in some external array
	void setVertexIndex(int firstIndex, int secondIndex, int thirdIndex)
	{
		firstIndex_ = firstIndex;
		secondIndex_ = secondIndex;
		thirdIndex_ = thirdIndex;
	}

protected:

  Vec3d first_, second_, third_;
  int index_;

	// Indeces of the three vertices in some
	// global array
	int firstIndex_, secondIndex_, thirdIndex_;
};

// makes the triangle list unique, using the "index_" field of the <TriangleClass>
template<class TriangleClass>
void makeUniqueList(vector<TriangleClass*> &tlist, vector<TriangleClass*> & uniqueList);

template<class TriangleClass>
void makeUniqueList(vector<TriangleClass*> &tlist); // overwrites tlist with the unique list

class TriangleWithCollisionInfo : public TriangleBasic
{
public:
  TriangleWithCollisionInfo(Vec3d first_g, Vec3d second_g, Vec3d third_g):
        TriangleBasic(first_g,second_g,third_g) { ComputeCollisionInfo(); }

  // squared 3d distance to a point
  // also returns the closest feature to the query point:
  //  0: vertex0
  //  1: vertex1
  //  2: vertex2
  //  3: edge among 01
  //  4: edge among 12
  //  5: edge among 20
  //  6: the face itself
  double distanceToPoint2(Vec3d point, int * closestFeature);
  double distanceToPoint2(Vec3d point, int * closestFeature, double * alpha, double * beta, double * gamma); // also returns the barycentric coordinates of the closest point

protected:

  // note: the following collision detection parameters are pre-computed with respect to a permuted set of triangle indices (a cycle)
    Mat3d Q;
    Vec3d x0;

    double sidea,sideb,sidec,area;
    Vec3d S1,S2,N11,N12,N21,N22;

  int permutation[6]; // used internally so that first triangle edge is always the longest
  // transformation(x) equals Q * x + x0;

  void ComputeCollisionInfo();
};


class TriangleWithCollisionInfoAndPseudoNormals : public TriangleWithCollisionInfo
{
public:
  // pseudoNormals is a caller-provided list of 7 normals, indexed the same as closest features below
  TriangleWithCollisionInfoAndPseudoNormals(Vec3d first_g, Vec3d second_g, Vec3d third_g, Vec3d * pseudoNormals);

  inline Vec3d pseudoNormal(int closestFeature) { return pseudoNormal_[closestFeature]; }
  inline Vec3d interpolatedVertexPseudoNormal(double alpha, double beta, double gamma) { return norm(alpha * pseudoNormal_[0] + beta * pseudoNormal_[1] + gamma * pseudoNormal_[2]); }

  // for vertices, returns the vertex itself
  // for edges, returns the midpoint of the edge (for the purposes of distance sign test, this could be any point on the edge)
  // for faces, returns the face centroid
  inline Vec3d pseudoClosestPosition(int closestFeature) { return pseudoClosestPosition_[closestFeature]; }

protected:

  Vec3d pseudoNormal_[7];
  Vec3d pseudoClosestPosition_[7];

};

#endif
