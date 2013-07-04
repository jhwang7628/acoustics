//==============================================================================
//  FieldBoundingBox.h
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
using namespace std;

#include "float.h"

#include "boundingBox.h"
#include "triangle.h"

#include "assert.h"

template<class Triangle> 
FieldBoundingBox::FieldBoundingBox(vector<Triangle> & tripool)
{
  // set bmin_, bmax_

  bmin_ =  Vec3d(+DBL_MAX,+DBL_MAX,+DBL_MAX);
  bmax_ =  Vec3d(-DBL_MAX,-DBL_MAX,-DBL_MAX);
                                                                                                                                                             
  for(unsigned int i=0; i < tripool.size(); i++) // over all vertices
  {
    Vec3d p;

    p = tripool[i].first();
                                                                                                                                                             
    if (p[0] < bmin_[0])
      bmin_[0] = p[0];
    if (p[0] > bmax_[0])
      bmax_[0] = p[0];
                                                                                                                                                             
    if (p[1] < bmin_[1])
      bmin_[1] = p[1];
    if (p[1] > bmax_[1])
      bmax_[1] = p[1];
                                                                                                                                                             
    if (p[2] < bmin_[2])
      bmin_[2] = p[2];
    if (p[2] > bmax_[2])
      bmax_[2] = p[2];

    p = tripool[i].second();
                                                                                                                                                             
    if (p[0] < bmin_[0])
      bmin_[0] = p[0];
    if (p[0] > bmax_[0])
      bmax_[0] = p[0];
                                                                                                                                                             
    if (p[1] < bmin_[1])
      bmin_[1] = p[1];
    if (p[1] > bmax_[1])
      bmax_[1] = p[1];
                                                                                                                                                             
    if (p[2] < bmin_[2])
      bmin_[2] = p[2];
    if (p[2] > bmax_[2])
      bmax_[2] = p[2];

    p = tripool[i].third();
                                                                                                                                                             
    if (p[0] < bmin_[0])
      bmin_[0] = p[0];
    if (p[0] > bmax_[0])
      bmax_[0] = p[0];
                                                                                                                                                             
    if (p[1] < bmin_[1])
      bmin_[1] = p[1];
    if (p[1] > bmax_[1])
      bmax_[1] = p[1];
                                                                                                                                                             
    if (p[2] < bmin_[2])
      bmin_[2] = p[2];
    if (p[2] > bmax_[2])
      bmax_[2] = p[2];
  }

  updateData();
}

void FieldBoundingBox::regularize()
{
  Vec3d center_ = 0.5 * (bmin_ + bmax_);
  Vec3d halfside_ = 0.5* (bmax_ - bmin_);
                                                                                                                                                             
  double maxHalf = halfside_[0];
  if (halfside_[1] > maxHalf)
    maxHalf = halfside_[1];
  if (halfside_[2] > maxHalf)
    maxHalf = halfside_[2];
                                                                                                                                                             
  Vec3d cubeHalf_ = Vec3d(maxHalf,maxHalf,maxHalf);

  bmin_ = center_ - cubeHalf_;
  bmax_ = center_ + cubeHalf_;

  updateData();

}

void FieldBoundingBox::updateData()
{
  center_ = 0.5 * (bmin_ + bmax_);
  halfSides_ = 0.5 * (bmax_ - bmin_);
}

void FieldBoundingBox::verifyBox()
{
  if (!((bmin_[0] <= bmax_[0]) && (bmin_[1] <= bmax_[1]) && (bmin_[2] <= bmax_[2])))
    printf("Error: inconsistent bounding box.\n");
}

// should this be turned into a self-modifying function?
void FieldBoundingBox::expand(double expansionFactor)
{

  bmin_ = center_ - expansionFactor * halfSides_;
  bmax_ = center_ + expansionFactor * halfSides_;

  updateData();
}

bool FieldBoundingBox::lineSegmentIntersection(Vec3d segmentStart, Vec3d segmentEnd, Vec3d * intersection)
{
  /*
    Jernej Barbic, CMU, 2005
    adapted from:

    Fast Ray-Box Intersection
    by Andrew Woo
    from "Graphics Gems", Academic Press, 1990
 */
                                                                                                                                                             
  #define NUMDIM  3
  #define RIGHT   0
  #define LEFT    1
  #define MIDDLE  2

  Vec3d dir = segmentEnd - segmentStart;

  bool inside = true;
  char quadrant[NUMDIM];
  register int i;
  int whichPlane;
  double maxT[NUMDIM];
  double candidatePlane[NUMDIM];

  /* Find candidate planes; this loop can be avoided if
  rays cast all from the eye(assume perpsective view) */
  for (i=0; i<NUMDIM; i++)
          if(segmentStart[i] < bmin_[i]) {
                  quadrant[i] = LEFT;
                  candidatePlane[i] = bmin_[i];
                  inside = false;
          }else if (segmentStart[i] > bmax_[i]) {
                  quadrant[i] = RIGHT;
                  candidatePlane[i] = bmax_[i];
                  inside = false;
          }else   {
                  quadrant[i] = MIDDLE;
          }
                                                                                                                                                             
  /* Ray origin inside bounding box */
  if(inside)      {
          *intersection = segmentStart;
          return (true);
  }
                                                                                                                                                             
                                                                                                                                                             
  /* Calculate T distances to candidate planes */
  for (i = 0; i < NUMDIM; i++)
          if (quadrant[i] != MIDDLE && dir[i] !=0.)
                  maxT[i] = (candidatePlane[i]-segmentStart[i]) / dir[i];
          else
                  maxT[i] = -1.;
                                                                                                                                                             
  /* Get largest of the maxT's for final choice of intersection */
  whichPlane = 0;
  for (i = 1; i < NUMDIM; i++)
          if (maxT[whichPlane] < maxT[i])
                  whichPlane = i;
                                                                                                                                                             
  /* Check final candidate actually inside box */
  if (maxT[whichPlane] < 0.) return (false);
  if (maxT[whichPlane] > 1.) return (false); // remove this for ray

  for (i = 0; i < NUMDIM; i++)
          if (whichPlane != i) {
                  (*intersection)[i] = segmentStart[i] + maxT[whichPlane] *dir[i];
                  if ((*intersection)[i] < bmin_[i] || (*intersection)[i] > bmax_[i])
                          return (false);
          } else {
                  (*intersection)[i] = candidatePlane[i];
          }
  return (true);                          /* ray hits box */
}

void FieldBoundingBox::print()
{
  cout << bmin_ << " " << bmax_ << endl;
}

template FieldBoundingBox::FieldBoundingBox(vector<TriangleBasic> & tripool);
template FieldBoundingBox::FieldBoundingBox(vector<TriangleWithCollisionInfo> & tripool);
template FieldBoundingBox::FieldBoundingBox(vector<TriangleWithCollisionInfoAndPseudoNormals> & tripool);
