//==============================================================================
//  Sphere.h
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

#include "boundingBox.h"
#include "sphere.h"

// does the sphere intersect the bounding box
bool Sphere::doesBoundingBoxIntersect(FieldBoundingBox & box) const
{
  Vec3d bmin,bmax;
  bmin = box.bmin();
  bmax = box.bmax();

  double d;

  #define COORDINATE_TEST(i)\
  d = bmin[i] - center[i];\
  if (d > 0)\
    dmin += d * d;\
  d = center[i] - bmax[i];\
  if (d > 0)\
    dmin += d * d;

  double dmin = 0;
  COORDINATE_TEST(0)
  COORDINATE_TEST(1)
  COORDINATE_TEST(2)

  return (dmin <= radius * radius);
}
    

