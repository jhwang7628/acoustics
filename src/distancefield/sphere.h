#ifndef __SPHERE_H__
#define __SPHERE_H__

#ifdef WIN32
  #include <windows.h>
#endif


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
#include "minivector.h"

class Sphere
{
public:

  Sphere(Vec3d center_g, double radius_g) : center(center_g), radius(radius_g) {}

  Sphere(double x, double y, double z, double radius_g) : center(Vec3d(x,y,z)), radius(radius_g) {}

  // does the sphere intersect the bounding box
  bool doesBoundingBoxIntersect(FieldBoundingBox & box) const;
    

private:

  Vec3d center;
  double radius;
};

#endif
