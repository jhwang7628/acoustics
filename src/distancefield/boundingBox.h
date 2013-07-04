#ifndef __BOUNDINGBOX_H__
#define __BOUNDINGBOX_H__

//==============================================================================
//  FieldBoundingBox.h
//  Author: Jernej Barbic <barbic@cs.cmu.edu>
//
//==============================================================================

#ifdef WIN32
  #include <windows.h>
#endif

#include <vector>
#include <queue>
#include <string>
#include <map>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "minivector.h"

//class TriangleBasic;

class FieldBoundingBox
{
public:

  FieldBoundingBox(Vec3d bmin_g=Vec3d(0,0,0), Vec3d bmax_g=Vec3d(1,1,1)): bmin_(bmin_g), bmax_(bmax_g) { updateData();}
  template<class Triangle> FieldBoundingBox(std::vector<Triangle> & tripool);

  // accessors
  Vec3d bmin() const { return bmin_;}
  Vec3d bmax() const { return bmax_;}

  Vec3d center() { return center_;}
  Vec3d halfSides() { return halfSides_;}

  double diameter() { return len(bmax_-bmin_); }

  // mutators
  void setbmin(Vec3d bmin_g) { bmin_=bmin_g; updateData();}
  void setbmin(double x, double y, double z) { bmin_=Vec3d(x,y,z); updateData();}
  void setbmax(Vec3d bmax_g) { bmax_=bmax_g; updateData();}
  void setbmax(double x, double y, double z) { bmax_=Vec3d(x,y,z); updateData();}

  // sanity check bmin <= bmax
  void verifyBox();

  // expands from the center 
  // factor of 1.0 indicates no expansion
  void expand(double expansionFactor);
  void regularize(); // converts the box into one with all sides equal

  bool lineSegmentIntersection(Vec3d segmentStart, Vec3d segmentEnd, Vec3d * intersection);

  void print();

protected:

  void updateData(); // updates center and half-sides

  Vec3d center_, halfSides_;

  Vec3d bmin_,bmax_;
};

#endif
