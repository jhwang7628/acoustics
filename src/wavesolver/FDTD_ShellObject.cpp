#include "wavesolver/FDTD_ShellObject.h"
//##############################################################################
// Function _BuildDistanceField
//##############################################################################
void FDTD_ShellObject::
Initialize()
{
    // find bounding box 
    Point3d  b_l, b_h, b_c; 
    _mesh->bounding_box(b_l, b_h); 
    b_c = (b_l + b_h)/2.0;
    const Vector3d h = (b_h - b_l)/2.0*1.05; // enlarge by 105%
    b_l = b_c - h; 
    b_h = b_c + h; 
    _bbox = BoundingBox(b_l, b_h); 
    _init = true; 
}
  
//##############################################################################
// Function DistanceToMesh
//##############################################################################
REAL FDTD_ShellObject:: 
DistanceToMesh(const Vector3d &position)
{
    if (!_init)
        Initialize();
    if (!_bbox.Inside(position))
        return std::numeric_limits<REAL>::max(); 
    const auto p_o = WorldToObjectPoint(position);
    Vector3d bp, pp; //boundary, projected
    int tri_idx; 
    REAL d; 
#pragma omp critical
    d = _mesh->ComputeClosestPointOnMesh(p_o, pp, tri_idx, pp);
    return d;
}
