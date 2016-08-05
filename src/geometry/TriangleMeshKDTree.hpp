#ifndef TRIANGLE_MESH_KD_TREE_H
#define TRIANGLE_MESH_KD_TREE_H 
#include <kdtree.h>
#include <geometry/TriangleMesh.hpp> 

//##############################################################################
// This class will put the triangle mesh class into a kd-tree and support
// subsequent query about nearest neighbours. 
//
// Currently, the kd-tree building and querying is done via external library
// VLFeat (http://www.vlfeat.org/index.html). Documentation for the kd-tree
// operations can be found here: http://www.vlfeat.org/api/kdtree.html 
//##############################################################################
template <typename T> 
class TriangleMeshKDTree : public TriangleMesh<T>
{
};

#endif
