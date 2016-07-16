#ifndef TET_MESH_INDEX_TO_SURFACE_MESH_H 
#define TET_MESH_INDEX_TO_SURFACE_MESH_H 
#include <map>
#include <string> 

//##############################################################################
// This class handles the index mapping between the tetrahedron mesh and
// surface mesh that is stored in the pipeline. 
//
// The tet mesh (*.tet) is produced by class isostuffer. It is also used by 
// elasticity solver and eigen solver, which means that the file *.modes is
// stored according to the vertices of the tet mesh. To render the modes
// correctly on the surface mesh, which is used throughout the FDTD simulator,
// a mapping needs to be constructed. This mapping is produced and stored by
// the class elasticity_solver and has extension *.geo.txt. It comes from an 
// std::map. The key is tet index and value is surface index. 
//
// The file format for *.geo.txt is prescribed in file
// tools/elasticity/elasticity_solver.cpp 
//
// <int: map_size>
// <int: tet_index> <int: surf_index> <float: surf_normal_x> <float: surf_normal_y> <float: surf_normal_z> <area> 
// <int: tet_index> <int: surf_index> <float: surf_normal_x> <float: surf_normal_y> <float: surf_normal_z> <area> 
// .
// .
// .
//
//##############################################################################
class TetMeshIndexToSurfaceMesh
{
    private: 
        std::map<int, int> _indexMap;

    public: 
        TetMeshIndexToSurfaceMesh(){}

        inline bool KeyExists(const int &tetIndex){return (_indexMap.find(tetIndex)!=_indexMap.end() ? true : false);}
        inline int N_surfaceVertices(){return _indexMap.size();}
        inline int GetSurfaceIndex(const int &tetIndex){return _indexMap.at(tetIndex);} 
        void ReadFromGeoFile(const std::string &geoFile); 
};

#endif
