#include <fstream>
#include <config.h>
#include <geometry/TetMeshIndexToSurfaceMesh.h> 

//##############################################################################
//##############################################################################
void TetMeshIndexToSurfaceMesh::
ReadFromGeoFile(const std::string &geoFile)
{
    std::ifstream inFile(geoFile); 
    if (!inFile)
        throw std::runtime_error("**ERROR** *.geo file cannot be opened: "+geoFile); 

    // read number of vertices
    int N_vertices = 0; 
    inFile >> N_vertices; 

    // read all indices
    int tetIndex, surfIndex; 
    REAL nx, ny, nz, area; 
    for (int idx=0; idx<N_vertices; ++idx) 
    {
        inFile >> tetIndex >> surfIndex >> nx >> ny >> nz >> area; 
        if (_indexMap.find(tetIndex) == _indexMap.end())
            _indexMap[tetIndex] = surfIndex; 
        else 
            throw std::runtime_error("**ERROR** Two identical tet mesh index found in the geofile: "+geoFile); 
    }
}
