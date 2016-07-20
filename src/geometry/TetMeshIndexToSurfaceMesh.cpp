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
        // push the (tet, surf) into _indexMap
        if (_indexMap.find(tetIndex) == _indexMap.end())
            _indexMap[tetIndex] = surfIndex; 
        else 
            throw std::runtime_error("**ERROR** Two identical tet mesh index found in the geofile: "+geoFile); 

        // push the (surf, tet) into _indexMap
        if (_inverseIndexMap.find(surfIndex) == _inverseIndexMap.end())
            _inverseIndexMap[surfIndex] = tetIndex; 
        else 
            throw std::runtime_error("**ERROR** Two identical surface mesh index found in the geofile: "+geoFile); 
    }
}
