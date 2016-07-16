#include <wavesolver/FDTD_RigidSoundObject.h>

//##############################################################################
//##############################################################################
void FDTD_RigidSoundObject::
GetVertexModeValues(const int &modeIndex, Eigen::VectorXd &modeValues)
{
    // first get the entire vector defined on volumetric mesh
    Eigen::VectorXd allModeValues; 
    GetVolumeVertexModeValues(modeIndex, allModeValues); 

    // reorder the entries such that we return the mode values on the surface
    // mesh indices 
    const int N_volumeIndices = allModeValues.size(); 
    const int N_surfaceIndices = _mesh->num_vertices(); 
    const auto &map = _tetMeshIndexToSurfaceMesh; 

    modeValues.resize(N_surfaceIndices); 
    for (int vol_idx=0; vol_idx<N_volumeIndices; ++vol_idx)
    {
        if (!map->KeyExists(vol_idx))
            continue; 

        // key exists, meaning there is a corresponding surface mesh idx
        const int surf_idx = map->GetSurfaceIndex(vol_idx); 
        modeValues(surf_idx) = allModeValues(vol_idx); 
    }
}
