#include <wavesolver/FDTD_RigidSoundObject.h>

//##############################################################################
//##############################################################################
void FDTD_RigidSoundObject::
SetVertexModeValues(const int &modeIndex)
{
    // first get the entire vector defined on volumetric mesh
    Eigen::VectorXd allModeValues; 
    GetVolumeVertexModeValues(modeIndex, allModeValues); 

    // reorder the entries such that we return the mode values on the surface
    // mesh indices 
    const int N_volumeIndices = allModeValues.size(); 
    const int N_surfaceIndices = _mesh->num_vertices(); 
    const auto &map = _tetMeshIndexToSurfaceMesh; 

    _activeModeValues.resize(N_surfaceIndices); 
    for (int vol_idx=0; vol_idx<N_volumeIndices; ++vol_idx)
    {
        if (!map->KeyExists(vol_idx))
            continue; 

        // key exists, meaning there is a corresponding surface mesh idx
        const int surf_idx = map->GetSurfaceIndex(vol_idx); 
        _activeModeValues(surf_idx) = allModeValues(vol_idx); 
    }

    _activeModeAttributes.modeMax = _activeModeValues.maxCoeff(); 
    _activeModeAttributes.modeMin = _activeModeValues.minCoeff(); 
    _activeModeAttributes.absMax  = std::max<REAL>(fabs(_activeModeAttributes.modeMax), fabs(_activeModeAttributes.modeMin)); 
}

//##############################################################################
//##############################################################################
void FDTD_RigidSoundObject::
GetVertexModeValues(const int &modeIndex, Eigen::VectorXd &modeValues)
{
    modeValues = _activeModeValues; 
}

//##############################################################################
//##############################################################################
void FDTD_RigidSoundObject::
GetVertexModeValuesNormalized(const int &modeIndex, Eigen::VectorXd &modeValues)
{
    modeValues = _activeModeValues / _activeModeAttributes.absMax; 
}
