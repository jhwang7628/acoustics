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

//##############################################################################
//##############################################################################
void FDTD_RigidSoundObject::
GetForceInModalSpace(const ImpactRecord &record, Eigen::VectorXd &forceInModalSpace)
{
    // surface mesh vertex id
    const int &surfaceVertexID = record.appliedVertex; 
    const Vector3d &force = record.impactVector; 
    
    // obtain tet mesh vertex id from the mapping
    const int vertexID = _tetMeshIndexToSurfaceMesh->GetTetIndex(surfaceVertexID); 

    if (vertexID >= ModalAnalysisObject::N_vertices())
        throw std::runtime_error("**ERROR** vertexID"+std::to_string(vertexID)+" out of bounds. Total #vertices = "+std::to_string(ModalAnalysisObject::N_vertices())); 
    const int startIndex = vertexID*3; 
    forceInModalSpace.resize(_eigenVectors.cols());
    forceInModalSpace = _eigenVectors.row(startIndex+0)*force.x 
                      + _eigenVectors.row(startIndex+1)*force.y 
                      + _eigenVectors.row(startIndex+2)*force.z; 
}

//##############################################################################
// The current implementation might be slow, but this should be used only for
// debugging.
//##############################################################################
void FDTD_RigidSoundObject:: 
GetModalDisplacement(const int &mode, Eigen::VectorXd &displacement)
{
    // first get the entire vector defined on volumetric mesh
    Eigen::VectorXd tetDisplacement; 
    Eigen::VectorXd qMode(_qNew.size()); 
    qMode.setZero(); 
    qMode(mode) = _qNew(mode); 
    GetVolumeVertexDisplacement(qMode, tetDisplacement); 

    // reorder the entries such that we return the mode values on the surface
    // mesh indices 
    const int N_volumeIndices = tetDisplacement.size()/3; 
    const int N_surfaceIndices = _mesh->num_vertices(); 
    const auto &map = _tetMeshIndexToSurfaceMesh; 

    _activeModeValues.resize(N_surfaceIndices); 
    for (int vol_idx=0; vol_idx<N_volumeIndices; ++vol_idx)
    {
        if (!map->KeyExists(vol_idx))
            continue; 

        // key exists, meaning there is a corresponding surface mesh idx
        const int surf_idx = map->GetSurfaceIndex(vol_idx); 
        const int flatIndex = vol_idx*3; 
        const Vector3d &normals = _mesh->normals().at(surf_idx); 
        const REAL modeValue = tetDisplacement(flatIndex+0) * normals.x
                             + tetDisplacement(flatIndex+1) * normals.y
                             + tetDisplacement(flatIndex+2) * normals.z;

        _activeModeValues(surf_idx) = modeValue; 
    }

    _activeModeAttributes.modeMax = _activeModeValues.maxCoeff(); 
    _activeModeAttributes.modeMin = _activeModeValues.minCoeff(); 
    _activeModeAttributes.absMax  = std::max<REAL>(fabs(_activeModeAttributes.modeMax), fabs(_activeModeAttributes.modeMin)); 
    displacement = _activeModeValues / 1e-5; 
}

//##############################################################################
// Advance the modal ODEs for N_steps
//##############################################################################
void FDTD_RigidSoundObject::
AdvanceModalODESolvers(const int &N_steps)
{
    // retrieve impact records (forces) within the time range 
    std::vector<ImpactRecord> impactRecords; 
    Eigen::VectorXd forceTimestep, forceBuffer; 
    for (int ts_idx=0; ts_idx<N_steps; ++ts_idx)
    {
        const REAL tsTimeStart = _time; 
        const REAL tsTimeStop  = _time + _ODEStepSize; 
        GetForces(tsTimeStart, tsTimeStop, impactRecords); 
        forceTimestep.setZero(N_Modes()); 
        const int N_records = impactRecords.size(); 
        // for each impact record within this timestep, concatenate the forces
        for (int rec_idx=0; rec_idx<N_records; ++rec_idx) 
        {
            GetForceInModalSpace(impactRecords.at(rec_idx), forceBuffer); 
            forceTimestep += forceBuffer; 
        }
        // step the system using force computed
        for (int mode_idx=0; mode_idx<N_Modes(); ++mode_idx) 
            _modalODESolvers.at(mode_idx)->StepSystem(_qOld(mode_idx), _qNew(mode_idx), forceTimestep(mode_idx)); 
        _time += _ODEStepSize;
        std::cout << "time = " << _time << std::endl; 
    }
}
