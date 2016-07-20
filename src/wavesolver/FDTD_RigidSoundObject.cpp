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
    forceInModalSpace = _eigenVectors.row(startIndex+0)*force.x 
                      + _eigenVectors.row(startIndex+1)*force.y 
                      + _eigenVectors.row(startIndex+2)*force.z; 
}

//##############################################################################
// Advance the modal ODEs for N_steps
//##############################################################################
void FDTD_RigidSoundObject::
AdvanceModalODESolvers(const int &N_steps)
{
    // retrieve impact records (forces) within the time range 
    std::vector<ImpactRecord> impactRecords; 
    Eigen::VectorXd force; 
    //_q.resize(N_Modes());
    //Eigen::VectorXd _q(N_Modes()); 
    for (int ts_idx=0; ts_idx<N_steps; ++ts_idx)
    {
        const REAL tsTimeStart = _time; 
        const REAL tsTimeStop  = _time + _ODEStepSize; 
        GetForces(tsTimeStart, tsTimeStop, impactRecords); 

        std::cout << "------\n";
        std::cout << impactRecords.size() << "\n";

        const int N_records = impactRecords.size(); 
        for (int rec_idx=0; rec_idx<N_records; ++rec_idx) 
        {
            GetForceInModalSpace(impactRecords.at(rec_idx), force); 
            std::cout << "force=" << force.transpose() << std::endl;


            for (int mode_idx=0; mode_idx<N_Modes(); ++mode_idx) 
                _modalODESolvers.at(mode_idx)->StepSystem(_qOld(mode_idx), _qNew(mode_idx), force(mode_idx)); 
        }
        // DEBUG FIXME
        std::cout << ">>> " << tsTimeStart << ": " << _qNew.transpose() << std::endl;

        _time += _ODEStepSize;
    }

//    // get impulse record 
//    Eigen::VectorXd forceInModalSpace; // should have length N_modes
//    GetForceInModalSpace(vertexID, impulse, forceInModalSpace); 
//    for (int mode_idx=0; mode_idx<N_Modes(); ++mode_idx) 
//    {
//        _modalODESolvers.at(mode_idx)->StepSystem(forceInModalSpace(mode_idx)); 
//    }
}
