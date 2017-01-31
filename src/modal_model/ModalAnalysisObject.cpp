#include <modal_model/ModalAnalysisObject.h> 

//##############################################################################
//##############################################################################
void ModalAnalysisObject::
ReadModeFromFile()
{
    if (!_modeFileSet) {
        std::cerr << "**WARNING** Mode file not set yet. Returning." << std::endl;
        return; }
    // first read the modes into ModeData struct
    ModeData modeData; 
    modeData.read(_modeFile.c_str()); 

    // fit data into Eigen matrix for later fast query
    const int N_modes = modeData.numModes(); 
    const int N_DOF   = modeData.numDOF(); 
    _eigenValues.resize(N_modes); 
    _eigenVectors.resize(N_DOF, N_modes); 
   
    // copy the matrices
    for (int m_idx=0; m_idx<N_modes; ++m_idx) 
    {
        _eigenValues(m_idx) = modeData.omegaSquared(m_idx); 
        std::vector<REAL> &mode = modeData.mode(m_idx); 
        for (int d_idx=0; d_idx<N_DOF; ++d_idx) 
            _eigenVectors(d_idx, m_idx) = mode.at(d_idx); 
    }
}

////##############################################################################
////##############################################################################
//void ModalAnalysisObject::
//GetForceInModalSpace(const int &vertexID, const Eigen::Vector3d &impulse, Eigen::VectorXd &forceInModalSpace)
//{
//    if (vertexID >= N_vertices())
//        throw std::runtime_error("**ERROR** vertexID"+std::to_string(vertexID)+" out of bounds. Total #vertices = "+std::to_string(N_vertices())); 
//    const int startIndex = vertexID*3; 
//    forceInModalSpace = _eigenVectors.row(startIndex+0)*impulse(0) 
//                      + _eigenVectors.row(startIndex+1)*impulse(1) 
//                      + _eigenVectors.row(startIndex+2)*impulse(2); 
//}

//##############################################################################
// This function computes u = U q. Depending on whether culling has been
// performed, the dimension of the returned u can be 3N_V or N_S. 
//##############################################################################
void ModalAnalysisObject::
GetVolumeVertexDisplacement(const Eigen::VectorXd &q, Eigen::VectorXd &u)
{
    if (q.size() != N_Modes()) 
        throw std::runtime_error("**ERROR** Input reduced q has wrong dimension: "+std::to_string(q.size())+"->"+std::to_string(N_Modes())); 

    if (_idType == VOL_UNCULL)
        u = _eigenVectors * q; 
    else 
        u = _eigenVectorsNormal * q; // used the projected if available. 
}

//##############################################################################
//##############################################################################
void ModalAnalysisObject::
GetVolumeVertexModeValues(const int &modeIndex, Eigen::VectorXd &modeValues)
{
    if (modeIndex >= N_Modes()) 
        throw std::runtime_error("**ERROR** mode index "+std::to_string(modeIndex)+" out of bounds"); 
    const int N = N_vertices(); 
    modeValues.resize(N); 
    for (int v_idx=0; v_idx<N; ++v_idx) 
    {
        const int d_idx = v_idx*3; 
        modeValues(v_idx) = pow(_eigenVectors(d_idx+0, modeIndex), 2)
                          + pow(_eigenVectors(d_idx+1, modeIndex), 2)
                          + pow(_eigenVectors(d_idx+2, modeIndex), 2); 
        modeValues(v_idx) = sqrt(modeValues(v_idx)); 
    }
}

//##############################################################################
//##############################################################################
void ModalAnalysisObject::
SetODESolverTime(const REAL &time)
{
    _time = time; 
    for (auto &odeSolver : _modalODESolvers)
        odeSolver->SetODECurrentTime(time);
}

//##############################################################################
//##############################################################################
void ModalAnalysisObject::
Initialize(const REAL &ODEStepSize, const std::string &modeFile, std::shared_ptr<ModalMaterial> materialPtr)
{
    _material = materialPtr; 
    _ODEStepSize = ODEStepSize; 
    SetModeFile(modeFile); 
    ReadModeFromFile(); 
    InitializeModalODESolvers(materialPtr); 
}

//##############################################################################
//##############################################################################
void ModalAnalysisObject::
InitializeModalODESolvers(std::shared_ptr<ModalMaterial> materialPtr)
{
    _modalODESolvers.resize(N_Modes()); 
    for (int mode_idx=0; mode_idx<N_Modes(); ++mode_idx) 
    {
        const REAL omegaSquared = _eigenValues(mode_idx) * materialPtr->inverseDensity;  // scaled it back, see Eq. 9 DyRT
        _modalODESolvers.at(mode_idx) = std::make_shared<ModalODESolver>(); 
        _modalODESolvers.at(mode_idx)->Initialize(materialPtr, omegaSquared, _ODEStepSize);  
    }
}

//##############################################################################
//##############################################################################
void ModalAnalysisObject:: 
InitializeSparseModalEncoder()
{
    // NOTE always uses projected modal matrix
    _modalAccEncoder = std::make_shared<SparseModalEncoder>(_eigenVectorsNormal);
}

//##############################################################################
// This function performs two main operations: 
//  1) remove nonsurface entries and reindex the eigen vectors. 
//  2) premultiply each vertex with the normal such that the number of rows
//     for eigenvectors are 3N->N, where N is number of surface vertices
//
//  At the end of the function, eigenvectors have number of rows reduced down to N_S
//  from 3N_V, where N_S is number of surface vertices and N_V is number  of
//  volume vertices.
//##############################################################################
void ModalAnalysisObject::
CullNonSurfaceModeShapes(std::shared_ptr<TetMeshIndexToSurfaceMesh> idMapPtr, std::shared_ptr<TriangleMesh<REAL> > meshPtr)
{
    // after modes are read, this function should only be called once
    if (_idType == SURF_CULLED) 
        return; 

    const int N_surfaceVertices = idMapPtr->N_surfaceVertices(); 
    const int N_volumeVertices  = N_vertices(); 
    const int N_modes = _eigenVectors.cols(); 
    Eigen::MatrixXd culledEigenVectors(N_surfaceVertices*3, N_modes); 
    _eigenVectorsNormal.resize(N_surfaceVertices, N_modes); 
    for (int vol_idx=0; vol_idx<N_volumeVertices; ++vol_idx)
    {
        if (!idMapPtr->KeyExists(vol_idx)) // skip vol indices that are interior
            continue; 
        const int surf_idx = idMapPtr->GetSurfaceIndex(vol_idx);
        const Vector3d &normal = meshPtr->normals().at(surf_idx); 
        // copy the block of eigenvectors and project to normal direction
        culledEigenVectors.row(surf_idx*3 + 0) = _eigenVectors.row(vol_idx*3 + 0); 
        culledEigenVectors.row(surf_idx*3 + 1) = _eigenVectors.row(vol_idx*3 + 1); 
        culledEigenVectors.row(surf_idx*3 + 2) = _eigenVectors.row(vol_idx*3 + 2); 
        _eigenVectorsNormal.row(surf_idx) = _eigenVectors.row(vol_idx*3 + 0) * normal.x
                                          + _eigenVectors.row(vol_idx*3 + 1) * normal.y
                                          + _eigenVectors.row(vol_idx*3 + 2) * normal.z;
    }

    // apply changes to member fields
    _eigenVectors = culledEigenVectors; 
    _tetMeshIndexToSurfaceMesh = idMapPtr; 
    _idType = SURF_CULLED; 
}

//##############################################################################
//##############################################################################
void ModalAnalysisObject::
CullHighFrequencyModes(const int &modesToKeep)
{
    if (modesToKeep >= N_Modes()) 
        return;

    // conservative resizing eigenvalues etc
    _eigenValues.conservativeResize(modesToKeep);
    _eigenVectors.conservativeResize(Eigen::NoChange, modesToKeep); 
    if (IDTypeIsSurf())
        _eigenVectorsNormal.conservativeResize(Eigen::NoChange, modesToKeep); 

    // reinitialize modal ode if material has been defined.
    if (_material)
    {
        _modalODESolvers.clear(); 
        InitializeModalODESolvers(_material);
    }
}

//##############################################################################
//##############################################################################
void ModalAnalysisObject::
CullHighFrequencyModes(const REAL &frequenciesToKeep)
{
    const int N_modes = N_Modes(); 
    int cullModeStart = -1;
    for (int mode_idx=0; mode_idx<N_modes; ++mode_idx)
    {
        if (GetModeFrequency(mode_idx) > frequenciesToKeep)
        {
            cullModeStart = mode_idx; 
            break;
        }
    }
    if (cullModeStart >= 0)
        CullHighFrequencyModes(cullModeStart);
    std::cout << *this << std::endl;
}

//##############################################################################
//##############################################################################
std::ostream &operator <<(std::ostream &os, const ModalAnalysisObject &object) 
{
    os << "--------------------------------------------------------------------------------\n" 
       << "Class ModalAnalysisObject\n" 
       << "--------------------------------------------------------------------------------\n"
       << " number of modes   : " << object.N_Modes() << "\n"
       << " number of DOF     : " << object.N_DOF() << "\n"
       << " ........................\n";
    if (object.N_Modes() > 0)
    {
    os << " highest frequency : " << object.GetModeFrequency(object.N_Modes()-1) << "\n" 
       << " lowest frequency  : " << object.GetModeFrequency(0) << "\n";
    }
    os << "--------------------------------------------------------------------------------" 
       << std::flush; 
    return os; 
}
