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
    std::cout << modeData << std::endl;

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
//##############################################################################
void ModalAnalysisObject::
GetVolumeVertexDisplacement(const Eigen::VectorXd &q, Eigen::VectorXd &u)
{
    if (q.size() != N_Modes()) 
        throw std::runtime_error("**ERROR** Input reduced q has wrong dimension: "+std::to_string(q.size())+"->"+std::to_string(N_Modes())); 
    u = _eigenVectors * q; 
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
std::ostream &operator <<(std::ostream &os, const ModalAnalysisObject &object) 
{
    os << "--------------------------------------------------------------------------------\n" 
       << "Class ModalAnalysisObject\n" 
       << "--------------------------------------------------------------------------------\n"
       << " number of modes: " << object.N_Modes() << "\n"
       << " number of DOF  : " << object.N_DOF() << "\n"
       << "--------------------------------------------------------------------------------" 
       << std::flush; 
    return os; 
}
