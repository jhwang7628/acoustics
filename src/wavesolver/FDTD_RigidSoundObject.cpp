#include <wavesolver/FDTD_RigidSoundObject.h>
#include <wavesolver/Wavesolver_ConstantsAndTypes.h>
#include <utils/SimpleTimer.h>
#include <geometry/Point3.hpp>
#include <macros.h>

//##############################################################################
//##############################################################################
void FDTD_RigidSoundObject::
Initialize() // TODO! should manage all parent class initialization
{
    assert(FDTD_RigidObject::_tetMeshIndexToSurfaceMesh && FDTD_RigidObject::_mesh); 
    CullNonSurfaceModeShapes(FDTD_RigidObject::_tetMeshIndexToSurfaceMesh, FDTD_RigidObject::_mesh); 
    CullHighFrequencyModes(MODE_SHAPE_CUTOFF_FREQ);

    // Copy the fixed coefficient for each ODE for vectorized IIR
    const int N_modes = N_Modes();  
    _coeff_qNew.resize(N_modes); 
    _coeff_qOld.resize(N_modes); 
    _coeff_Q.resize(N_modes); 
    for (int mode_idx=0; mode_idx<N_modes; ++mode_idx) 
    {
        ModalODESolverPtr &odeSolverPtr = _modalODESolvers.at(mode_idx); 
        _coeff_qNew(mode_idx) = odeSolverPtr->GetCoefficient_qNew(); 
        _coeff_qOld(mode_idx) = odeSolverPtr->GetCoefficient_qOld(); 
        _coeff_Q(mode_idx) = odeSolverPtr->GetCoefficient_Q(); 
    }
    InitializeModeVectors(); 
}

//##############################################################################
//##############################################################################
void FDTD_RigidSoundObject::
InitializeModeVectors()
{
    const int N_modes = N_Modes(); 
    _qOld.setZero(N_modes); 
    _qNew.setZero(N_modes); 
    _qOldDot.setZero(N_modes); 
    _qOldDDot.setZero(N_modes); 
}

//##############################################################################
//##############################################################################
void FDTD_RigidSoundObject::
GetVertexModeValuesNormalized(const int &modeIndex, Eigen::VectorXd &modeValues)
{
    // first get the entire vector defined on volumetric mesh
    Eigen::VectorXd allModeValues; 
    GetVolumeVertexModeValues(modeIndex, allModeValues); 

    // reorder the entries such that we return the mode values on the surface
    // mesh indices 
    const int N_volumeIndices = allModeValues.size(); 
    const int N_surfaceIndices = _mesh->num_vertices(); 
    const auto &map = FDTD_RigidObject::_tetMeshIndexToSurfaceMesh; 

    modeValues.resize(N_surfaceIndices); 
    for (int vol_idx=0; vol_idx<N_volumeIndices; ++vol_idx)
    {
        if (!map->KeyExists(vol_idx))
            continue; 

        // key exists, meaning there is a corresponding surface mesh idx
        const int surf_idx = map->GetSurfaceIndex(vol_idx); 
        modeValues(surf_idx) = allModeValues(vol_idx); 
    }

    const REAL modeMax = modeValues.maxCoeff(); 
    const REAL modeMin = modeValues.minCoeff(); 
    const REAL absMax  = std::max<REAL>(fabs(modeMax), fabs(modeMin)); 
    modeValues = modeValues / absMax; 
}

//##############################################################################
// This function computes Q = U^T f. Note that whether culling was performed or
// not, the transformation should use _eigenVectors, which has dimension
// 3N_V before culling, and 3N_S after culling.
//##############################################################################
void FDTD_RigidSoundObject::
GetForceInModalSpace(const ImpactRecord &record, Eigen::VectorXd &forceInModalSpace)
{
    // surface mesh vertex id
    const Vector3d &force = record.impactVector; 
    const int &vertexID = record.appliedVertex; // FIXME debug

    if (vertexID >= ModalAnalysisObject::N_vertices())
        throw std::runtime_error("**ERROR** vertexID "+std::to_string(vertexID)+" out of bounds. Total #vertices = "+std::to_string(ModalAnalysisObject::N_vertices())); 
    const int startIndex = vertexID*3; 
    forceInModalSpace.resize(_eigenVectors.cols());
    forceInModalSpace = _eigenVectors.row(startIndex+0)*force.x 
                      + _eigenVectors.row(startIndex+1)*force.y 
                      + _eigenVectors.row(startIndex+2)*force.z; 
}

//##############################################################################
// Get the normal displacement
//##############################################################################
void FDTD_RigidSoundObject::
GetModalDisplacementAux(const int &mode, Eigen::VectorXd &displacement)
{
    assert(IDTypeIsSurf()); 

    _timer_substep_q2u[0].Start(); 
    if (mode >= 0)
    {
        // this implementation copies the vector, should only use when
        // debugging
        Eigen::VectorXd qMode(_qNew.size()); 
        qMode.setZero(); 
        qMode(mode) = _qNew(mode); 
        GetVolumeVertexDisplacement(qMode, displacement); 
    }
    else 
    {
        GetVolumeVertexDisplacement(_qNew, displacement); 
    }
    _timer_substep_q2u[0].Pause(); 
    displacement /= 1E-5;
}

//##############################################################################
//##############################################################################
void FDTD_RigidSoundObject:: 
GetModalDisplacement(const int &mode, Eigen::VectorXd &displacement)
{
    GetModalDisplacementAux(mode, displacement);
}

//##############################################################################
//##############################################################################
void FDTD_RigidSoundObject:: 
GetModalDisplacement(Eigen::VectorXd &displacement)
{
    GetModalDisplacementAux(-1, displacement);
}

//##############################################################################
// Advance the modal ODEs for N_steps
//##############################################################################
void FDTD_RigidSoundObject::
AdvanceModalODESolvers(const int &N_steps)
{
    for (int ts_idx=0; ts_idx<N_steps; ++ts_idx)
    {
        // retrieve impact records (forces) within the time range 
        std::vector<ImpactRecord> impactRecords; 
        Eigen::VectorXd forceTimestep, forceBuffer; 
        _timer_substep_advanceODE[0].Start(); 
        const REAL tsTimeStart = _time; 
        const REAL tsTimeStop  = _time + _ODEStepSize; 
        GetForces(tsTimeStart, tsTimeStop, impactRecords); 
        _timer_substep_advanceODE[0].Pause(); 

        _timer_substep_advanceODE[1].Start(); 
        forceTimestep.setZero(N_Modes()); 
        const int N_records = impactRecords.size(); 
        // for each impact record within this timestep, concatenate the forces
        for (int rec_idx=0; rec_idx<N_records; ++rec_idx) 
        {
            GetForceInModalSpace(impactRecords.at(rec_idx), forceBuffer); 
            forceTimestep += forceBuffer; 
        }
        _timer_substep_advanceODE[1].Pause(); 

        // step the system using force computed
        _timer_substep_advanceODE[2].Start(); 
        Eigen::VectorXd q = _coeff_qNew.cwiseProduct(_qNew) + _coeff_qOld.cwiseProduct(_qOld) + _coeff_Q.cwiseProduct(forceTimestep); 

        // compute velocity and acceleration using finite-difference. Notice
        // that this corresponds to the state at qOld after updates
        _qOldDot  = (q - _qOld) / (2.0*_ODEStepSize); 
        _qOldDDot = (q + _qOld - 2.0*_qNew) / pow(_ODEStepSize,2); 

        _qOld = _qNew; 
        _qNew = q; 

        _time += _ODEStepSize;
        std::cout << "modal ode time = " << _time << std::endl; 
        _timer_substep_advanceODE[2].Pause();
    }
}

//##############################################################################
// Advance the modal ODEs for N_steps with logging file
//##############################################################################
void FDTD_RigidSoundObject::
AdvanceModalODESolvers(const int &N_steps, std::ofstream &of_displacement, std::ofstream &of_q)
{
    {
        const int N_rows = N_steps; 
        const int N_cols = _mesh->num_vertices(); 
        of_displacement.write((char*)&N_rows, sizeof(int)); 
        of_displacement.write((char*)&N_cols, sizeof(int)); 
    }
    {
        const int N_rows = N_steps; 
        const int N_cols = N_Modes(); 
        of_q.write((char*)&N_rows, sizeof(int)); 
        of_q.write((char*)&N_cols, sizeof(int)); 
    }
    for (int ts_idx=0; ts_idx<N_steps; ++ts_idx)
    {
        Eigen::VectorXd displacements; 
        _timer_mainstep[0].Start(); 
        AdvanceModalODESolvers(1);
        _timer_mainstep[0].Pause(); 
        _timer_mainstep[1].Start(); 
        GetModalDisplacement(displacements);
        _timer_mainstep[1].Pause(); 
        _timer_mainstep[2].Start(); 
        of_displacement.write((char*)displacements.data(), sizeof(double)*displacements.size()); 
        of_q.write((char*)_qNew.data(), sizeof(double)*_qNew.size()); 
        _timer_mainstep[2].Pause(); 
    }
    std::cout << "Total Timing: \n"
              << " Advance ODEs      : " << _timer_mainstep[0].Duration()           << " sec\n"
              << "  Interpolate force: " << _timer_substep_advanceODE[0].Duration() << " sec\n"
              << "  Transform force  : " << _timer_substep_advanceODE[1].Duration() << " sec\n"
              << "  Step ODE         : " << _timer_substep_advanceODE[2].Duration() << " sec\n"
              << " q->u              : " << _timer_mainstep[1].Duration()           << " sec\n"
              << "  Get volume u     : " << _timer_substep_q2u[0].Duration()        << " sec\n"
              << "  Map to surface   : " << _timer_substep_q2u[1].Duration()        << " sec\n"
              << " Write data to disk: " << _timer_mainstep[2].Duration()           << " sec\n"; 
    const REAL scale = 1.0 / ((REAL)N_steps/1000.); 
    std::cout << "Average Timing: \n"
              << " Advance ODEs      : " << _timer_mainstep[0].Duration()           * scale << " ms\n"
              << "  Interpolate force: " << _timer_substep_advanceODE[0].Duration() * scale << " ms\n"
              << "  Transform force  : " << _timer_substep_advanceODE[1].Duration() * scale << " ms\n"
              << "  Step ODE         : " << _timer_substep_advanceODE[2].Duration() * scale << " ms\n"
              << " q->u              : " << _timer_mainstep[1].Duration()           * scale << " ms\n"
              << "  Get volume u     : " << _timer_substep_q2u[0].Duration()        * scale << " ms\n"
              << "  Map to surface   : " << _timer_substep_q2u[1].Duration()        * scale << " ms\n"
              << " Write data to disk: " << _timer_mainstep[2].Duration()           * scale << " ms\n"; 
}

//##############################################################################
// Brute force looping for now
//##############################################################################
REAL FDTD_RigidSoundObject::
SampleModalDisplacement(const Vector3d &samplePoint, const Vector3d &sampleNormal, const REAL &sampleTime) // TODO!
{
    std::cerr << "**WARNING** return zero for modal displacement since time does not match. interpolation not implemented.\n";
    return 0.0;
}

//##############################################################################
// Brute force looping for now
//##############################################################################
REAL FDTD_RigidSoundObject::
SampleModalVelocity(const Vector3d &samplePoint, const Vector3d &sampleNormal, const REAL &sampleTime) // TODO!
{
    std::cerr << "**WARNING** return zero for modal velocity since time does not match. interpolation not implemented.\n";
    return 0.0;
}

//##############################################################################
// Brute force looping for now
//##############################################################################
REAL FDTD_RigidSoundObject::
SampleModalAcceleration(const Vector3d &samplePoint, const Vector3d &sampleNormal, const REAL &sampleTime)
{
    int closestIndex = -1;
    REAL closestDistance = std::numeric_limits<REAL>::max(); 
    if (EQUAL_FLOATS(sampleTime, _time-_ODEStepSize))
    {
        const std::vector<Point3<REAL> > &vertices = _mesh->vertices(); 
        const int N_vertices = vertices.size(); 
        for (int vert_idx=0; vert_idx<N_vertices; ++vert_idx)
        {
            const REAL distance = (vertices.at(vert_idx) - samplePoint).normSqr();
            if (distance < closestDistance)
            {
                closestIndex = vert_idx; 
                closestDistance = distance; 
            }
        }
        const REAL sampledValue = _eigenVectorsNormal.row(closestIndex).dot(_qOldDDot); 
        return sampledValue; 
    }
    std::cerr << "**WARNING** return zero for modal acceleration since time does not match. interpolation not implemented.\n";
    return 0.0;
}

