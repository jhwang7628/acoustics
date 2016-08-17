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

    _q_p.setZero(N_modes); 
    _q_c.setZero(N_modes); 
    _q_n.setZero(N_modes); 
    _q_nn.setZero(N_modes); 

    _qDot_c_minus.setZero(N_modes); 
    _qDDot_c.setZero(N_modes); 
    _qDDot_c_plus.setZero(N_modes); 
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
        Eigen::VectorXd qMode(_q_c.size()); 
        qMode.setZero(); 
        qMode(mode) = _q_c(mode); 
        GetVolumeVertexDisplacement(qMode, displacement); 
    }
    else 
    {
        GetVolumeVertexDisplacement(_q_c, displacement); 
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
        const REAL h2 = pow(_ODEStepSize, 2);
        _timer_substep_advanceODE[2].Start(); 
        _q_nn = _coeff_qNew.cwiseProduct(_q_n) + _coeff_qOld.cwiseProduct(_q_c) + _coeff_Q.cwiseProduct(forceTimestep);
        _qDot_c_minus = (_q_c - _q_p) / _ODEStepSize;
        _qDDot_c = (_q_n + _q_p - 2.0*_q_c) / h2;
        _qDDot_c_plus = 0.5 * ((_q_nn + _q_c - 2.0*_q_n)/h2 + _qDDot_c);
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
        of_q.write((char*)_q_c.data(), sizeof(double)*_q_c.size()); 
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
// This function should be called after acoustic time stepping has been
// completed.
//##############################################################################
void FDTD_RigidSoundObject::
UpdateQPointers()
{
    if (IsModalObject())
    {
        _q_p = _q_c; 
        _q_c = _q_n; 
        _q_n = _q_nn; 

        _time += _ODEStepSize;
    }
}

//##############################################################################
// Brute force looping for now
//##############################################################################
REAL FDTD_RigidSoundObject::
SampleModalDisplacement(const Vector3d &samplePoint, const Vector3d &sampleNormal, const REAL &sampleTime) // TODO!
{
    std::cerr << "**WARNING** Modal displacement sampling is not implemented.\n";
    return 0.0;
}

//##############################################################################
// Brute force looping for now
//##############################################################################
REAL FDTD_RigidSoundObject::
SampleModalVelocity(const Vector3d &samplePoint, const Vector3d &sampleNormal, const REAL &sampleTime) // TODO!
{
    int closestIndex = -1;
    REAL closestDistance = std::numeric_limits<REAL>::max(); 
    // transform the sample point to object frame
    const Eigen::Vector3d samplePointObject_e = _modelingTransformInverse * Eigen::Vector3d(samplePoint.x, samplePoint.y, samplePoint.z); 
    const Vector3d samplePointObject(samplePointObject_e[0], samplePointObject_e[1], samplePointObject_e[2]);  

    // nearest neighbour lookup
    const std::vector<Point3<REAL> > &vertices = _mesh->vertices(); 
    const int N_vertices = vertices.size(); 
    for (int vert_idx=0; vert_idx<N_vertices; ++vert_idx)
    {
        const REAL distance = (vertices.at(vert_idx) - samplePointObject).normSqr();
        if (distance < closestDistance)
        {
            closestIndex = vert_idx; 
            closestDistance = distance; 
        }
    }

    // evaluate sample value
    REAL sampledValue; 
    if (EQUAL_FLOATS(sampleTime, _time - 1.5*_ODEStepSize)) // sample at current time
        sampledValue = _eigenVectorsNormal.row(closestIndex).dot(_qDot_c_minus); 
    else
        throw std::runtime_error("**ERROR** Queried timestamp unexpected for modal velocity sampling. Double check.");

    return sampledValue;
}

//##############################################################################
// Brute force looping for now
//##############################################################################
REAL FDTD_RigidSoundObject::
SampleModalAcceleration(const Vector3d &samplePoint, const Vector3d &sampleNormal, const REAL &sampleTime)
{
    int closestIndex = -1;
    REAL closestDistance = std::numeric_limits<REAL>::max(); 
    // transform the sample point to object frame
    const Eigen::Vector3d samplePointObject_e = _modelingTransformInverse * Eigen::Vector3d(samplePoint.x, samplePoint.y, samplePoint.z); 
    const Vector3d samplePointObject(samplePointObject_e[0], samplePointObject_e[1], samplePointObject_e[2]);  

    // nearest neighbour lookup
    const std::vector<Point3<REAL> > &vertices = _mesh->vertices(); 
    const int N_vertices = vertices.size(); 
    for (int vert_idx=0; vert_idx<N_vertices; ++vert_idx)
    {
        const REAL distance = (vertices.at(vert_idx) - samplePointObject).normSqr();
        if (distance < closestDistance)
        {
            closestIndex = vert_idx; 
            closestDistance = distance; 
        }
    }

    // evaluate sample values
    REAL sampledValue; 
    if (EQUAL_FLOATS(sampleTime, _time-_ODEStepSize)) // sample at current time
        sampledValue = _eigenVectorsNormal.row(closestIndex).dot(_qDDot_c); 
    else if (EQUAL_FLOATS(sampleTime, _time-0.5*_ODEStepSize))
        sampledValue = _eigenVectorsNormal.row(closestIndex).dot(_qDDot_c_plus); 
    else
        throw std::runtime_error("**ERROR** Queried timestamp unexpected for modal acceleration sampling. Double check.");

    return sampledValue;
}

//##############################################################################
// Brute force looping for now
//##############################################################################
REAL FDTD_RigidSoundObject::
PerfectHarmonics_SampleModalVelocity(const int &mode, const Vector3d &samplePoint, const Vector3d &sampleNormal, const REAL &sampleTime)
{
    int closestIndex = -1;
    REAL closestDistance = std::numeric_limits<REAL>::max(); 
    // transform the sample point to object frame
    const Eigen::Vector3d samplePointObject_e = _modelingTransformInverse * Eigen::Vector3d(samplePoint.x, samplePoint.y, samplePoint.z); 
    const Vector3d samplePointObject(samplePointObject_e[0], samplePointObject_e[1], samplePointObject_e[2]);  

    // nearest neighbour lookup
    const std::vector<Point3<REAL> > &vertices = _mesh->vertices(); 
    const int N_vertices = vertices.size(); 
    for (int vert_idx=0; vert_idx<N_vertices; ++vert_idx)
    {
        const REAL distance = (vertices.at(vert_idx) - samplePointObject).normSqr();
        if (distance < closestDistance)
        {
            closestIndex = vert_idx; 
            closestDistance = distance; 
        }
    }

    const REAL omega = 2.0 * M_PI * GetModeFrequency(mode); 
    REAL sampledValue; 
    if (EQUAL_FLOATS(sampleTime, _time - 1.5*_ODEStepSize)) // sample at current time
        sampledValue = _eigenVectorsNormal(closestIndex, mode) * (-omega * sin(omega * sampleTime)); 
    else
        throw std::runtime_error("**ERROR** Queried timestamp unexpected for modal velocity sampling. Double check.");

    return sampledValue;
}

//##############################################################################
// Brute force looping for now
//##############################################################################
REAL FDTD_RigidSoundObject::
PerfectHarmonics_SampleModalAcceleration(const int &mode, const Vector3d &samplePoint, const Vector3d &sampleNormal, const REAL &sampleTime)
{
    int closestIndex = -1;
    REAL closestDistance = std::numeric_limits<REAL>::max(); 
    // transform the sample point to object frame
    const Eigen::Vector3d samplePointObject_e = _modelingTransformInverse * Eigen::Vector3d(samplePoint.x, samplePoint.y, samplePoint.z); 
    const Vector3d samplePointObject(samplePointObject_e[0], samplePointObject_e[1], samplePointObject_e[2]);  

    // nearest neighbour lookup
    const std::vector<Point3<REAL> > &vertices = _mesh->vertices(); 
    const int N_vertices = vertices.size(); 
    for (int vert_idx=0; vert_idx<N_vertices; ++vert_idx)
    {
        const REAL distance = (vertices.at(vert_idx) - samplePointObject).normSqr();
        if (distance < closestDistance)
        {
            closestIndex = vert_idx; 
            closestDistance = distance; 
        }
    }

#ifdef DEBUG
#ifdef USE_OPENMP
#pragma omp critical
#endif
        {
            const Point3<REAL> &vertex = vertices.at(closestIndex); 
            const Eigen::Vector3d vertexWorld_e = _modelingTransform * Eigen::Vector3d(vertex.x, vertex.y, vertex.z); 
            const Vector3d vertexWorld(vertexWorld_e[0], vertexWorld_e[1], vertexWorld_e[2]); 
            // write these special points for debugging purpose 
            _debugArrowStart.push_back(vertexWorld); 
            _debugArrowNormal.push_back(-vertexWorld + samplePoint); 
        }
#endif

    const REAL omega = 2.0 * M_PI * GetModeFrequency(mode); 
    REAL sampledValue; 
    if (EQUAL_FLOATS(sampleTime, _time-_ODEStepSize))
        sampledValue = _eigenVectorsNormal(closestIndex, mode) * (-omega*omega * cos(omega * sampleTime));
    else if (EQUAL_FLOATS(sampleTime, _time-0.5*_ODEStepSize))
        sampledValue = _eigenVectorsNormal(closestIndex, mode) * (-omega*omega * cos(omega * sampleTime)); 
    else
        throw std::runtime_error("**ERROR** Queried timestamp unexpected for modal acceleration sampling. Double check.");

    return sampledValue;
}

//##############################################################################
//##############################################################################
void FDTD_RigidSoundObject::
PrintAllVelocity(const std::string &filename, const int &mode) const
{
    const std::vector<Point3<REAL> > &vertices = _mesh->vertices(); 
    const int N_vertices = vertices.size(); 
    const REAL omega = 2.0 * M_PI * GetModeFrequency(mode); 
    std::ofstream of(filename.c_str()); 
    for (int v_idx=0; v_idx<N_vertices; ++v_idx)
    {
        std::cout << "vid=" << v_idx << ", velocity BC = " << _eigenVectorsNormal(v_idx, mode) * (-omega) << std::endl;
        of << v_idx << " " << _eigenVectorsNormal(v_idx, mode) * (-omega) << std::endl;
    }
    of.close();
}
