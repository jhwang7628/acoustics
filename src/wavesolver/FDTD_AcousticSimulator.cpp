#include <wavesolver/FDTD_AcousticSimulator.h>
#include <wavesolver/SimWorldAuxDS.h>
#include <io/FDTD_ListenShell.hpp>
#include <geometry/BoundingBox.h>
#include <utils/IO/IO.h>
#include <macros.h>
//#define DEBUG_HARMONIC_SOURCE
// uncomment this if you want the Modal ODEs to jump directly to the start time
// without the effects of earlier impulses
//#define JUMP_MODAL_ODES_TO_START


// THESE TWO ARE DEPRECATED -- USE "has_modal_source" and "has_acc_noise_source"
// in the config xml file instead
//#define ACC_NOISE_SOURCE
//#define MODAL_SOURCE
//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
_ParseSolverSettings()
{
    _parser = std::make_shared<ImpulseResponseParser>(_configFile);
    _acousticSolverSettings = std::make_shared<PML_WaveSolver_Settings>();
    _sceneObjects = std::make_shared<FDTD_Objects>();

    _parser->GetSolverSettings(_acousticSolverSettings);
    _parser->GetObjects(_acousticSolverSettings, _sceneObjects);
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
_SetBoundaryConditions() // NOTE: decrecated, not called anymore
{
    auto &objects = _sceneObjects->GetRigidObjects();
    for (auto &m : objects)
    {
        RigidObjectPtr objectPtr = m.second;
#ifdef MODAL_SOURCE
        // add modal vibrational source
        VibrationalSourcePtr sourcePtr(new ModalVibrationalSource(objectPtr));
        objectPtr->AddVibrationalSource(sourcePtr);
#endif
#ifdef ACC_NOISE_SOURCE
        // add acceleration noise source
        VibrationalSourcePtr anSourcePtr(new AccelerationNoiseVibrationalSource(objectPtr));
        objectPtr->AddVibrationalSource(anSourcePtr);
#endif
#ifdef DEBUG_HARMONIC_SOURCE
        // add debug harmonic source
        const REAL omega = 2.0*M_PI*1500.0;
        //const REAL omega = 2.0*M_PI*objectPtr->GetModeFrequency(0);
        const REAL phase = 0.0;
        VibrationalSourcePtr hSourcePtr(new HarmonicVibrationalSource(objectPtr, omega, phase));
        //VibrationalSourcePtr hSourcePtr(new HarmonicVibrationalSource(objectPtr, omega, phase, 1000., 0));
        //std::dynamic_pointer_cast<HarmonicVibrationalSource>(hSourcePtr)->SetRange(0., 1./1500.);
        objectPtr->AddVibrationalSource(hSourcePtr);
#endif
        //objectPtr->TestObjectBoundaryCondition();
    }
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
_SetPressureSources()
{
    std::vector<PressureSourcePtr> &scenePressureSources = _sceneObjects->GetPressureSources();
   _parser->GetPressureSources(_acousticSolverSettings->soundSpeed, scenePressureSources);
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
_SetListeningPoints()
{
    Vector3Array &listeningPoints = (_acousticSolverSettings->useShell ?
            _owner->listen->speakers :
            _acousticSolverSettings->listeningPoints);
    if (_acousticSolverSettings->useShell && _owner)
    {
        ListeningUnit::refShell = FDTD_ListenShell<REAL>::Build(_acousticSolverSettings->refShellFile);
        const REAL r_out = _owner->GetBoundingBox().minlength()/2.0
                         - (_acousticSolverSettings->PML_width+1)*_acousticSolverSettings->cellSize;
        const REAL r_inn = r_out - _acousticSolverSettings->spacing;
        assert(r_out>0. && r_inn>0.);
        _owner->listen->outShell =
            std::make_shared<FDTD_ListenShell<REAL>>(
                    ListeningUnit::refShell.get(),
                    _owner->BoundingBoxCenter(),
                    r_out);
        _owner->listen->innShell =
            std::make_shared<FDTD_ListenShell<REAL>>(
                    ListeningUnit::refShell.get(),
                    _owner->BoundingBoxCenter(),
                    r_inn);
        _owner->UpdateSpeakers();
    }
    else
    {
        _parser->GetListeningPoints(listeningPoints);
    }

    _acousticSolverSettings->listeningFile = CompositeFilename("listening_pressure_%6u.dat");
}

//##############################################################################
//##############################################################################
std::string FDTD_AcousticSimulator::
CompositeFilename(const std::string filename_)
{
    char buffer[512];
    buffer[0] = 0;

    if(_acousticSolverSettings->timeParallel)
    {
        snprintf(buffer, 512, "%05d_", _acousticSolverSettings->indTimeChunks);
    }

    std::string filename =
        (_simulatorID ? (*_simulatorID) + "_" + buffer + filename_ : filename_);
    snprintf(buffer, 512, _acousticSolverSettings->outputPattern.c_str(),
             filename.c_str());
    return std::string(buffer);
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
_SaveSolverSettings(const std::string &filename)
{
    std::ofstream of(filename.c_str());
    of << *_acousticSolver << std::endl;
    of << _acousticSolver->GetGrid() << std::endl;
    _sceneObjects->PrintAllSources(of);

    std::string xmlFilename_s = CompositeFilename("input_control_file.xml");
    std::cout << "copy file from " << _configFile
              << " to " << xmlFilename_s << std::endl;
    IO::CopyFile(_configFile, xmlFilename_s);
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
_SavePressureCellPositions(const std::string &filename)
{
    const Tuple3i &pFieldDivisions = GetGrid().pressureFieldDivisions();
    const int N_cells = _acousticSolver->numCells();
    Eigen::MatrixXd vertexPosition(N_cells,3);
    int count = 0;
    Vector3d vPosition;
    for (int kk=0; kk<pFieldDivisions[2]; kk++)
        for (int jj=0; jj<pFieldDivisions[1]; jj++)
            for (int ii=0; ii<pFieldDivisions[0]; ii++)
            {
                Tuple3i  vIndex(ii, jj, kk);
                vPosition = _acousticSolver->fieldPosition(vIndex);
                vertexPosition(count, 0) = vPosition[0];
                vertexPosition(count, 1) = vPosition[1];
                vertexPosition(count, 2) = vPosition[2];
                count ++;
            }

    try
    {
        IO::writeMatrixX<double>(vertexPosition, filename.c_str(), IO::BINARY);
    }
    catch (const std::runtime_error &e)
    {
    }
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
_SaveVelocityCellPositions(const std::string &filename, const int &dim)
{
    const Tuple3i &divisions = _acousticSolver->velocityFieldDivisions(dim);
    const int N_cells = _acousticSolver->numVelocityCells(dim);
    Eigen::MatrixXd vertexPosition(N_cells,3);
    int count = 0;
    Vector3d vPosition;
    for (int kk=0; kk<divisions.z; kk++)
        for (int jj=0; jj<divisions.y; jj++)
            for (int ii=0; ii<divisions.x; ii++)
            {
                Tuple3i  vIndex(ii, jj, kk);
                vPosition = _acousticSolver->velocityFieldPosition(vIndex, dim);
                vertexPosition(count, 0) = vPosition[0];
                vertexPosition(count, 1) = vPosition[1];
                vertexPosition(count, 2) = vPosition[2];
                count ++;
            }

    try
    {
        IO::writeMatrixX<double>(vertexPosition, filename.c_str(), IO::BINARY);
    }
    catch (const std::runtime_error &e)
    {
    }
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
_SaveListeningPositions(const std::string &filename)
{
    Vector3Array &listeningPoints = (_acousticSolverSettings->useShell ?
            _owner->listen->speakers :
            _acousticSolverSettings->listeningPoints);
    const int N = listeningPoints.size();
    if (N==0) return;
    Eigen::MatrixXd listeningPoints_eigen(N, 3);
    for (int ii=0; ii<N; ii++)
    {
        listeningPoints_eigen(ii, 0) = listeningPoints.at(ii).x;
        listeningPoints_eigen(ii, 1) = listeningPoints.at(ii).y;
        listeningPoints_eigen(ii, 2) = listeningPoints.at(ii).z;
    }

    try
    {
        IO::writeMatrixX<double>(listeningPoints_eigen, filename.c_str(), IO::BINARY);
    }
    catch (const std::runtime_error &e)
    {
    }

    if (_owner && _acousticSolverSettings->useShell &&
        _owner->listen->mode == ListeningUnit::MODE::SHELL)
    {
        std::string file_shell;
        file_shell = CompositeFilename("out_shell_config.txt");
        {
            std::ofstream ofs(file_shell.c_str());
            _owner->listen->outShell->WriteToFile(ofs);
            ofs.close();
        }
        file_shell = CompositeFilename("inn_shell_config.txt");
        {
            std::ofstream ofs(file_shell.c_str());
            _owner->listen->innShell->WriteToFile(ofs);
            ofs.close();
        }
    }
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
_SavePressureTimestep(const std::string &filename)
{
    const Tuple3i &pFieldDivisions = GetGrid().pressureFieldDivisions();
    const int N_cells = _acousticSolver->numCells();
    int count = 0;
    std::shared_ptr<Eigen::MatrixXd> vertexPressure(new Eigen::MatrixXd(N_cells, 1));
    for (int kk=0; kk<pFieldDivisions[2]; kk++)
        for (int jj=0; jj<pFieldDivisions[1]; jj++)
            for (int ii=0; ii<pFieldDivisions[0]; ii++)
            {
                Tuple3i  vIndex( ii, jj, kk );
                VECTOR vPressure;
                _acousticSolver->vertexPressure( vIndex, vPressure );
                (*vertexPressure)(count, 0) = vPressure[0];
                count ++;
            }
    try
    {
        IO::writeMatrixX<double>(*vertexPressure, filename.c_str(), IO::BINARY);
    }
    catch (const std::runtime_error &e)
    {
    }
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
_SaveVelocityTimestep(const std::string &filename, const int &dim)
{
    const Tuple3i &divisions = _acousticSolver->velocityFieldDivisions(dim);
    const int N_cells = _acousticSolver->numVelocityCells(dim);
    int count = 0;
    std::shared_ptr<Eigen::MatrixXd> vertexVelocity_eigen(new Eigen::MatrixXd(N_cells, 1));
    for (int kk=0; kk<divisions.z; kk++)
        for (int jj=0; jj<divisions.y; jj++)
            for (int ii=0; ii<divisions.x; ii++)
            {
                Tuple3i  vIndex( ii, jj, kk );
                VECTOR vVelocity;
                _acousticSolver->vertexVelocity( vIndex, dim, vVelocity );
                (*vertexVelocity_eigen)(count, 0) = vVelocity[0];
                count ++;
            }
    try
    {
        IO::writeMatrixX<double>(*vertexVelocity_eigen, filename.c_str(), IO::BINARY);
    }
    catch (const std::runtime_error &e)
    {
    }
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
_SaveListeningData(const std::string &filename)
{
    Vector3Array &listeningPoints = (_owner ?
            _owner->listen->speakers :
            _acousticSolverSettings->listeningPoints);
    if (listeningPoints.size()>0)
    {
        Eigen::MatrixXd data;
        _acousticSolver->FetchPressureData(listeningPoints, data);
        try
        {
            IO::writeMatrixX<double>(data, filename.c_str(), IO::BINARY);
        }
        catch (const std::runtime_error &e)
        {
        }

#if DEBUG_ANALYTICAL_ACC_NOISE == 1
        const int N_points = listeningPoints.size();
        const std::string filenameAnalytical = filename + std::string(".analytical");
        std::vector<RigidSoundObjectPtr> &objects = _sceneObjects->GetRigidSoundObjects();
        for (auto &soundObject : objects)
        {
            const std::string filenameAnalyticalObject = filenameAnalytical + std::string(".") + soundObject->GetMeshName();
            Eigen::MatrixXd dataAnalytical(N_points, 1);
            for (int p_idx=0; p_idx<N_points; ++p_idx)
            {
                const REAL p_AN = soundObject->EvaluateAccelerationNoiseAnalytical(listeningPoints.at(p_idx), _simulationTime, _acousticSolverSettings->airDensity, _acousticSolverSettings->soundSpeed, 0.05);
                dataAnalytical(p_idx, 0) = p_AN;
            }
            try
            {
                IO::writeMatrixX<double>(dataAnalytical, filenameAnalyticalObject.c_str(), IO::BINARY);
            }
            catch (const std::runtime_error &e)
            {
            }
        }
#endif
    }

}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
_SaveModalFrequencies(const std::string &filename)
{
    auto &objects = _sceneObjects->GetRigidObjects();
    for (const auto &m : objects)
    {
        if (m.second->Type() == RIGID_SOUND_OBJ)
        {
            const auto &object = std::dynamic_pointer_cast<FDTD_RigidSoundObject>(m.second);
            if (object->IsModalObject())
            {
                const std::string objFilename = filename + "_" + object->GetMeshName();
                std::ofstream of(objFilename.c_str());
                of << setprecision(20);
                const int N_modes = object->N_Modes();
                for (int m_idx=0; m_idx<N_modes; ++m_idx)
                {
                    of << object->GetModeFrequency(m_idx) << "\n";
                }
                of.close();
            }
        }
    }
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
_SaveSimulationSnapshot(const std::string &numbering)
{
    const auto &solver = GetSolver();
#ifdef USE_COLLOCATED
    MATRIX p_pml[3];
    MATRIX p_pml_full;
    MATRIX v_pml[3];
    FloatArray p_gc;
    MATRIX p_collocated[3];
    int p_collocated_ind = -1;
    solver->GetAllSimulationData(p_pml, p_pml_full, v_pml, p_gc, p_collocated, p_collocated_ind);
    std::string f_p_pml, f_p_pml_full, f_v_pml, f_p_gc, f_p_collocated, f_p_collocated_ind;
    for (int dim=0; dim<3; ++dim)
    {
        std::ostringstream oss;
        oss << std::setw(1) << dim;
        f_p_pml             = CompositeFilename("snapshot_"+numbering+"_p_pml_"+oss.str()+".dat");
        f_v_pml             = CompositeFilename("snapshot_"+numbering+"_v_pml_"+oss.str()+".dat");
        f_p_collocated      = CompositeFilename("snapshot_"+numbering+"_p_collocated_"+oss.str()+".dat");
        p_pml[dim].write(f_p_pml.c_str());
        v_pml[dim].write(f_v_pml.c_str());
        p_collocated[dim].write(f_p_collocated.c_str());
    }
    f_p_pml_full        = CompositeFilename("snapshot_"+numbering+"_p_pml_full.dat");
    f_p_collocated_ind  = CompositeFilename("snapshot_"+numbering+"_p_collocated_ind.dat");
    p_pml_full.write(f_p_pml_full.c_str());
    {
        std::ofstream of(f_p_collocated_ind.c_str());
        of << p_collocated_ind;
        of.close();
    }
#endif
    ++_snapshotIndex;
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
InitializeSolver()
{
    assert(CanInitializeSolver());
    _ParseSolverSettings();
    auto p = _acousticSolverSettings->solverControlPolicy;
    assert( p && p->type == "static");
    auto policy = std::dynamic_pointer_cast<Static_Policy>(p);
    BoundingBox solverBox(_acousticSolverSettings->cellSize,
                          policy->cellDivisions,
                          policy->domainCenter);
    InitializeSolver(solverBox, _acousticSolverSettings);
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
InitializeSolver(const BoundingBox &solverBox,
                 const PML_WaveSolver_Settings_Ptr &settings)
{
    assert(CanInitializeSolver());
    if (_acousticSolverSettings != settings) _acousticSolverSettings = settings;

    // initialize solver and set various things
    _acousticSolver = std::make_shared<PML_WaveSolver>(solverBox, settings, _sceneObjects);
    _SetListeningPoints();
    //_SetBoundaryConditions();
    _SetPressureSources();

    const REAL startTime = settings->fastForwardToEventTime;  // default: 0
    ResetStartTime(startTime);

    // save settings
    SaveSolverConfig();
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
ResetStartTime(const REAL &startTime)
{
    std::cout << "\nReset system start time to " << startTime << std::endl;
    _simulationTime = startTime;
    if (_acousticSolverSettings->rigidsimDataRead)
    {
        AnimateObjects(_simulationTime);
        auto &objects = _sceneObjects->GetRigidObjects();
        for (auto &m : objects)
            if (m.second->Type() == RIGID_SOUND_OBJ &&
                std::dynamic_pointer_cast<FDTD_RigidSoundObject>(m.second)->Animated())
                std::dynamic_pointer_cast<FDTD_RigidSoundObject>(m.second)->ResetUnionBox();
    }
    _acousticSolver->Reinitialize_PML_WaveSolver(_acousticSolverSettings->useMesh, startTime);
    auto &objects = _sceneObjects->GetRigidObjects();
    if( _sceneObjects->AreEvalsDisabled() )
        _sceneObjects->EnableAllEvals();

    for (auto &m : objects)
    {
        if (m.second->Type() == RIGID_SOUND_OBJ)
        {
            auto object = std::dynamic_pointer_cast<FDTD_RigidSoundObject>(m.second);
            if (object->IsModalObject())
            {
                // take one step to update q Vectors. Note that this way all modal odes are one
                // step further than the solver because we need derivatives along
                // the way.
#ifdef JUMP_MODAL_ODES_TO_START
                object->SetODESolverTime(startTime);
#else
                REAL earliestImpulseTime = object->GetFirstImpulseTime();
                if ( earliestImpulseTime < startTime )
                {
                    REAL timeStepSize = object->GetODEStepSize();
                    int nTimeSteps = std::ceil((startTime - earliestImpulseTime) / timeStepSize );
                    // Calculate the start time by setting it back nTimeSteps timesteps.
                    // Don't do it this way, because we need to match the accumulated machine error in the ODE stepper:
                    // REAL ODEStartTime = startTime - nTimeSteps * timeStepSize;

                    // calculate ODE start time by subtracting nTimeSteps timesteps (this accurately models machine error)
                    REAL ODEStartTime = startTime;
                    for( int ind = 0; ind < nTimeSteps; ind++)
                        ODEStartTime -= timeStepSize;

                    std::cout << "Setting time back by " << nTimeSteps << " timesteps to " << ODEStartTime << std::endl;

                    object->SetODESolverTime(ODEStartTime);
                    for( int ind = 0; ind < nTimeSteps; ind++ )
                    {
                        object->AdvanceModalODESolvers(1);
                        object->UpdateQPointers();
                    }
                }
                else
                {
                    object->SetODESolverTime(startTime);
                }
#endif
                object->AdvanceModalODESolvers(1);
                object->UpdateQPointers();
            }
        }
    }
    _acousticSolver->GetGrid().ResetCellHistory(true);
    std::cout << "\nDone resetting start time!" << std::endl;
}

//##############################################################################
//##############################################################################
bool FDTD_AcousticSimulator::
RunForSteps(const int &N_steps)
{
    bool continueStepping = true;
    for (int step_idx=0; step_idx<N_steps; ++step_idx)
    {
        // first step all modal ode, so that its one step ahead of main acoustic
        // simulator. this is needed because central difference is used for
        // velocity and accleration estimates.
        const REAL odeTime = _sceneObjects->AdvanceAllModalODESolvers(1);

        // step acoustic equations
        continueStepping = _acousticSolver->stepSystem();
        PostStepping(odeTime);
    }
    return continueStepping;
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
Run()
{
    bool continueStepping = true;
    while(continueStepping)
    {
        // first step all modal ode, so that its one step ahead of main acoustic
        // simulator. this is needed because central difference is used for
        // velocity and accleration estimates.
        const REAL odeTime = _sceneObjects->AdvanceAllModalODESolvers(1);

        // step acoustic equations
        continueStepping = _acousticSolver->stepSystem();
        PostStepping(odeTime);
    }
    EndStepping();
}

//##############################################################################
//##############################################################################
bool FDTD_AcousticSimulator::
RunHalfStep(const int &flag)
{
    if (flag == 0)
    {
        _sceneObjects->AdvanceAllModalODESolvers(1);
        _acousticSolver->stepSystemHalf(flag);
    }
    else
    {
        _acousticSolver->stepSystemHalf(flag);
        PostStepping(-99);
    }
    return true;
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
PostStepping(const REAL &odeTime)
{
    auto &settings = _acousticSolverSettings;
    // file saving
    if (_stepIndex % settings->timeSavePerStep == 0)
    {
        const int timeIndex = _stepIndex/settings->timeSavePerStep;
        std::ostringstream oss;
        oss << std::setw(8) << std::setfill('0') << timeIndex;
        const std::string filenameProbe = CompositeFilename("data_listening_"+oss.str()+".dat");
        // uncomment if want to write file for each step
        //_SaveListeningData(filenameProbe);
        //if (settings->writePressureFieldToDisk)
        //{
        //    const std::string filenameField = CompositeFilename("data_pressure_"+oss.str()+".dat");
        //    _SavePressureTimestep(filenameField);
        //    // uncomment if want to store velocities
        //    //for (int dim=0; dim<3; ++dim)
        //    //{
        //    //    const std::string filenameVelocityField = CompositeFilename("velocity_"+std::to_string(dim)+"_"+oss.str()+".dat");
        //    //    _SaveVelocityTimestep(filenameVelocityField, dim);
        //    //}
        //}
    }
    std::cout << "Acoustic simulator time = " << _simulationTime
              << "; step index= " << _stepIndex
              << "; Modal ODE time = " << odeTime << std::endl;
    _stepIndex ++;
    _simulationTime += settings->timeStepSize;
    //TestMoveObjects();

#ifdef DEBUG
    _acousticSolver->PrintAllFieldExtremum();
#endif
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
EndStepping()
{
    _SaveSimulationSnapshot(std::to_string(_snapshotIndex));
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
PreviewStepping(const uint &speed)
{
    auto &settings = _acousticSolverSettings;
    std::cout << "Acoustic simulator time = " << _simulationTime << std::endl;
    _stepIndex += speed;
    _simulationTime += settings->timeStepSize*(REAL)speed;
    AnimateObjects(_simulationTime);
}


//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
SaveSolverConfig()
{
    const std::string solverSettings_s = CompositeFilename("solver_settings.txt");
    const std::string vertexPosition_s = CompositeFilename("pressure_vertex_position.dat");
    const std::string listeningPosition_s = CompositeFilename("listening_position.dat");
    const std::string modalFrequencies_s = CompositeFilename("modal_frequencies.txt");
    _SaveSolverSettings(solverSettings_s);
    _SavePressureCellPositions(vertexPosition_s);
    for (int dim=0; dim<3; ++dim)
    {
        const string velocityVertexPosition_s = CompositeFilename("velocity_"+std::to_string(dim)+"_vertex_position.dat");
        _SaveVelocityCellPositions(velocityVertexPosition_s, dim);
    }
    _SaveListeningPositions(listeningPosition_s);

    _SaveModalFrequencies(modalFrequencies_s);
}

//##############################################################################
// This function animates objects in the scene. Silent return if no stored
// rigidsim data found.
//##############################################################################
void FDTD_AcousticSimulator::
AnimateObjects(const REAL newTime)
{
    _sceneObjects->AnimateObjects(newTime);
}

//##############################################################################
// Function MoveSimBox
//   This function moves the sim box
//##############################################################################
bool FDTD_AcousticSimulator::
SetFieldCenter(const Vector3d &center)
{
    auto &field = GetGrid().pressureField();
    Tuple3i offset;
    const Vector3d nowCenter = _owner->BoundingBoxCenter();
    for (int d=0; d<3; ++d)
    {
        offset[d] = (int)((center[d] - nowCenter[d])
            /_acousticSolverSettings->cellSize);
    }
    const int l1 = abs(offset[0]) + abs(offset[1]) + abs(offset[2]);
    if (l1 > 0)
    {
        _acousticSolver->ScheduleMoveBox(offset);
        return true;
    }
    return false;
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
TestAllComponents()
{
    std::cout << "TEST BOUNDARY CONDITIONS\n";
    const int N_objects = _sceneObjects->N();
    for (int index=0; index<N_objects; ++index)
    {
        RigidObjectPtr objectPtr = _sceneObjects->GetPtr(index);
        objectPtr->TestObjectBoundaryCondition();
    }

    std::cout << "TEST PRESSURE SOURCES\n";
    const int N_sources = _sceneObjects->N_sources();
    for (int index=0; index<N_sources; ++index)
    {
        PressureSourcePtr &sourcePtr = _sceneObjects->GetPressureSourcePtr(index);
        sourcePtr->PrintSourceInfo(std::cout);
    }
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
TestMoveObjects()
{
    const int N_objects = _sceneObjects->N();
    for (int ii=0; ii<N_objects; ii++)
    {
        FDTD_RigidObject &animatedObject = _sceneObjects->Get(ii);
        animatedObject.ApplyTranslation(0.0, -2.e-5, 0.0);
        //animatedObject.ApplyTranslation(0.0, -5.6689E-6, 0.0);
        //animatedObject.PrintBoundingBox();
        //animatedObject.PrintTransformation();
    }
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
TestAnimateObjects(const int &N_steps)
{
    auto &settings = _acousticSolverSettings;
    for (int ii=0; ii<N_steps; ++ii)
    {
        _simulationTime += settings->timeStepSize;
        AnimateObjects(_simulationTime);
    }
}
