#include <wavesolver/FDTD_AcousticSimulator.h> 
#include <utils/IO/IO.h>
#include <macros.h>

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
_ParseSolverSettings()
{
    _parser = std::make_shared<ImpulseResponseParser>(_configFile);
    _acousticSolverSettings = std::make_shared<PML_WaveSolver_Settings>(); 
    _sceneObjects = std::make_shared<FDTD_Objects>();
    _parser->GetObjects(_sceneObjects); 
    _parser->GetSolverSettings(_acousticSolverSettings); 
    _canInitializeSolver = true;
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
_SetBoundaryConditions()
{
    const int N_objects = _sceneObjects->N();
    for (int index=0; index<N_objects; ++index)
    {
        RigidSoundObjectPtr objectPtr = _sceneObjects->GetPtr(index);
        // add modal vibrational source // FIXME debug
        //VibrationalSourcePtr sourcePtr(new ModalVibrationalSource(objectPtr)); 
        //objectPtr->AddVibrationalSource(sourcePtr); 

        // add acceleration noise source
        VibrationalSourcePtr anSourcePtr(new AccelerationNoiseVibrationalSource(objectPtr)); 
        objectPtr->AddVibrationalSource(anSourcePtr);

        // add debug harmonic source
        //const REAL omega = 2.0*M_PI*500.0;
        //const REAL phase = 0.0;
        ////VibrationalSourcePtr sourcePtr(new HarmonicVibrationalSource(objectPtr, omega, phase, 100., 0.301997)); 
        //VibrationalSourcePtr sourcePtr(new HarmonicVibrationalSource(objectPtr, omega, phase)); 
        //objectPtr->AddVibrationalSource(sourcePtr); 
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
    Vector3Array &listeningPoints = _acousticSolverSettings->listeningPoints; 
    _parser->GetListeningPoints(listeningPoints); 
    if (listeningPoints.size() == 0) 
    {
        _acousticSolverSettings->listening = false; 
    }
    else 
    {
        _acousticSolverSettings->listening = true; 
        _acousticSolverSettings->listeningFile = _CompositeFilename("listening_pressure_%6u.dat"); 
    }
}

//##############################################################################
//##############################################################################
std::string FDTD_AcousticSimulator::
_CompositeFilename(const std::string filename)
{
    char buffer[512];
    snprintf(buffer, 512, _acousticSolverSettings->outputPattern.c_str(), filename.c_str()); 
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

    std::string xmlFilename_s = _CompositeFilename("input_control_file.xml"); 
    IO::CopyFile(_configFile, xmlFilename_s); 
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
_SavePressureCellPositions(const std::string &filename)
{
    const int N_cellsEachDimension = _acousticSolverSettings->cellDivisions;
    const int N_cells = _acousticSolver->numCells();
    Eigen::MatrixXd vertexPosition(N_cells,3); 
    int count = 0;
    Vector3d vPosition; 
    for (int kk=0; kk<N_cellsEachDimension; kk++)
        for (int jj=0; jj<N_cellsEachDimension; jj++)
            for (int ii=0; ii<N_cellsEachDimension; ii++)
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
    Vector3Array &listeningPoints = _acousticSolverSettings->listeningPoints; 
    const int N = listeningPoints.size(); 
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
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
_SavePressureTimestep(const std::string &filename)
{
    const int N_cellsEachDimension = _acousticSolverSettings->cellDivisions;
    const int N_cells = _acousticSolver->numCells();
    int count = 0;
    std::shared_ptr<Eigen::MatrixXd> vertexPressure(new Eigen::MatrixXd(N_cells, 1)); 
    for (int kk=0; kk<N_cellsEachDimension; kk++)
        for (int jj=0; jj<N_cellsEachDimension; jj++)
            for (int ii=0; ii<N_cellsEachDimension; ii++)
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
    if (_acousticSolverSettings->listening)
    {
        Vector3Array &listeningPoints = _acousticSolverSettings->listeningPoints; 
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
    auto &objects = _sceneObjects->GetRigidSoundObjects(); 
    for (const auto &object : objects)
    {
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

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
InitializeSolver()
{
    if (!_canInitializeSolver)
        _ParseSolverSettings();

    // if rigidsim data exists, read and apply them to objects before
    // initializing solver.
    const auto &settings = _acousticSolverSettings; 
    if (settings->rigidsimDataRead)
    {
        _sceneObjectsAnimator = std::make_shared<FDTD_RigidObject_Animator>(); 
        _sceneObjectsAnimator->ReadAllKinematics(settings->fileDisplacement, settings->fileVelocity, settings->fileAcceleration); 
        AnimateObjects(); // apply the transformation right away
    }
    else 
    {
        std::cerr << "**WARNING** no rigidsim data read\n";
    }

    // initialize solver and set various things
    _SetListeningPoints(); 
    _acousticSolver = std::make_shared<PML_WaveSolver>(_acousticSolverSettings, _sceneObjects); 
    _SetBoundaryConditions();
    _SetPressureSources();

    REAL startTime = 0.0; 
    // if no pressure sources found, get the earliest impact event and reset/shift all solver time to that event
    if (!_sceneObjects->HasExternalPressureSources() && _acousticSolverSettings->fastForwardToEarliestImpact)
        startTime = _sceneObjects->GetEarliestImpactEvent() - _acousticSolverSettings->timeStepSize; 
    else 
        startTime = _acousticSolverSettings->fastForwardToEventTime; 
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
    _stepIndex = 0; 
    _simulationTime = startTime; 
    if (_acousticSolverSettings->rigidsimDataRead)
    {
        AnimateObjects();
        auto &objects = _sceneObjects->GetRigidSoundObjects(); 
        for (auto &object : objects) 
            if (object->IsModalObject())
                object->ResetUnionBox();
    }
    _acousticSolver->Reinitialize_PML_WaveSolver(_acousticSolverSettings->useMesh, startTime); 
    auto &objects = _sceneObjects->GetRigidSoundObjects(); 
    for (auto &object : objects) 
    {
        if (object->IsModalObject())
        {
            // take one step to update q Vectors. Note that this way all modal odes are one
            // step further than the solver because we need derivatives along
            // the way.
            object->SetODESolverTime(startTime);
            object->AdvanceModalODESolvers(1); 
            object->UpdateQPointers(); 
        }
    }
    _acousticSolver->GetGrid().ResetCellHistory(true);
}

//##############################################################################
//##############################################################################
bool FDTD_AcousticSimulator::
RunForSteps(const int &N_steps)
{
    bool continueStepping = true; 
    auto &settings = _acousticSolverSettings; 

    for (int step_idx=0; step_idx<N_steps; ++step_idx)
    {
        // first step all modal ode, so that its one step ahead of main acoustic
        // simulator. this is needed because central difference is used for
        // velocity and accleration estimates. 
        const REAL odeTime = _sceneObjects->AdvanceAllModalODESolvers(1); 

        // step acoustic equations
        continueStepping = _acousticSolver->stepSystem();
        if (_stepIndex % settings->timeSavePerStep == 0)
        {
            const int timeIndex = _stepIndex/settings->timeSavePerStep; 
            std::ostringstream oss; 
            oss << std::setw(8) << std::setfill('0') << timeIndex; 
            const std::string filenameProbe = _CompositeFilename("data_listening_"+oss.str()+".dat"); 
            _SaveListeningData(filenameProbe);
            if (settings->writePressureFieldToDisk)
            {
                const std::string filenameField = _CompositeFilename("data_pressure_"+oss.str()+".dat"); 
                _SavePressureTimestep(filenameField); 
                // uncomment if want to store velocities
                //for (int dim=0; dim<3; ++dim) 
                //{
                //    const std::string filenameVelocityField = _CompositeFilename("velocity_"+std::to_string(dim)+"_"+oss.str()+".dat"); 
                //    _SaveVelocityTimestep(filenameVelocityField, dim); 
                //}
            }
        }
        // update modal vectors for the next time step
        for (int obj_idx=0; obj_idx<_sceneObjects->N(); ++obj_idx)
            _sceneObjects->GetPtr(obj_idx)->UpdateQPointers(); 
        std::cout << "Acoustic simulator time = " << _simulationTime << "; Modal ODE time = " << odeTime << std::endl;
        _stepIndex ++;
        _simulationTime += settings->timeStepSize; 

        AnimateObjects(); 

        if (!continueStepping)
            break; 
    }
    return continueStepping; 
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
Run()
{
    bool continueStepping = true; 
    auto &settings = _acousticSolverSettings; 

    int count = 0;
    while(continueStepping) 
    {
        // first step all modal ode, so that its one step ahead of main acoustic
        // simulator. this is needed because central difference is used for
        // velocity and accleration estimates. 
        const REAL odeTime = _sceneObjects->AdvanceAllModalODESolvers(1); 

        // step acoustic equations
        continueStepping = _acousticSolver->stepSystem();
        if (_stepIndex % settings->timeSavePerStep == 0)
        {
            const int timeIndex = _stepIndex/settings->timeSavePerStep; 
            std::ostringstream oss; 
            oss << std::setw(8) << std::setfill('0') << timeIndex; 
            const std::string filenameProbe = _CompositeFilename("data_listening_"+oss.str()+".dat"); 
            _SaveListeningData(filenameProbe);
            if (settings->writePressureFieldToDisk)
            {
                const std::string filenameField = _CompositeFilename("data_pressure_"+oss.str()+".dat"); 
                _SavePressureTimestep(filenameField); 
                // uncomment if want to store velocities
                //for (int dim=0; dim<3; ++dim) 
                //{
                //    const std::string filenameVelocityField = _CompositeFilename("velocity_"+std::to_string(dim)+"_"+oss.str()+".dat"); 
                //    _SaveVelocityTimestep(filenameVelocityField, dim); 
                //}
            }
        }
        // update modal vectors for the next time step
        for (int obj_idx=0; obj_idx<_sceneObjects->N(); ++obj_idx)
            _sceneObjects->GetPtr(obj_idx)->UpdateQPointers(); 

        std::cout << "Acoustic simulator time = " << _simulationTime << "; Modal ODE time = " << odeTime << std::endl;
        _stepIndex ++;
        _simulationTime += settings->timeStepSize; 

        AnimateObjects(); 

        count ++; 

#ifdef DEBUG
        _acousticSolver->PrintAllFieldExtremum();
#endif
    }
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
SaveSolverConfig()
{
    const std::string solverSettings_s = _CompositeFilename("solver_settings.txt"); 
    const std::string vertexPosition_s = _CompositeFilename("pressure_vertex_position.dat"); 
    const std::string listeningPosition_s = _CompositeFilename("listening_position.dat"); 
    const std::string modalFrequencies_s = _CompositeFilename("modal_frequencies.txt");
    _SaveSolverSettings(solverSettings_s);
    _SavePressureCellPositions(vertexPosition_s); 
    for (int dim=0; dim<3; ++dim) 
    {
        const string velocityVertexPosition_s = _CompositeFilename("velocity_"+std::to_string(dim)+"_vertex_position.dat"); 
        _SaveVelocityCellPositions(velocityVertexPosition_s, dim); 
    }
    if (_acousticSolverSettings->listening)
        _SaveListeningPositions(listeningPosition_s); 

    _SaveModalFrequencies(modalFrequencies_s); 
}

//##############################################################################
// This function animates objects in the scene. Silent return if no stored
// rigidsim data found.
//##############################################################################
void FDTD_AcousticSimulator::
AnimateObjects()
{
    if (_sceneObjectsAnimator) 
    {
        Vector3d displacement; 
        Vector3d rotationAxis; 
        REAL rotationAngle;
        Quaternion<REAL> quaternion; 
        for (int obj_idx=0; obj_idx<_sceneObjects->N(); ++obj_idx)
        {
            const int rigidsimObjectID = std::stoi(_sceneObjects->GetMeshName(obj_idx)); 
            if (_sceneObjects->GetPtr(rigidsimObjectID)->IsModalObject())
            {
                _sceneObjectsAnimator->GetObjectDisplacement(rigidsimObjectID, _simulationTime, displacement, quaternion); 
                rotationAngle = quaternion.toAxisRotR(rotationAxis); 
                _sceneObjects->GetPtr(obj_idx)->SetTransform(displacement.x, displacement.y, displacement.z, rotationAngle, rotationAxis.x, rotationAxis.y, rotationAxis.z); 

#ifdef DEBUG_PRINT
                std::cout << "object " << obj_idx << " has translation = " << _sceneObjects->GetPtr(obj_idx)->GetTranslation().transpose() << std::endl;
#endif
            }
        }
    }
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
        AnimateObjects(); 
    }
}
