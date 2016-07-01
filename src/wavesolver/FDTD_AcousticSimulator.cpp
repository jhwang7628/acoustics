#include <wavesolver/FDTD_AcousticSimulator.h> 
#include <utils/IO/IO.h>

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
_GetSolverSettings()
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
    // FIXME debug: for now, attach a harmonic vibrational to all objects in the scene
    const int N_objects = _sceneObjects->N();
    const REAL omega = 2.0*M_PI*500.0;
    const REAL phase = 0.0;
    for (int index=0; index<N_objects; ++index)
    {
        RigidObjectPtr objectPtr = _sceneObjects->GetPtr(index);
        VibrationalSourcePtr sourcePtr(new HarmonicVibrationalSource(objectPtr, omega, phase)); 
        objectPtr->AddVibrationalSource(sourcePtr); 
        //objectPtr->TestObjectBoundaryCondition();
    }
    // TODO {
    // parser-based 
    // } TODO 
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

    IO::writeMatrixX<double>(vertexPosition, filename.c_str(), IO::BINARY);
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

    IO::writeMatrixX<double>(listeningPoints_eigen, filename.c_str(), IO::BINARY); 
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

    IO::writeMatrixX<double>(*vertexPressure, filename.c_str(), IO::BINARY);
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
        IO::writeMatrixX<double>(data, filename.c_str(), IO::BINARY);
    }
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
InitializeSolver()
{
    if (!_canInitializeSolver)
        _GetSolverSettings();

    // setup listening points
    _SetListeningPoints(); 
    // initialize solver
    _acousticSolver = std::make_shared<PML_WaveSolver>(*_acousticSolverSettings, _sceneObjects); 
    // setup source objects in the scene
    _SetBoundaryConditions();
    _SetPressureSources();
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
Run()
{
    bool continueStepping = true; 
    int stepIndex = 0; 
    SaveSolverConfig();

    while(continueStepping) 
    {
        continueStepping = _acousticSolver->stepSystem();
        if (stepIndex % _acousticSolverSettings->timeSavePerStep == 0)
        {
            const int timeIndex = stepIndex/_acousticSolverSettings->timeSavePerStep; 
            std::ostringstream oss; 
            oss << std::setw(6) << std::setfill('0') << timeIndex; 
            const std::string filenameField = _CompositeFilename("pressure_"+oss.str()+".dat"); 
            _SavePressureTimestep(filenameField); 
            const std::string filenameProbe = _CompositeFilename("listening_"+oss.str()+".dat"); 
            _SaveListeningData(filenameProbe);
        }
#ifdef DEBUG
        _acousticSolver->PrintAllFieldExtremum();
#endif
        stepIndex ++;
        _simulationTime += _acousticSolverSettings->timeStepSize; 

        // debug FIXME
        if (stepIndex > 20)
            TestMoveObjects();
    }
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
SaveSolverConfig()
{
    const string solverSettings_s = _CompositeFilename("solver_settings.txt"); 
    const string vertexPosition_s = _CompositeFilename("vertex_position.dat"); 
    const string listeningPosition_s = _CompositeFilename("listening_position.dat"); 
    _SaveSolverSettings(solverSettings_s);
    _SavePressureCellPositions(vertexPosition_s); 
    if (_acousticSolverSettings->listening)
        _SaveListeningPositions(listeningPosition_s); 
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