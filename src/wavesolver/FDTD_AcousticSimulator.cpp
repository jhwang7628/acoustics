#include <wavesolver/FDTD_AcousticSimulator.h> 

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
    // TODO debug: for now, attach a harmonic vibrational to all objects in the scene
    const int N_objects = _sceneObjects->N();
    const REAL omega = 2.0*M_PI*500.0;
    const REAL phase = 0.0;
    for (int index=0; index<N_objects; ++index)
    {
        RigidObjectPtr objectPtr = _sceneObjects->GetPtr(index);
        VibrationalSourcePtr sourcePtr(new HarmonicVibrationalSource(objectPtr, omega, phase)); 
        //VibrationalSourcePtr sourcePtr = std::make_unique<HarmonicVibrationalSource>(objectPtr, omega, phase); 
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
InitializeSolver()
{
    if (!_canInitializeSolver)
        _GetSolverSettings();

    _acousticSolver = std::make_shared<PML_WaveSolver>(*_acousticSolverSettings, _sceneObjects); 
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
    while(continueStepping) 
    {
        continueStepping = _acousticSolver->stepSystem();
        if (stepIndex % _acousticSolverSettings->timeSavePerStep == 0)
        {
            SaveSolverResult(); 
        }
        stepIndex++;
    }
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
SaveSolverResult()
{
    std::cout << "save routine\n"; 
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
