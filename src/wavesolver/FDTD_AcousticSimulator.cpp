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
        objectPtr->TestObjectBoundaryCondition();
    }
    // TODO {
    // parser-based 
    // } TODO 
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
        //continueStepping = _acousticSolver->stepSystem(_boundaryConditionEvaluator);  //TODO 
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
