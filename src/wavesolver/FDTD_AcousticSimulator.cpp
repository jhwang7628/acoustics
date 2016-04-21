#include <wavesolver/FDTD_AcousticSimulator.h> 

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
GetBasicSolverSettings()
{
    _parser = ImpulseResponseParser(_configFile)



    _canInitializeSolver = true;
}


//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator::
InitializeSolver()
{
    if (!_canInitializeSolver) 
        throw std::runtime_error("**ERROR** settings for the solver not ready"); 
}
