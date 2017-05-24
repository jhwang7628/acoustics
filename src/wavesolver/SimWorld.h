#ifndef SIM_WORLD_H 
#define SIM_WORLD_H 
#include <set>
#include "logging/logging.h"
#include "parser/ImpulseResponseParser.h" 
#include "wavesolver/PML_WaveSolver_Settings.h" 
#include "wavesolver/FDTD_AcousticSimulator.h" 

//##############################################################################
// Class SimWorld 
//   This class manages a world in R3, and interactions between simulation boxes
//   as well as the sound displays. 
//
//   The interactions of simulation boxes can be classified as follows: 
//    - MOVE: 
//    - SPLIT: 
//    - MERGE: 
//    - FILL: 
//##############################################################################
class SimWorld
{
    struct ActiveSimUnit
    {
        FDTD_AcousticSimulator_Ptr simulator; 
        FDTD_Objects_Ptr           objects; 
    }; 
    using ActiveSimUnit_UPtr = std::unique_ptr<ActiveSimUnit>; 

    FDTD_Objects_Ptr                _objectCollections; 
    std::set<ActiveSimUnit_UPtr>     _simUnits; 
    PML_WaveSolver_Settings_Ptr     _simulatorSettings; 

public: 
    void Build(ImpulseResponseParser_Ptr &parser); 
};

#endif
