#ifndef SIM_WORLD_H 
#define SIM_WORLD_H 
#include <set>
#include "config.h"
#include "logging/logging.h"
#include "geometry/BoundingBox.h" 
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

    struct State
    {
        REAL time = 0.; 
    } _state; 

    FDTD_Objects_Ptr                _objectCollections; 
    std::set<ActiveSimUnit_UPtr>    _simUnits; 
    PML_WaveSolver_Settings_Ptr     _simulatorSettings; 

public: 
    // Getters 
    inline const FDTD_Objects_Ptr &GetSceneObjects(){return _objectCollections;} 
    inline REAL GetWorldTime(){return _state.time;}
    inline auto GetSolverSettings(){return _simulatorSettings;}
    std::vector<BoundingBox> GetSolverBBoxs(); 
    // 
    void Build(ImpulseResponseParser_Ptr &parser); 
    bool StepWorld(); 
};
using SimWorld_UPtr = std::unique_ptr<SimWorld>; 

#endif
