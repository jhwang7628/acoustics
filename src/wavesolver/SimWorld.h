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
// Struct ActiveSimUnit
//##############################################################################
struct ActiveSimUnit
{
    FDTD_AcousticSimulator_Ptr simulator; 
    FDTD_Objects_Ptr           objects; 
    Vector3d                   boxCenter; 
    bool                       viewerUpdate = false; 
}; 
using ActiveSimUnit_UPtr = std::unique_ptr<ActiveSimUnit>; 
using ActiveSimUnit_Ptr = std::shared_ptr<ActiveSimUnit>; 

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
    struct State
    {
        REAL time = 0.25; 
    } _state; 

    FDTD_Objects_Ptr            _objectCollections; 
    std::set<ActiveSimUnit_Ptr> _simUnits; 
    PML_WaveSolver_Settings_Ptr _simulatorSettings; 

public: 
    // Getters 
    inline const FDTD_Objects_Ptr &GetSceneObjects(){return _objectCollections;} 
    inline REAL GetWorldTime(){return _state.time;}
    inline auto GetSolverSettings(){return _simulatorSettings;}
    inline const auto &GetActiveSimUnits(){return _simUnits;}
    std::vector<BoundingBox> GetSolverBBoxs(); 
    // 
    void Build(ImpulseResponseParser_Ptr &parser); 
    void UpdateObjectState(const REAL &time); 
    bool StepWorld(); 
    void PreviewStepping(const int &previewSpeed);
};
using SimWorld_UPtr = std::unique_ptr<SimWorld>; 

#endif
