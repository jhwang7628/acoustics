#ifndef SIM_WORLD_H 
#define SIM_WORLD_H 
#include <set>
#include "config.h"
#include "logging/logging.h"
#include "geometry/BoundingBox.h" 
#include "parser/ImpulseResponseParser.h" 
#include "wavesolver/PML_WaveSolver_Settings.h" 
#include "wavesolver/FDTD_AcousticSimulator.h" 
#include "wavesolver/AudioOutput.h" 

//##############################################################################
// Struct ListeningUnit
//##############################################################################
struct ListeningUnit
{
    // speaker/mics
    Vector3Array        speakers; 
    static Vector3Array microphones;
    static REAL DelayLineScaling(const Vector3d &spk, const Vector3d &mic)
    {return 1./(spk-mic).length();}
}; 
using ListeningUnit_UPtr = std::unique_ptr<ListeningUnit>; 

//##############################################################################
// Struct ActiveSimUnit
//##############################################################################
struct ActiveSimUnit
{
    // main components
    FDTD_AcousticSimulator_Ptr simulator; 
    FDTD_Objects_Ptr           objects; 
    ListeningUnit_UPtr         listen; 

    // topology
    REAL                       lowerRadiusBound; 
    REAL                       upperRadiusBound; 
    Vector3d                   boxCenter; 
    bool                       boxCenterChanged = false; 

    // helper
    Vector3Array &UpdateSpeakers(); 
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
        REAL time = 0.26; 
    } _state; 

    FDTD_Objects_Ptr            _objectCollections; 
    std::set<ActiveSimUnit_Ptr> _simUnits; 
    PML_WaveSolver_Settings_Ptr _simulatorSettings; 

public: 
    SimWorld() = default;
    ~SimWorld()
    {
        AudioOutput::instance()->CloseStream(); 
    }

    struct WorldRasterizer
    {
        REAL cellSize; 
        Vector3d worldCenter = Vector3d(0,0,0); 
        inline Tuple3i rasterize(const Vector3d &pos)
        {
            return Tuple3i((int)((pos[0]-worldCenter[0])/cellSize), 
                           (int)((pos[1]-worldCenter[1])/cellSize),
                           (int)((pos[2]-worldCenter[2])/cellSize));
        }
        inline Vector3d cellCenter(const Tuple3i &indices)
        {
            return Vector3d(((REAL)indices[0]-0.5)*cellSize,
                            ((REAL)indices[1]-0.5)*cellSize,
                            ((REAL)indices[2]-0.5)*cellSize) + worldCenter; 
        }
    };
    static WorldRasterizer rasterizer; 

    // Getters 
    inline const FDTD_Objects_Ptr &GetSceneObjects(){return _objectCollections;} 
    inline void SetWorldTime(const REAL &t){_state.time=t;}
    inline REAL GetWorldTime(){return _state.time;}
    inline auto GetSolverSettings(){return _simulatorSettings;}
    inline const auto &GetActiveSimUnits(){return _simUnits;}
    std::vector<std::pair<ActiveSimUnit_Ptr,BoundingBox>> GetSolverBBoxs(); 
    // 
    void Build(ImpulseResponseParser_Ptr &parser); 
    void UpdateObjectState(const REAL &time); 
    bool StepWorld(); 
    void PreviewStepping(const int &previewSpeed);
};
using SimWorld_UPtr = std::unique_ptr<SimWorld>; 

#endif
