#ifndef SIM_WORLD_AUX_DS_H
#define SIM_WORLD_AUX_DS_H
#include "config.h" 
#include "wavesolver/PML_WaveSolver_Settings.h" 
#include "wavesolver/FDTD_AcousticSimulator.h" 
#include "wavesolver/FDTD_Objects.h"
#include "geometry/BoundingBox.h"
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
    int      divisions; 
    REAL     lowerRadiusBound; 
    REAL     upperRadiusBound; 
    Vector3d unitCenter; // might not be the pressure field center
    bool     boxCenterChanged = false; 

    // ids
    std::string               *unitID; 

    // helper
    Vector3Array &UpdateSpeakers(); 
    const BoundingBox GetBoundingBox()
    {return simulator->GetSolver()->GetSolverBBox();}
    const Vector3d BoundingBoxCenter()
    {return GetBoundingBox().center();}
}; 
using ActiveSimUnit_UPtr = std::unique_ptr<ActiveSimUnit>; 
using ActiveSimUnit_Ptr = std::shared_ptr<ActiveSimUnit>; 

//##############################################################################
// Struct BoundaryCell
//##############################################################################
struct BoundaryCell
{
    ActiveSimUnit_Ptr owner; 
    Tuple3i           cellIndices; 
}; 
using BoundaryCell_Ptr = std::shared_ptr<BoundaryCell>; 

#endif
