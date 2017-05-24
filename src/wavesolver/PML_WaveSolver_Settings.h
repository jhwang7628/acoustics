#ifndef PML_WAVESOLVER_SETTINGS_H 
#define PML_WAVESOLVER_SETTINGS_H 
#include <TYPES.h>
#include <wavesolver/Wavesolver_ConstantsAndTypes.h> 

//##############################################################################
// Stores the various settings for class PML_WaveSolver
//##############################################################################
struct PML_WaveSolver_Settings
{
    // Physical parameters
    REAL    soundSpeed; 
    REAL    airDensity; 

    // discretization settings 
    REAL    cellSize; 
    int     cellDivisions; 
    REAL    timeEnd; 
    REAL    timeStepSize; 
    int     timeSavePerStep; 
    Vector3d domainCenter;  

    // IO settings
    std::string     outputPattern; 

    // Boundary settings
    REAL PML_width; 
    REAL PML_strength;

    // Optional settings, mostly switches
    int  boundaryConditionPreset; // 0: no wall. 1: wall on +x, +y, +z, 2: wall on all but +z
    bool useMesh; 
    bool useGhostCell; 

    bool            listening; 
    std::string     listeningFile; 
    Vector3Array    listeningPoints; 

    bool            writePressureFieldToDisk; 

    // rigidsim data files, used to animate objects
    bool            rigidsimDataRead = false; 
    std::string     fileDisplacement; 
    std::string     fileVelocity; 
    std::string     fileAcceleration;
    bool            fastForwardToEarliestImpact;  
    REAL            fastForwardToEventTime; // in effect only if fastForwardToEarliestImpact = false

    // ghost cell settings
    int FV_boundarySubdivision;

    // additional options
    bool validateUsingFBem; 
};
using PML_WaveSolver_Settings_Ptr = std::shared_ptr<PML_WaveSolver_Settings>; 

#endif
