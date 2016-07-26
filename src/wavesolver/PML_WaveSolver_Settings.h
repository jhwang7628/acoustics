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

    // IO settings
    std::string     outputPattern; 

    // Boundary settings
    REAL    PML_width; 
    REAL    PML_strength;

    // Optional settings, mostly switches
    bool            useMesh; 
    bool            cornellBoxBoundaryCondition;
    bool            useGhostCell; 

    bool            listening; 
    std::string     listeningFile; 
    Vector3Array    listeningPoints; 

    bool            writePressureFieldToDisk; 
};

#endif
