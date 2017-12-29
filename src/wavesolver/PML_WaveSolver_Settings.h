#ifndef PML_WAVESOLVER_SETTINGS_H 
#define PML_WAVESOLVER_SETTINGS_H 
#include <TYPES.h>
#include <wavesolver/Wavesolver_ConstantsAndTypes.h> 


//##############################################################################
// Forward declaration
//##############################################################################
struct Solver_Control_Policy; 
struct Static_Policy; 
struct Dynamic_Policy; 
 
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
    REAL    timeEnd; 
    REAL    timeStepSize; 
    int     timeSavePerStep; 
    int     numberTimeSteps = -1; 

    // IO settings
    std::string     outputPattern; 

    // Boundary settings
    REAL PML_width; 
    REAL PML_strength;

    // listening shells 
    bool useShell = false; 
    std::string refShellFile; 
    REAL spacing; 

    // damping due to viscosity of air
    // see paper: computing room acoustics with CUDA - 3D FDTD schemes with boundary losses and viscosity
    // alpha = 0.0 results in lossless media
    bool useAirViscosity = false;
    REAL alpha;

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

    // kinematic scripting files from blender, this should be disjoint with rigidsimDataRead
    bool            kinFileExists = false; 
    std::map<std::string, KinematicsMetadata> objKinematicsMetadata; 
    bool            onlyObjSequence = false; 

    bool            fastForwardToEarliestImpact;  
    REAL            fastForwardToEventTime; // in effect only if fastForwardToEarliestImpact = false

    // time parallelization 
    bool            timeParallel = false;
    int             indTimeChunks = 0;
    int             numTimeChunks = 1;
    REAL            modalParallelErrorTol;// = 1e-7;
    REAL            modalMinDampingPeriods;// = 2;
    REAL            overlapTime; // = 6e-3
    REAL            stopBoundaryAccTime; // to be set

    // ghost cell settings
    int FV_boundarySubdivision;

    // additional options
    bool validateUsingFBem; 
    enum BoundaryHandling
    {
        RASTERIZE = 0,
        FULLY_COUPLED = 1
    } boundaryHandlingType;

    std::shared_ptr<Solver_Control_Policy> solverControlPolicy; 
};
using PML_WaveSolver_Settings_Ptr = std::shared_ptr<PML_WaveSolver_Settings>; 

//##############################################################################
// Struct Solver_Control_Policy
//##############################################################################
struct Solver_Control_Policy
{
    std::string type; 
    virtual ~Solver_Control_Policy(){}
};
using Solver_Control_Policy_Ptr = std::shared_ptr<Solver_Control_Policy>;

//##############################################################################
// Struct Static_Policy
//##############################################################################
struct Static_Policy : public Solver_Control_Policy
{
    int      cellDivisions; 
    Vector3d domainCenter; 
};

//##############################################################################
// Struct Dynamic_Policy
//##############################################################################
struct Dynamic_Policy : public Solver_Control_Policy
{
    int padding; 
};
#endif
