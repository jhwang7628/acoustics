#ifndef FDTD_TYPES_H
#define FDTD_TYPES_H

#include <Eigen/Dense>
#include <parser/Parser.h>
#include <math/MLSModeInterpolator.hpp>

//##############################################################################
// Forward Declaration
//##############################################################################
class FDTD_RigidObject;
class FDTD_RigidSoundObject;
class VibrationalSource;
class PressureSource;
class PML_WaveSolver_Settings;

//##############################################################################
// Macros
//##############################################################################
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288419716939937510582
#endif

#define PRINT_EIGEN_VECTOR3(x) std::cout << x[0] << ", " << x[1] << ", " << x[2] << std::endl
#define RAD_2_DEG(x) (x*180.0/M_PI)
#define DEG_2_RAD(x) (x/180.0*M_PI)

#define FOR_ALL_3D_GRID_VECTOR3(start,range,ii,jj,kk) \
    for (ii=start.x; ii<start.x+range.x; ++ii) \
        for (jj=start.y; jj<start.y+range.y; ++jj) \
            for (kk=start.z; kk<start.z+range.z; ++kk)

//##############################################################################
// Options
//##############################################################################
#define USE_ADF // use adaptive distance field
#define USE_COLLOCATED
//#define NO_INTERNAL_DAMPING
//#define USE_FV
//#define USE_MERGE_SIM_BOXES // TODO this does not work well, instability grows

// disable the sharp interface and use fully rasterized approach with
// staggered acceleration sampling. Seems to give better stability in some cases
//#define FULL_RASTERIZE
#ifdef FULL_RASTERIZE
#define ONLY_SAMPLE_CENTER
#define STAGGER_ACC_SAMPLE
#endif


//##############################################################################
// Typedefs
//##############################################################################
//typedef std::vector<Source> SourceVector;
typedef std::shared_ptr<FDTD_RigidObject> RigidObjectPtr;
typedef std::shared_ptr<FDTD_RigidSoundObject> RigidSoundObjectPtr;
typedef std::shared_ptr<VibrationalSource> VibrationalSourcePtr;
typedef std::shared_ptr<PressureSource> PressureSourcePtr;
typedef std::vector<VibrationalSourcePtr>::iterator SourceIterator;
typedef std::shared_ptr<PML_WaveSolver_Settings> PML_WaveSolver_Settings_Ptr;

typedef MLSModeInterpolator<REAL,3,1> T_MLS;
typedef T_MLS::MLSPoint MLSPoint;
typedef T_MLS::MLSVal MLSVal;
typedef T_MLS::MLSMatrix MLSMatrix;
typedef T_MLS::MLSVector MLSVector;
typedef Eigen::aligned_allocator<MLSPoint> P_ALLOCATOR;
typedef Eigen::aligned_allocator<MLSVal> V_ALLOCATOR;

//##############################################################################
// Constants
//##############################################################################
const double SMALL_NUM = 1E-12;
const double sqrt_2_pi = sqrt(2.0*M_PI);
const double sqrt_2 = sqrt(2.0);
const double DISTANCE_TOLERANCE = 0;
const double GAUSSIAN_CHECK_BOUND = 2.0; // for pressure sources only check within
const double AABB_CHECK_TOLERANCE_SCALE = 1.1; // scaling factor applied when checking AABB inside
const double MODE_SHAPE_CUTOFF_FREQ = 20000.0;
const double GHOST_CELL_ENTRY_THRESHOLD = 1E-7;
const double IMPULSE_VEL_THRESHOLD = 0.05;
const double INTERPOLATION_DIFF_TOL = 20.0;

const REAL D_INF = std::numeric_limits<REAL>::infinity();
const int I_INF = std::numeric_limits<int>::infinity();

#ifdef USE_ADF
//const double ADF_SUBDIVIDE_RADIUS = 0.0025*sqrt(3.0); // should be h/4
const double ADF_SUBDIVIDE_RADIUS = 0.0005*sqrt(3.0);
const int ADF_MAX_OCTREE_LEVELS = 9;
const double ADF_ERROR_TOLERANCE = 0.000001;
#endif

const int GHOST_CELL_JACOBI_MAX_ITERATION = 200;
//const REAL KD_NEAREST_TOLERANCE = 1E-6;
const REAL KD_NEAREST_TOLERANCE = 0.;
const REAL TRI_NORMAL_PUSH_DIST = 0.0005;
const double CONDN_LIM = 25;

// debug
const REAL CONTACT_TIMESCALE_SCALE = 1.;

//##############################################################################
// Enum
//##############################################################################
enum VibrationalSourceType { HarmonicSource=0 };

//##############################################################################
// Structs
//##############################################################################
struct KinematicsMetadata
{
    std::string objId;
    std::string fileKinematics;
    REAL stepSize;
};

struct ChunkPartitionParam
{
    bool adaptive = false;
    int N_0;
    int N_maxc;
    REAL L_z;
    std::string outFile;
};
using ChunkPartitionParam_Ptr =
    std::shared_ptr<ChunkPartitionParam>;

//##############################################################################
// Debugging flags, settings
//##############################################################################
#define DEBUG_ANALYTICAL_ACC_NOISE 0
#define DEBUG_PERFECT_MODAL_HARMONICS 1 // replace q(t) by cos(omega t)
#define DEBUG_WRITE_REFLECTION_ARROWS_INTERVAL -1
//##############################################################################
#endif
