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

// use adaptive distance field
#define USE_ADF
            
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

typedef MLSModeInterpolator<double,3,3> MLSInterpolatorType; 
typedef MLSInterpolatorType::MLSPoint MLSPoint; 
typedef MLSInterpolatorType::MLSVal MLSVal; 
typedef MLSInterpolatorType::MLSMatrix MLSMatrix; 
typedef MLSInterpolatorType::MLSVector MLSVector;
typedef Eigen::aligned_allocator<MLSPoint> P_ALLOCATOR; 
typedef Eigen::aligned_allocator<MLSVal> V_ALLOCATOR; 

//##############################################################################
// Constants
//##############################################################################
//const double sqrt_pi_over_2 = sqrt(M_PI/2.0); 
const double sqrt_2_pi = sqrt(2.0*M_PI); 
const double sqrt_2 = sqrt(2.0); 
const double DISTANCE_TOLERANCE = 0;
const double GAUSSIAN_CHECK_BOUND = 3.0; // for pressure sources only check within
const double AABB_CHECK_TOLERANCE_SCALE = 1.1; // scaling factor applied when checking AABB inside
const double MODE_SHAPE_CUTOFF_FREQ = 44100.0;

const REAL D_INF = std::numeric_limits<REAL>::infinity(); 
const int I_INF = std::numeric_limits<int>::infinity(); 

#ifdef USE_ADF
const double ADF_SUBDIVIDE_RADIUS = 0.0025*sqrt(3.0);
const int ADF_MAX_OCTREE_LEVELS = 9; 
const double ADF_ERROR_TOLERANCE = 0.00001;
#endif

//##############################################################################
// Enum
//##############################################################################
enum VibrationalSourceType { HarmonicSource=0 };

//##############################################################################
#endif 
