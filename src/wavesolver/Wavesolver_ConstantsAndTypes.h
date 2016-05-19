#ifndef FDTD_TYPES_H 
#define FDTD_TYPES_H 

#include <parser/Parser.h>
//#include <wavesolver/FDTD_RigidObject.h> 
//#include <wavesolver/VibrationalSource.h> 

//##############################################################################
// Forward Declaration
//##############################################################################
class FDTD_RigidObject; 
class VibrationalSource; 
class PressureSource; 

//##############################################################################
// Macros
//##############################################################################
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288419716939937510582
#endif
                                                                                
#define PRINT_EIGEN_VECTOR3(x) std::cout << x[0] << ", " << x[1] << ", " << x[2] << std::endl

#define FOR_ALL_3D_GRID_VECTOR3(start,range,ii,jj,kk) \
    for (ii=start.x; ii<start.x+range.x; ++ii) \
        for (jj=start.y; jj<start.y+range.y; ++jj) \
            for (kk=start.z; kk<start.z+range.z; ++kk)
            
//##############################################################################
// Typedefs 
//##############################################################################
//typedef std::vector<Source> SourceVector; 
typedef std::shared_ptr<FDTD_RigidObject> RigidObjectPtr; 
typedef std::unique_ptr<VibrationalSource> VibrationalSourcePtr; 
typedef std::unique_ptr<PressureSource> PressureSourcePtr; 
typedef std::vector<VibrationalSourcePtr>::iterator SourceIterator; 

//##############################################################################
// Constants
//##############################################################################
//const double sqrt_pi_over_2 = sqrt(M_PI/2.0); 
const double sqrt_2_pi = sqrt(2.0*M_PI); 
const double sqrt_2 = sqrt(2.0); 
//const double DISTANCE_TOLERANCE = 0.0;
const double GAUSSIAN_CHECK_BOUND = 3.0; // for pressure sources only check within
const double AABB_CHECK_TOLERANCE_SCALE = 1.1; // scaling factor applied when checking AABB inside

const REAL D_INF = std::numeric_limits<REAL>::infinity(); 
const int I_INF = std::numeric_limits<int>::infinity(); 

//##############################################################################
// Enum
//##############################################################################
enum VibrationalSourceType { HarmonicSource=0 };

//##############################################################################
#endif 
