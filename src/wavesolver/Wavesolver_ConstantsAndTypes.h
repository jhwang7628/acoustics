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
//##############################################################################
#define PRINT_EIGEN_VECTOR3(x) std::cout << x[0] << ", " << x[1] << ", " << x[2] << std::endl

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
const double DISTANCE_TOLERANCE = 1E-5; 
const double GAUSSIAN_CHECK_BOUND = 3.0; // for pressure sources only check within
const REAL D_INF = std::numeric_limits<REAL>::infinity(); 

//##############################################################################
// Enum
//##############################################################################
enum VibrationalSourceType { HarmonicSource=0 };

//##############################################################################
#endif 
