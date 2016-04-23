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

//##############################################################################
// Macros
//##############################################################################
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288419716939937510582
#endif

//##############################################################################
// Typedefs 
//##############################################################################
//typedef std::vector<Source> SourceVector; 
typedef std::shared_ptr<FDTD_RigidObject> RigidObjectPtr; 
typedef std::unique_ptr<VibrationalSource> VibrationalSourcePtr; 
typedef std::vector<VibrationalSourcePtr>::iterator SourceIterator; 

//##############################################################################
// Constants
//##############################################################################
//const double sqrt_pi_over_2 = sqrt(M_PI/2.0); 
const double sqrt_2_pi = sqrt(2.0*M_PI); 
const double sqrt_2 = sqrt(2.0); 
const double DISTANCE_TOLERANCE = 1E-12; 

//##############################################################################
// Enum
//##############################################################################
enum VibrationalSourceType { HarmonicSource=0 };

//##############################################################################
#endif 
