#ifndef IMPULSE_RESPONSE_TYPES_AND_CONSTANTS
#define IMPULSE_RESPONSE_TYPES_AND_CONSTANTS

#include "parser/Parser.h"


#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288419716939937510582
#endif

//const double sqrt_pi_over_2 = sqrt(M_PI/2.0); 
const double sqrt_2_pi = sqrt(2.0*M_PI); 
const double sqrt_2 = sqrt(2.0); 

typedef Parser::ImpulseResponseParms::VolumetricSource Source; 
typedef std::vector<Source> SourceVector; 

#endif 
