#ifndef WAVESOLVER_VOLUMETRIC_SOURCE_H 
#define WAVESOLVER_VOLUMETRIC_SOURCE_H 

#include <TYPES.h> 
#include <wavesolver/Source.h> 
#include <linearalgebra/Vector3.hpp> 

//##############################################################################
// Volumetric source
//##############################################################################
class VolumetricSource : public Source
{

    // TODO 
    public: 

        // position of the impulse response source
        Vector3d    position;  

        // width of the impulse response source in space and time 
        REAL        widthTime; 
        REAL        widthSpace;  // optional, default to sound_speed*width_time

        // offset of the source release. 
        REAL        offsetTime;  // optional, default to 0.0

        // extra normalization constant for this source 
        REAL        normalizeConstant; //optional, default to 1/(sqrt(2pi)*width_space)^3, so that integral of gaussian is normalized to 1

        // if true then flip sign of gaussian
        bool        flipSign; // optional, default to false
};

#endif 
