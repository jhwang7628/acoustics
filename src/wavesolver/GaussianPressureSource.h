#ifndef GAUSSIAN_PRESSURE_SOURCE_H 
#define GAUSSIAN_PRESSURE_SOURCE_H 

#include <TYPES.h> 
#include <config.h> 
#include <linearalgebra/Vector3.hpp> 
#include <wavesolver/PressureSource.h> 
#include <wavesolver/Wavesolver_ConstantsAndTypes.h> 

//##############################################################################
// Gaussian point pressure source
//##############################################################################
class GaussianPressureSource : public PressureSource 
{
    private: 

        Vector3d    _sourcePosition;
        
        // width of the impulse response source in space and time 
        REAL        _widthSpace;  // optional, default to sound_speed*width_tim
        REAL        _widthTime; 
        // offset of the source release. 
        REAL        _offsetTime;  // optional, default to 0.0
        // extra normalization constant for this source 
        REAL        _normalizeConstant; //optional, default to 1/(sqrt(2pi)*width_space)^3, so that integral of gaussian is normalized to 1
        // if true then flip sign of gaussian
        bool        _flipSign; // optional, default to false

    public: 
        GaussianPressureSource(){} 
        GaussianPressureSource(const Vector3d &sourcePosition, const REAL &widthSpace, const REAL &widthTime, const REAL &offsetTime, const REAL &normalizeConstant) 
            : _sourcePosition(sourcePosition),
              _widthSpace(widthSpace), 
              _widthTime(widthTime), 
              _offsetTime(offsetTime), 
              _normalizeConstant(normalizeConstant)
        {
            _bboxWorld = BoundingBox(sourcePosition-widthSpace*GAUSSIAN_CHECK_BOUND, sourcePosition+widthSpace*GAUSSIAN_CHECK_BOUND);
        }
             
        inline const Vector3d &SourcePosition() const {return _bboxWorld.centroid;} 
        virtual REAL Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time); 
        virtual void PrintSourceInfo(std::ostream &os) const; 
};

#endif 
