#include <wavesolver/GaussianPressureSource.h> 

//##############################################################################
// Use bounding box to check if evaulation is needed
//##############################################################################
REAL GaussianPressureSource::
Evaluate(const Vector3d &evaluatePosition, const Vector3d &normal, const REAL &time)
{
    if (!_bboxWorld.Inside(evaluatePosition.x, evaluatePosition.y, evaluatePosition.z))
        return 0.0; 

    const Vector3d &sourcePosition = SourcePosition(); 
    REAL value = (evaluatePosition.x - sourcePosition.x)*(evaluatePosition.x - sourcePosition.x) 
               + (evaluatePosition.y - sourcePosition.y)*(evaluatePosition.y - sourcePosition.y) 
               + (evaluatePosition.z - sourcePosition.z)*(evaluatePosition.z - sourcePosition.z); 
    value  = -value / (2.0*_widthSpace*_widthSpace); 
    value  = exp(value); 

    // integral of a gaussian in time: 
    // http://www.wolframalpha.com/input/?i=int%28+exp%28-%28x-a%29%5E2%2F2%2Fb%5E2+%29+%29+dx
#if 0 // uncomment if use staggered grid
    value *= erf((time - _offsetTime)/(sqrt_2*_widthTime)) - erf(-_offsetTime/(sqrt_2*_widthTime)); 
    value *= _normalizeConstant; // just to scale it up, the wave equation is linear
#else
    value *= exp(-pow((time-_offsetTime)/_widthTime, 2)/2.0); 
    value *= _normalizeConstant; // just to scale it up, the wave equation is linear
#endif

    return value; 
}

//##############################################################################
//##############################################################################
void GaussianPressureSource::
PrintSourceInfo(std::ostream &os) const
{
    os << "--------------------------------------------------------------------------------\n" 
       << "Class GaussianPressureSource\n" 
       << "--------------------------------------------------------------------------------\n"
       << " source position    : " << _sourcePosition.x << ", " << _sourcePosition.y << ", " << _sourcePosition.z << "\n"
       << " width space        : " << _widthSpace << "\n"
       << " width time         : " << _widthTime << "\n"
       << " offset time        : " << _offsetTime << "\n" 
       << " normalize constant : " << _normalizeConstant << "\n"
       << "--------------------------------------------------------------------------------\n";
}
