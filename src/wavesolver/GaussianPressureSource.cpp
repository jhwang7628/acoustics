#include <wavesolver/GaussianPressureSource.h> 

//##############################################################################
// Use bounding box to check if evaulation is needed
//##############################################################################
REAL GaussianPressureSource::
Evaluate(const Vector3d &evaluatePosition, const Vector3d &normal, const REAL &time, const int &hintTriangle)
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
#ifdef USE_COLLOCATED 
    value *= exp(-pow((time-_offsetTime)/_widthTime, 2)/2.0); 
    value *= 1./3.6836*_normalizeConstant; // just to match the erf impl
#else
    value *= erf((time - _offsetTime)/(sqrt_2*_widthTime)) - erf(-_offsetTime/(sqrt_2*_widthTime)); 
    value *= _normalizeConstant; // just to scale it up, the wave equation is linear
#endif

    return value; 
}

//##############################################################################
//##############################################################################
REAL GaussianPressureSource::
Evaluate(const int &vertexID, const Vector3d &vertexNormal, const REAL &time)
{
    throw std::runtime_error("**ERROR** Cannot sample GaussianPressureSource using vertexID");
    return 0.0;
}

//##############################################################################
//##############################################################################
Vector3d GaussianPressureSource::
Evaluate(const int &vertexID, const REAL &time)
{
    throw std::runtime_error("**ERROR** Cannot sample GaussianPressureSource using vertexID");
    return Vector3d();
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
