//////////////////////////////////////////////////////////////////////
// SampledFunction.cpp: Implementation
//
//////////////////////////////////////////////////////////////////////

#include "SampledFunction.h"

const InterpolationFunction::FunctionName SampledFunction::INTERPOLATION_TYPE
                                    = InterpolationFunction::INTERPOLATION_MITCHELL_NETRAVALI;

const int SampledFunction::INITIAL_LENGTH = 10;

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
SampledFunction::SampledFunction( int samplingRate )
    : _samplingRate( samplingRate ),
      _function( INITIAL_LENGTH, 0.0 )
{
    _h = 1.0 / (REAL)samplingRate;

    _interpolationFunction = new InterpolationMitchellNetravali( _h );
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
SampledFunction::SampledFunction()
    : _samplingRate( 44100 ),
      _function( INITIAL_LENGTH, 0.0 )
{
    _h = 1.0 / (REAL)_samplingRate;

    _interpolationFunction = new InterpolationMitchellNetravali( _h );
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
SampledFunction::SampledFunction( const SampledFunction &f )
{
    _samplingRate = f._samplingRate;
    _h = f._h;

    _interpolationFunction = new InterpolationMitchellNetravali( _h );
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
SampledFunction::~SampledFunction()
{
    delete _interpolationFunction;
}

//////////////////////////////////////////////////////////////////////
// Adds a sample to the given point in time
//////////////////////////////////////////////////////////////////////
void SampledFunction::addSample( REAL t, REAL sampleValue )
{
    int                        sampleMin;
    int                        sampleMax;

    sampleMin = (int)floor( ( t - _interpolationFunction->supportLength() ) / _h );
    sampleMax = (int)ceil( ( t + _interpolationFunction->supportLength() ) / _h );
    sampleMin = max( sampleMin, 0 );

    // Enlarge the function if necessary, padding new space with zeros
    while ( sampleMax >= _function.size() )
    {
        _function.resize( 2 * _function.size(), 0.0 );
    }

    for ( int sample_idx = sampleMin; sample_idx <= sampleMax; sample_idx++ )
    {
        REAL                     tSample = _h * (REAL)sample_idx;

        _function[ sample_idx ]
            += _interpolationFunction->evaluate( tSample, t ) * sampleValue;
    }
}

//////////////////////////////////////////////////////////////////////
// Add a sample distribued with the given interpolation
// function
//////////////////////////////////////////////////////////////////////
void SampledFunction::addSample( const InterpolationFunction *interp,
                                 REAL t, REAL sampleValue )
{
    int                        sampleMin;
    int                        sampleMax;

    sampleMin = (int)floor( ( t - interp->supportLength() ) / _h );
    sampleMax = (int)ceil( ( t + interp->supportLength() ) / _h );
    sampleMin = max( sampleMin, 0 );
    sampleMax = max( sampleMax, 0 );

    // Enlarge the function if necessary, padding new space with zeros
    while ( sampleMax >= (int)_function.size() )
    {
        _function.resize( 2 * _function.size(), 0.0 );
    }

    for ( int sample_idx = sampleMin; sample_idx <= sampleMax; sample_idx++ )
    {
        REAL                     tSample = _h * (REAL)sample_idx;

        _function[ sample_idx ] += interp->evaluate( tSample, t ) * sampleValue;
    }
}

//////////////////////////////////////////////////////////////////////
// Adds a signal defined with the given interpolation function
//////////////////////////////////////////////////////////////////////
void SampledFunction::addSignal( const InterpolationFunction *interp,
                                 const REAL *data, size_t dataSize,
                                 REAL startTime )
{
    REAL                       endTime;

    int                        startSample;
    int                        endSample;

    endTime = interp->h() * (REAL)dataSize + startTime + interp->supportLength();
    startTime -= interp->supportLength();

    startSample = (int)( startTime / _h );
    endSample = (int)( endTime / _h );

    startSample = max( startSample, 0 );

    // Enlarge the function if necessary, padding new space with zeros
    while ( endSample >= (int)_function.size() )
    {
        _function.resize( 2 * _function.size(), 0.0 );
    }

    for ( int sample_idx = startSample; sample_idx <= endSample; sample_idx++ ) {
        REAL                     t = _h * (REAL)sample_idx;

        _function[ sample_idx ] += EvaluateFunction( interp, data, dataSize,
                startTime, t );
    }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
SampledFunction& SampledFunction::operator=( const SampledFunction &f )
{
    clear();

    delete _interpolationFunction;

    _samplingRate = f._samplingRate;
    _h = f._h;

    _interpolationFunction = new InterpolationMitchellNetravali( _h ); 

    return *this;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
REAL SampledFunction::EvaluateFunction( const InterpolationFunction *interp,
                                        const REAL *data, size_t dataSize,
                                        REAL startTime, REAL t )
{
    REAL                       sampleNumReal;
    int                        sampleNum;
    int                        startSample;
    int                        endSample;

    REAL                       functionValue = 0.0;
    REAL                       sampleTime;

    sampleNumReal = ( t - startTime ) / interp->h();
    sampleNum = (int)floor( sampleNumReal );

    startSample = sampleNum - interp->support() + 1;
    endSample = sampleNum + interp->support();

    startSample = max( startSample, 0 );
    endSample = min( endSample, (int)dataSize - 1 );

    for ( int sample_idx = startSample; sample_idx <= endSample; sample_idx++ ) {
        sampleTime = startTime + interp->h() * (REAL)sample_idx;

        functionValue += interp->evaluate( t, sampleTime ) * data[ sample_idx ];
    }

    return functionValue;
}
