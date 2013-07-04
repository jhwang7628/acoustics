//////////////////////////////////////////////////////////////////////
// InterpolationFunction.cpp: Implementation
//
//////////////////////////////////////////////////////////////////////

#include "InterpolationFunction.h"

#include <utils/trace.h>

#include <boost/bind.hpp>

#include <iostream>

using namespace std;

//////////////////////////////////////////////////////////////////////
// InterpolationFunction implementation
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
InterpolationFunction::InterpolationFunction( REAL h,
                                              int support,
                                              REAL supportLength )
    : _h( h ),
      _support( support ),
      _supportLength( supportLength ),
      _normalizationCoefficient( 1.0 )
{
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
InterpolationFunction::~InterpolationFunction()
{
}

//////////////////////////////////////////////////////////////////////
// Looks up the enumerator function name
//////////////////////////////////////////////////////////////////////
InterpolationFunction::FunctionName
InterpolationFunction::GetFunctionName( string functionName )
{
    if ( functionName == "mitchell-netravali" )
    {
        return INTERPOLATION_MITCHELL_NETRAVALI;
    }
    else if ( functionName == "perlin-improved" )
    {
        return PERLIN_IMPROVED;
    }
    else if ( functionName == "bspline-cubic" )
    {
        return BSPLINE_CUBIC;
    }
    else
    {
        cerr << "Warning: Interpolation function name not found.  "
            "Defaulting to Mitchell-Netravali" << endl;

        return INTERPOLATION_MITCHELL_NETRAVALI;
    }
}

//////////////////////////////////////////////////////////////////////
// Returns the requested interpolation function
//////////////////////////////////////////////////////////////////////
InterpolationFunction *InterpolationFunction::GetFunction( FunctionName name,
                                                           REAL h,
                                                           REAL supportLength )
{
    switch ( name )
    {
        case INTERPOLATION_POLY6:
            {
                return new InterpolationPoly6( h, supportLength );
            }
        case INTERPOLATION_MITCHELL_NETRAVALI:
        default:
            {
                return new InterpolationMitchellNetravali( h );
            }
        case PERLIN_IMPROVED:
            {
                return new InterpolationImprovedPerlin( h );
            }
        case BSPLINE_CUBIC:
            {
                return new InterpolationBSplineCubic( h );
            }
        case CATMULL_ROM_CUBIC:
            {
                return new InterpolationCatmullRomCubic( h );
            }
    }
}

//////////////////////////////////////////////////////////////////////
// Returns a Function object for the specified interpolation function
//////////////////////////////////////////////////////////////////////
Function InterpolationFunction::GetFunctionContainer( FunctionName name,
                                                      REAL h,
                                                      REAL supportLength )
{
    Function       function;

    switch ( name )
    {
        case INTERPOLATION_MITCHELL_NETRAVALI:
            {
                InterpolationMitchellNetravali interpFunction( h );

                function.function() = boost::bind(
                        &InterpolationMitchellNetravali::evaluate,
                        interpFunction, _1, 0.0 );

                function.supportSamples() = 2;
                function.support() = 2.0 * h;

                break;
            }
        default:
            {
                TRACE_ASSERT( 0, "Function type not supported" );

                break;
            }
#if 0
        case PERLIN_IMPROVED:
            {
            }
        case BSPLINE_CUBIC:
            {
            }
        case CATMULL_ROM_CUBIC:
            {
            }
#endif
    }

    return function;
}

//////////////////////////////////////////////////////////////////////
// Finds a normalization coefficient for this function (so that
// evaluations at distances of _h sum up to 1).
//////////////////////////////////////////////////////////////////////
void InterpolationFunction::findNormalization()
{
    REAL           total = 0.0;

    for ( int i = -_support; i <= _support; i++ )
    {
        total += evaluate( _h * (REAL)i, 0.0 );
    }

    _normalizationCoefficient = 1.0 / total;

    cout << "Got normalization " << _normalizationCoefficient << endl;
}

//////////////////////////////////////////////////////////////////////
// InterpolationPoly6 implementation
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
InterpolationPoly6::InterpolationPoly6( REAL h, REAL supportLength )
    : InterpolationFunction( h, (int)ceil( supportLength ),
                             supportLength * h )
{
    TRACE_ASSERT( supportLength > 0, "Invalid support length" );
}

//////////////////////////////////////////////////////////////////////
// t is the evaluation time, tSample is the location in
// time of the sample with which this function is associated.
//////////////////////////////////////////////////////////////////////
REAL InterpolationPoly6::evaluate( REAL t, REAL tSample ) const
{
    REAL functionValue;

    functionValue = t - tSample;

    if ( abs( functionValue ) > _supportLength )
    {
        return 0.0;
    }

    functionValue /= _supportLength;
    functionValue *= functionValue;
    functionValue = 1.0 - functionValue;
    functionValue *= functionValue * functionValue;

    return functionValue;
}

//////////////////////////////////////////////////////////////////////
// Same convention as above - evaluates the derivative
//////////////////////////////////////////////////////////////////////
REAL InterpolationPoly6::derivative( REAL t, REAL tSample ) const
{
    REAL derivativeValue;

    derivativeValue = t - tSample;

    if ( abs( derivativeValue ) > _supportLength )
    {
        return 0.0;
    }

    derivativeValue /= _supportLength;
    derivativeValue *= derivativeValue;
    derivativeValue = 1.0 - derivativeValue;
    derivativeValue *= derivativeValue;
    derivativeValue *= -6.0;
    derivativeValue /= _supportLength * _supportLength;
    derivativeValue *= t - tSample;

    return derivativeValue;
}

//////////////////////////////////////////////////////////////////////
// Generates a plot of this interpolation function
//////////////////////////////////////////////////////////////////////
void InterpolationPoly6::plot( FloatArray &fPlot, FloatArray &dfPlot,
                               REAL div ) const
{
    int ndivs = (int)ceil( 2.0 * _supportLength / div );

    fPlot.resize( ndivs );
    dfPlot.resize( ndivs );

    for ( int i = 0; i < ndivs; i++ )
    {
        fPlot[i] = evaluate( -1.0 * _supportLength + div * (REAL)i, 0.0 );
        dfPlot[i] = derivative( -1.0 * _supportLength + div * (REAL)i, 0.0 );
    }
}

//////////////////////////////////////////////////////////////////////
// InterpolationMitchellNetravali implementation
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
InterpolationMitchellNetravali::InterpolationMitchellNetravali( REAL h )
    : InterpolationFunction( h, 2, 2.0 * h )
{
}

//////////////////////////////////////////////////////////////////////
// t is the evaluation time, tSample is the location in
// time of the sample with which this function is associated.
//////////////////////////////////////////////////////////////////////
REAL InterpolationMitchellNetravali::evaluate( REAL t, REAL tSample ) const
{
    REAL arg;
    REAL functionValue = 0.0;

    arg = abs( t - tSample ) / _h;

    if ( arg >= 2.0 )
    {
        return 0.0;
    }

    if ( arg < 1.0 )
    {
        arg = 1.0 - arg;

        functionValue = -15.0 * arg;
        functionValue = arg * ( 18.0 + functionValue );
        functionValue = arg * ( 9.0 + functionValue );
        functionValue = 2.0 + functionValue;
    }
    else
    {
        arg = 2.0 - arg;

        functionValue = 5.0 * arg;
        functionValue = arg * arg * ( functionValue - 3.0 );
    }

    functionValue /= 18.0;

    return functionValue;
}

//////////////////////////////////////////////////////////////////////
// Same convention as above - evaluates the derivative
//////////////////////////////////////////////////////////////////////
REAL InterpolationMitchellNetravali::derivative( REAL t, REAL tSample ) const
{
    REAL arg;
    REAL derivativeValue = 0.0;
    REAL argSign;

    arg = abs( t - tSample ) / _h;

    argSign = ( t - tSample ) > 0.0 ? 1.0 : -1.0;

    if ( arg >= 2.0 )
    {
        return 0.0;
    }

    if ( arg < 1.0 )
    {
        arg = 1.0 - arg;

        derivativeValue = -45.0 * arg;
        derivativeValue = arg * ( 36.0 + derivativeValue );
        derivativeValue = 9.0 + derivativeValue;
    }
    else
    {
        arg = 2.0 - arg;

        derivativeValue = 15.0 * arg;
        derivativeValue = arg * ( derivativeValue - 6.0 );
    }

    derivativeValue *= -1.0 * argSign;
    derivativeValue /= _h * 18.0;

    return derivativeValue;
}

//////////////////////////////////////////////////////////////////////
// Generates a plot of this interpolation function
//////////////////////////////////////////////////////////////////////
void InterpolationMitchellNetravali::plot( FloatArray &fPlot,
                                           FloatArray &dfPlot,
                                           REAL div ) const
{
    int ndivs = (int)ceil( 2.0 * _supportLength / div );

    fPlot.resize( ndivs );
    dfPlot.resize( ndivs );

    for ( int i = 0; i < ndivs; i++ )
    {
        fPlot[i] = evaluate( -1.0 * _supportLength + div * (REAL)i, 0.0 );
        dfPlot[i] = derivative( -1.0 * _supportLength + div * (REAL)i, 0.0 );
    }
}

//////////////////////////////////////////////////////////////////////
// InterpolationBSplineCubic implementation
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
    InterpolationBSplineCubic::InterpolationBSplineCubic( REAL h )
: InterpolationFunction( h, 2, 2.0 * h )
{
}

//////////////////////////////////////////////////////////////////////
// t is the evaluation time, tSample is the location in
// time of the sample with which this function is associated.
//////////////////////////////////////////////////////////////////////
REAL InterpolationBSplineCubic::evaluate( REAL t, REAL tSample ) const
{
    REAL arg;
    REAL functionValue = 0.0;

    arg = abs( t - tSample ) / _h;

    if ( arg >= 2.0 )
    {
        return 0.0;
    }

    if ( arg < 1.0 )
    {
        arg = 1.0 - arg;

        functionValue = -3.0 * arg;
        functionValue = arg * ( 3.0 + functionValue );
        functionValue = arg * ( 3.0 + functionValue );
        functionValue = 1.0 + functionValue;
    }
    else
    {
        arg = 2.0 - arg;

        functionValue = arg * arg * arg;
    }

    functionValue /= 6.0;

    return functionValue;
}

//////////////////////////////////////////////////////////////////////
// Same convention as above - evaluates the derivative
//////////////////////////////////////////////////////////////////////
REAL InterpolationBSplineCubic::derivative( REAL t, REAL tSample ) const
{
    REAL arg;
    REAL derivativeValue = 0.0;
    REAL argSign;

    arg = abs( t - tSample ) / _h;

    argSign = ( t - tSample ) > 0.0 ? 1.0 : -1.0;

    if ( arg >= 2.0 )
    {
        return 0.0;
    }

    if ( arg < 1.0 )
    {
        arg = 1.0 - arg;

        derivativeValue = -9.0 * arg;
        derivativeValue = arg * ( 6.0 + derivativeValue );
        derivativeValue = 3.0 + derivativeValue;
    }
    else
    {
        arg = 2.0 - arg;

        derivativeValue = 3.0 * arg * arg;
    }

    derivativeValue *= -1.0 * argSign;
    derivativeValue /= _h * 6.0;

    return derivativeValue;
}

//////////////////////////////////////////////////////////////////////
// Generates a plot of this interpolation function
//////////////////////////////////////////////////////////////////////
void InterpolationBSplineCubic::plot( FloatArray &fPlot,
        FloatArray &dfPlot,
        REAL div ) const
{
    int ndivs = (int)ceil( 2.0 * _supportLength / div );

    fPlot.resize( ndivs );
    dfPlot.resize( ndivs );

    for ( int i = 0; i < ndivs; i++ )
    {
        fPlot[i] = evaluate( -1.0 * _supportLength + div * (REAL)i, 0.0 );
        dfPlot[i] = derivative( -1.0 * _supportLength + div * (REAL)i, 0.0 );
    }
}

//////////////////////////////////////////////////////////////////////
// InterpolationCatmullRomCubic implementation
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
InterpolationCatmullRomCubic::InterpolationCatmullRomCubic( REAL h )
    : InterpolationFunction( h, 2, 2.0 * h )
{
}

//////////////////////////////////////////////////////////////////////
// t is the evaluation time, tSample is the location in
// time of the sample with which this function is associated.
//////////////////////////////////////////////////////////////////////
REAL InterpolationCatmullRomCubic::evaluate( REAL t, REAL tSample ) const
{
    REAL arg;
    REAL functionValue = 0.0;

    arg = abs( t - tSample ) / _h;

    if ( arg >= 2.0 )
    {
        return 0.0;
    }

    if ( arg < 1.0 )
    {
        arg = 1.0 - arg;

        functionValue = -3.0 * arg;
        functionValue = arg * ( 4.0 + functionValue );
        functionValue = arg * ( 1.0 + functionValue );
    }
    else
    {
        arg = 2.0 - arg;

        functionValue = arg - 1;
        functionValue = arg * arg * functionValue;
    }

    functionValue /= 2.0;

    return functionValue;
}

//////////////////////////////////////////////////////////////////////
// Same convention as above - evaluates the derivative
//////////////////////////////////////////////////////////////////////
REAL InterpolationCatmullRomCubic::derivative( REAL t, REAL tSample ) const
{
    REAL arg;
    REAL derivativeValue = 0.0;
    REAL argSign;

    arg = abs( t - tSample ) / _h;

    argSign = ( t - tSample ) > 0.0 ? 1.0 : -1.0;

    if ( arg >= 2.0 )
    {
        return 0.0;
    }

    if ( arg < 1.0 )
    {
        arg = 1.0 - arg;

        derivativeValue = -9.0 * arg;
        derivativeValue = arg * ( 8.0 + derivativeValue );
        derivativeValue = 1.0 + derivativeValue;
    }
    else
    {
        arg = 2.0 - arg;

        derivativeValue = 3.0 * arg;
        derivativeValue = arg * ( derivativeValue - 2.0 );
    }

    derivativeValue *= -1.0 * argSign;
    derivativeValue /= _h * 2.0;

    return derivativeValue;
}

//////////////////////////////////////////////////////////////////////
// Generates a plot of this interpolation function
//////////////////////////////////////////////////////////////////////
void InterpolationCatmullRomCubic::plot( FloatArray &fPlot,
                                         FloatArray &dfPlot,
                                         REAL div ) const
{
    int ndivs = (int)ceil( 2.0 * _supportLength / div );

    fPlot.resize( ndivs );
    dfPlot.resize( ndivs );

    for ( int i = 0; i < ndivs; i++ )
    {
        fPlot[i] = evaluate( -1.0 * _supportLength + div * (REAL)i, 0.0 );
        dfPlot[i] = derivative( -1.0 * _supportLength + div * (REAL)i, 0.0 );
    }
}

//////////////////////////////////////////////////////////////////////
// InterpolationImprovedPerlin implementation
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
InterpolationImprovedPerlin::InterpolationImprovedPerlin( REAL h )
    : InterpolationFunction( h, 2, 2.0 * h )
{
    findNormalization();
}

//////////////////////////////////////////////////////////////////////
// t is the evaluation time, tSample is the location in
// time of the sample with which this function is associated.
//////////////////////////////////////////////////////////////////////
REAL InterpolationImprovedPerlin::evaluate( REAL t, REAL tSample ) const
{
    REAL arg;
    REAL functionValue = 0.0;

#if 0
    arg = ( _supportLength - abs( t - tSample ) ) / _supportLength;
#endif
    arg = abs( t - tSample ) / _h;
    arg = 2.0 - arg;

    if ( abs( t - tSample ) >= _supportLength  )
    {
        return 0.0;
    }

#if 0
    functionValue = 6.0 * arg - 15.0;
    functionValue = 10.0 + arg * functionValue;
    functionValue *= arg * arg * arg;

    functionValue *= _normalizationCoefficient;
#endif

    functionValue = 3.0 * arg - 15.0;
    functionValue = 20.0 + arg * functionValue;
    functionValue *= arg * arg * arg;

    functionValue /= 32.0;

    return functionValue;
}

//////////////////////////////////////////////////////////////////////
// Same convention as above - evaluates the derivative
//////////////////////////////////////////////////////////////////////
REAL InterpolationImprovedPerlin::derivative( REAL t, REAL tSample ) const
{
    REAL arg;
    REAL derivativeValue = 0.0;
    REAL argSign;

#if 0
    arg = ( _supportLength - abs( t - tSample ) ) / _supportLength;
#endif
    arg = abs( t - tSample ) / _h;
    arg = 2.0 - arg;

    argSign = ( t - tSample ) > 0.0 ? 1.0 : -1.0;

    if ( abs( t - tSample ) >= _supportLength  )
    {
        return 0.0;
    }

#if 0
    derivativeValue = 30.0 * arg - 60.0;
    derivativeValue = 30.0 + arg * derivativeValue;
    derivativeValue *= arg * arg;

    derivativeValue *= -1.0 * argSign;
    derivativeValue /= _supportLength;

    derivativeValue *= _normalizationCoefficient;
#endif

    derivativeValue = 15.0 * arg - 60.0;
    derivativeValue = 60.0 + arg * derivativeValue;
    derivativeValue *= arg * arg;

    derivativeValue *= -1.0 * argSign;
    derivativeValue /= _h * 32.0;

    return derivativeValue;
}

//////////////////////////////////////////////////////////////////////
// Generates a plot of this interpolation function
//////////////////////////////////////////////////////////////////////
void InterpolationImprovedPerlin::plot( FloatArray &fPlot,
                                        FloatArray &dfPlot,
                                        REAL div ) const
{
    int ndivs = (int)ceil( 2.0 * _supportLength / div );

    fPlot.resize( ndivs );
    dfPlot.resize( ndivs );

    for ( int i = 0; i < ndivs; i++ )
    {
        fPlot[i] = evaluate( -1.0 * _supportLength + div * (REAL)i, 0.0 );
        dfPlot[i] = derivative( -1.0 * _supportLength + div * (REAL)i, 0.0 );
    }
}
