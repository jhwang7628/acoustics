//////////////////////////////////////////////////////////////////////
// InterpolationFunction.h: Interface for the InterpolationFunction
//                          class
//
//////////////////////////////////////////////////////////////////////

#ifndef INTERPOLATION_FUNCTION_H
#define INTERPOLATION_FUNCTION_H

#include <linearalgebra/MATRIX.h>
#include <linearalgebra/VECTOR.h>

#include <TYPES.h>

#include "Function.h"

#include <string>

//////////////////////////////////////////////////////////////////////
// InterpolationFunction class
//
// Interpolation functions used for evaluating and differentiating
// a discrete time signal between time samples
//////////////////////////////////////////////////////////////////////
class InterpolationFunction {
    public:
        InterpolationFunction( REAL h, int support, REAL supportLength );

        // Destructor
        virtual ~InterpolationFunction();

        // t is the evaluation time, tSample is the location in
        // time of the sample with which this function is associated.
        virtual REAL             evaluate( REAL t, REAL tSample ) const = 0;

        // Same convention as above - evaluates the derivative
        virtual REAL             derivative( REAL t, REAL tSample ) const = 0;

        // Generates a plot of this interpolation function
        virtual void             plot( FloatArray &fPlot, FloatArray &dfPlot,
                                       REAL div ) const = 0;

        // Accessors
        REAL                     h() const
        {
            return _h;
        }

        int                      support() const
        {
            return _support;
        }

        REAL                     supportLength() const
        {
            return _supportLength;
        }

        virtual const char      *name() const = 0;

        void setTimeScale( REAL h )
        {
            _h = h;
            _supportLength = _h * (REAL)_support;
        }

        enum FunctionName {
            INTERPOLATION_POLY6 = 0,
            INTERPOLATION_MITCHELL_NETRAVALI,
            PERLIN_IMPROVED,
            BSPLINE_CUBIC,
            CATMULL_ROM_CUBIC
        };

        // Looks up the enumerator function name
        static FunctionName GetFunctionName( std::string functionName );

        // Returns the requested interpolation function
        static InterpolationFunction *GetFunction( FunctionName name,
                                                   REAL h,
                                                   REAL supportLength = -1.0 );

        // Returns a Function object for the specified interpolation function
        static Function GetFunctionContainer( FunctionName name,
                                              REAL h,
                                              REAL supportLength = -1.0 );

    protected:
        // Finds a normalization coefficient for this function (so that
        // evaluations at distances of _h sum up to 1).
        void findNormalization();

    protected:
        // Time length between samples
        REAL                     _h;

        // Support of the interpolation function, expressed
        // in number of samples
        int                      _support;

        // Actual real-valued support length for the function
        REAL                     _supportLength;

        // Normalization coefficient
        REAL                     _normalizationCoefficient;

};

//////////////////////////////////////////////////////////////////////
// 6th degree polynomial interpolation function
//
// Let r = ( t - t_k )
// Let h = support length
// f(t) = ( 1 - (r / h)^2 )^3
//////////////////////////////////////////////////////////////////////
class InterpolationPoly6 : public InterpolationFunction
{
    public:
        // supportLength is expressed as a fraction of h
        InterpolationPoly6( REAL h, REAL supportLength );

        // Destructor
        virtual ~InterpolationPoly6() {}

        // t is the evaluation time, tSample is the location in
        // time of the sample with which this function is associated.
        virtual REAL             evaluate( REAL t, REAL tSample ) const;

        // Same convention as above - evaluates the derivative
        virtual REAL             derivative( REAL t, REAL tSample ) const;

        // Generates a plot of this interpolation function
        virtual void             plot( FloatArray &fPlot, FloatArray &dfPlot,
                                       REAL div ) const;

        virtual const char      *name() const
        {
            return "Degree 6 Polynomial";
        }

    protected:

};

//////////////////////////////////////////////////////////////////////
// Mitchell-Netravali Cubic
//
//                  / -15 (1-|t|)^3 + 18 (1-|t|)^2 + (1-|t|) + 2  (if |t| < 1)
// f(t) = (1/18) * |  5 (2-|t|)^3 - 3 (2-|t|)^2             (if 1 <= |t| <= 2)
//                  \ 0                                              otherwise
//////////////////////////////////////////////////////////////////////
class InterpolationMitchellNetravali : public InterpolationFunction
{
    public:
        // supportLength is expressed as a fraction of h
        InterpolationMitchellNetravali( REAL h );

        // Destructor
        virtual ~InterpolationMitchellNetravali() {}

        // t is the evaluation time, tSample is the location in
        // time of the sample with which this function is associated.
        virtual REAL             evaluate( REAL t, REAL tSample ) const;

        // Same convention as above - evaluates the derivative
        virtual REAL             derivative( REAL t, REAL tSample ) const;

        // Generates a plot of this interpolation function
        virtual void             plot( FloatArray &fPlot, FloatArray &dfPlot,
                                       REAL div ) const;

        virtual const char      *name() const
        {
            return "Mitchell-Netravali";
        }

    protected:

};

//////////////////////////////////////////////////////////////////////
// B-spline Cubic
//
//                 / -3 (1-|t|)^3 + 3 (1-|t|)^2 + 3 (1-|t|) + 1  (if |t| < 1)
// f(t) = (1/6) * |  (2-|t|)^3                             (if 1 <= |t| <= 2)
//                 \ 0                                              otherwise
//////////////////////////////////////////////////////////////////////
class InterpolationBSplineCubic : public InterpolationFunction
{
    public:
        // supportLength is expressed as a fraction of h
        InterpolationBSplineCubic( REAL h );

        // Destructor
        virtual ~InterpolationBSplineCubic() {}

        // t is the evaluation time, tSample is the location in
        // time of the sample with which this function is associated.
        virtual REAL             evaluate( REAL t, REAL tSample ) const;

        // Same convention as above - evaluates the derivative
        virtual REAL             derivative( REAL t, REAL tSample ) const;

        // Generates a plot of this interpolation function
        virtual void             plot( FloatArray &fPlot, FloatArray &dfPlot,
                                       REAL div ) const;

        virtual const char      *name() const
        {
            return "B-spline Cubic";
        }

    protected:

};

//////////////////////////////////////////////////////////////////////
// Catmull-Rom Cubic
//
//                 / -3 (1-|t|)^3 + 4 (1-|t|)^2 + (1-|t|)  (if |t| < 1)
// f(t) = (1/2) * |  (2-|t|)^3 - (2-|t|)^2           (if 1 <= |t| <= 2)
//                 \ 0                                        otherwise
//////////////////////////////////////////////////////////////////////
class InterpolationCatmullRomCubic : public InterpolationFunction
{
    public:
        // supportLength is expressed as a fraction of h
        InterpolationCatmullRomCubic( REAL h );

        // Destructor
        virtual ~InterpolationCatmullRomCubic() {}

        // t is the evaluation time, tSample is the location in
        // time of the sample with which this function is associated.
        virtual REAL             evaluate( REAL t, REAL tSample ) const;

        // Same convention as above - evaluates the derivative
        virtual REAL             derivative( REAL t, REAL tSample ) const;

        // Generates a plot of this interpolation function
        virtual void             plot( FloatArray &fPlot, FloatArray &dfPlot,
                                       REAL div ) const;

        virtual const char      *name() const
        {
            return "Catmull-Rom Cubic";
        }

    protected:

};

//////////////////////////////////////////////////////////////////////
// 5th degree polynomial from Improved Perlin Noise (2002 paper)
//
// f(t) = (1/32) * ( 3 (2-|t|)^5 - 15 (2-|t|)^4 + 20 (2-|t|)^3
//////////////////////////////////////////////////////////////////////
class InterpolationImprovedPerlin : public InterpolationFunction
{
    public:
        // supportLength is expressed as a fraction of h
        InterpolationImprovedPerlin( REAL h );

        // Destructor
        virtual ~InterpolationImprovedPerlin() {}

        // t is the evaluation time, tSample is the location in
        // time of the sample with which this function is associated.
        virtual REAL             evaluate( REAL t, REAL tSample ) const;

        // Same convention as above - evaluates the derivative
        virtual REAL             derivative( REAL t, REAL tSample ) const;

        // Generates a plot of this interpolation function
        virtual void             plot( FloatArray &fPlot, FloatArray &dfPlot,
                                       REAL div ) const;

        virtual const char      *name() const
        {
            return "Improved Perlin";
        }

    protected:

};

#endif
