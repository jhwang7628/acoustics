//////////////////////////////////////////////////////////////////////
// SampledFunction.h: Interface for the SampledFunction class
//
//////////////////////////////////////////////////////////////////////

#ifndef SAMPLED_FUNCTION_H
#define SAMPLED_FUNCTION_H

#include <TYPES.h>

#include "InterpolationFunction.h"

//////////////////////////////////////////////////////////////////////
// SampledFunction class
//
// Models a function sampled at some rate, with support for adding
// new samples at non-aligned positions
//////////////////////////////////////////////////////////////////////
class SampledFunction {
    public:
        SampledFunction( int samplingRate );
        SampledFunction();
        SampledFunction( const SampledFunction &f );

        // Destructor
        virtual ~SampledFunction();

        // Adds a sample to the given point in time
        void addSample( REAL t, REAL sampleValue );

        // Add a sample distribued with the given interpolation
        // function
        void addSample( const InterpolationFunction *interp,
                        REAL t, REAL sampleValue );

        // Adds a signal defined with the given interpolation function
        void addSignal( const InterpolationFunction *interp,
                        const REAL *data, size_t dataSize,
                        REAL startTime );

        const FloatArray &data() const { return _function; }

        void clear()
        {
            _function.clear();
            _function.resize( INITIAL_LENGTH, 0.0 );
        }

        SampledFunction &operator=( const SampledFunction &f );

        static REAL EvaluateFunction( const InterpolationFunction *interp,
                                      const REAL *data, size_t dataSize,
                                      REAL startTime, REAL t );

        int samplingRate() const
        {
            return _samplingRate;
        }

    protected:

    private:
        int                      _samplingRate;
        REAL                     _h;

        static const InterpolationFunction::FunctionName  INTERPOLATION_TYPE;
        static const int                                  INITIAL_LENGTH;

        InterpolationFunction   *_interpolationFunction;

        // The signal itself
        FloatArray               _function;

};

#endif
