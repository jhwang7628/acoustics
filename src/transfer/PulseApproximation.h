//////////////////////////////////////////////////////////////////////
// PulseApproximation.h: Interface for the PulseApproximation class
//
//////////////////////////////////////////////////////////////////////

#ifndef PULSE_APPROXIMATION_H
#define PULSE_APPROXIMATION_H

#include <geometry/Point3.hpp>

#include <linearalgebra/MATRIX.h>

#include <math/SampledFunction.h>

#include <TYPES.h>

#include <vector>

#include "RadialApproximation.h"

//////////////////////////////////////////////////////////////////////
// PulseApproximation class
//
// Simple data-driven approximation of far field acceleration noise
// from a single object
//////////////////////////////////////////////////////////////////////
class PulseApproximation : public RadialApproximation {
    public:
        PulseApproximation( const std::vector<std::vector<FloatArray> > &fieldData,
                            const Point3d &fieldCenter, REAL fieldRadius,
                            int fieldResolution,
                            REAL fieldTimeScale, REAL fieldSupport,
                            int fieldSampleRate, int field_id,
                            REAL c = 343.0 );

        // Construct by reading from a file
        PulseApproximation( const std::string &filename, REAL c = 343.0 );

        PulseApproximation() {}

        // Destructor
        virtual ~PulseApproximation();

        // Saves the field to disk
        virtual void write( const std::string &filename ) const;

        // Assuming the object undergoes a half-sine acceleration pulse
        // of the given length, add the contribution of this pulse to
        // the given listening position
        virtual bool addSinePulse( Point3d listeningPosition,
                                   REAL pulseStartTime,
                                   REAL pulseLength, REAL pulseAmplitude,
                                   SampledFunction &outputSignal,
                                   REAL scale = 1.0 );

        // Evaluates the field stored by this approximation at the given
        // listening position
        virtual void evaluatePANField( Point3d listeningPosition,
                                       SampledFunction &outputSignal );

        // Resamples this field on only one octant of angle space
        PulseApproximation *resampleOctant( int octantResolution ) const;

        // Resamples the field to a fixed set of points
        PulseApproximation *resamplePoints( const Vector3Array &newDirections ) const;

        virtual void *data()
        {
            return (void *)&_fieldData;
        }

        virtual int signalLength() const
        {
            return ( _fieldData.size() > 0 ) ? _fieldData[ 0 ].size() : 0;
        }

        virtual CompactRadialApproximation *buildCompactApproximation( 
                                                    AccelerationSet &allFields,
                                                    const std::string *directionSetFile = NULL )
        {
            printf( "Not implemented\n" );
            abort();
            return NULL;
        }

    public:

        static void ReadAccelerationSet( const std::string &filePrefix,
                                         AccelerationSet &fields );

    protected:
        // Reads field data from disk
        virtual size_t read( const std::string &filename );

        // Gets the given direction field at the given sample
        virtual REAL sampleDirection( int direction_idx, int sample_idx,
                int term_idx, int field_idx ) const
        {
            // FIXME: a bit of a hack here - use term_idx as an offset
            return _fieldData[ direction_idx + term_idx ][ sample_idx ];
        }

    private:
        std::vector<FloatArray>  _fieldData;

};

#endif
