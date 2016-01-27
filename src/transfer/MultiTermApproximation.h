//////////////////////////////////////////////////////////////////////
// MultiTermApproximation.h: Interface for the MultiTermApproximation
//                           class
//
//////////////////////////////////////////////////////////////////////

#ifndef MULTI_TERM_APPROXIMATION_H
#define MULTI_TERM_APPROXIMATION_H

#include <geometry/Point3.hpp>

#include <linearalgebra/MATRIX.h>
#include <linearalgebra/SPARSE_MATRIX.h>

#include <math/SampledFunction.h>

#include <utils/trace.h>

#include <TYPES.h>

#include <vector>

#include "RadialApproximation.h"

//////////////////////////////////////////////////////////////////////
// Version of the above class that stores everything in one set
// so that we can save the cost of computing interpolation
// coefficients, etc.
//////////////////////////////////////////////////////////////////////
class CompactMultiTermApproximation : public CompactRadialApproximation
{
    private:
        typedef std::vector<FloatArray>  FieldTerm;

    public:
        CompactMultiTermApproximation( std::vector<FieldTerm> *fieldTerms[ NUM_ACCEL_DIRECTIONS ],
                                       int baseSample, int signalLength, int startSample,
                                       const Point3d &fieldCenter, REAL fieldRadius,
                                       int fieldResolution,
                                       REAL fieldTimeScale, REAL fieldSupport,
                                       int fieldSampleRate,
                                       REAL c = 343.0 );

        virtual ~CompactMultiTermApproximation();

        // Assuming the object undergoes a half-sine acceleration pulse
        // of the given length, add the contribution of this pulse to
        // the given listening position
        virtual bool addSinePulse( Point3d listeningPosition,
                                   REAL pulseStartTime,
                                   REAL pulseLength,
                                   REAL pulseAmplitude[ NUM_ACCEL_DIRECTIONS ],
                                   SampledFunction &outputSignal,
                                   REAL scale = 1.0 );

        // Gets the given direction field at the given sample
        virtual REAL sampleDirection( int direction_idx, int sample_idx,
                int term_idx, int field_idx ) const
        {
            return _fields[ field_idx ]->at( term_idx )[ direction_idx ][ sample_idx ];
        }

    protected:

    private:
        std::vector<FieldTerm>  *_fields[ NUM_ACCEL_DIRECTIONS ];

        int                      _baseSample;
#if 0
        int                      _signalLength;
#endif
        int                      _startSample;

};

//////////////////////////////////////////////////////////////////////
// MultiTermApproximation class
//
// Simple data-driven approximation of far field acceleration noise
// from a single object
//////////////////////////////////////////////////////////////////////
class MultiTermApproximation : public RadialApproximation {
    public:
        // Provide the file name prefix from which we can load
        // the precomputed data at each radial shell
        MultiTermApproximation( const FloatArray &radii, int numTerms,
                                const std::string &shellDataFilePrefix,
                                const Point3d &fieldCenter,
                                int fieldResolution,
                                REAL fieldTimeScale, REAL fieldSupport,
                                int fieldSampleRate, int field_id,
                                REAL svdTolerance,
                                REAL c = 343.0 );

        // Construct by reading from a file
        MultiTermApproximation( const std::string &filename, REAL c = 343.0 );

        MultiTermApproximation()
        {
        }

        // Destructor
        virtual ~MultiTermApproximation();

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

        virtual void *data()
        {
            return (void *)&_fieldTerms;
        }

        virtual int signalLength() const
        {
            return _signalLength;
        }

        virtual CompactRadialApproximation *buildCompactApproximation( 
                                                AccelerationSet &allFields,
                                                const std::string *directionSetFile = NULL );

        // Builds a new multi-term approximation by resampling this one to a new
        // angular resolution
        MultiTermApproximation *resampleApproximation( int newResolution ) const;

    public:

        static void ReadAccelerationSet( const std::string &filePrefix,
                                         AccelerationSet &fields );

        // Builds the least squares system matrix for a multi-term expansion
        static IndexRange BuildMultiTermSystem( const FloatArray &inputRadii,
                                                int numSamples, int numTerms,
                                                int fieldSampleRate,
                                                MATRIX &systemMatrix,
                                                REAL c = 343.0 );

        // Reorders field data so that all radii for a single shell
        // are put together in the same column
        static void BuildMultiTermRHS( const MATRIX &fieldData, MATRIX &rhs, int numRadii );

        // When building a multi-term expansion, adds matrix entries for
        // a single radius value
        static void AddRadiusTerms( REAL radius, int numSamples,
                                    int radiusExponent,
                                    int fieldSampleRate,
                                    SPARSE_MATRIX &systemMatrix,
                                    int base_column,
                                    int start_row,
                                    REAL c = 343.0 );

        static void AddRadiusTerms( REAL radius, int numSamples,
                                    int radiusExponent,
                                    int fieldSampleRate,
                                    MATRIX &systemMatrix,
                                    int base_column,
                                    int start_row, int start_col,
                                    REAL c = 343.0 );

        // Helper function for the above
        static void GetRadiusIndexRange( REAL radius, int numSamples,
                                         int fieldSampleRate,
                                         int &startIndex, int &endIndex,
                                         REAL c = 343.0 );

    protected:
        // Reads field data from disk
        virtual size_t read( const std::string &filename );

        // Gets the given direction field at the given sample
        virtual REAL sampleDirection( int direction_idx, int sample_idx,
                int term_idx, int field_idx ) const
        {
            return _fieldTerms[ term_idx ][ direction_idx ][ sample_idx ];
        }

    private:
        // Initializes the content of _fieldTerms
        //
        // Performs a least squares solve to find the best solution for
        // the given number of terms
        void initFieldTerms( const std::string &shellDataFilePrefix,
                             const FloatArray &radii, int numTerms,
                             int field_id, REAL svdTolerance );

        // Helper function for init field terms.  Once the least squares
        // solve is done, this packs the matrix data from the solution in
        // to vectors in _fieldTerms
        void extractFieldTerms( const MATRIX &solvedData,
                                const MATRIX &systemMatrix,
                                const IndexRange &columnRange,
                                int numTerms );

    private:
        typedef std::vector<FloatArray>  FieldTerm;

        std::vector<FieldTerm>   _fieldTerms;

        // The "base" sample index.  That is, sample zero is assumed to have
        // this (presumably negative) index.  This is because the pulse at the
        // origin may have had to start before time zero to reach the field
        // points by the desired time.
        int                      _baseSample;
        int                      _signalLength;

        // The starting sample (ie. we assume that all samples between
        // zero and this one have zero value)
        int                      _startSample;

};

#endif
