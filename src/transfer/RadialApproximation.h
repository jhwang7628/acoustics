//////////////////////////////////////////////////////////////////////
// RadialApproximation.h: Interface for the RadialApproximation class
//
//////////////////////////////////////////////////////////////////////

#ifndef RADIAL_APPROXIMATION_H
#define RADIAL_APPROXIMATION_H

#include <geometry/Point3.hpp>

#include <linearalgebra/Vector3.hpp>

#include <math/SampledFunction.h>

#include <TYPES.h>

#include <vector>

//////////////////////////////////////////////////////////////////////
// Models a field defined as a function of angular direction
// surrounding a point in space (alternately, we can view this as
// modelling a function stored on the surface of a sphere)
//
// We use this base class for handling things like radial
// interpolation, etc.
//////////////////////////////////////////////////////////////////////
class RadialField {
    public:
        RadialField( const Point3d &fieldCenter, REAL fieldRadius,
                     int fieldResolution );

        RadialField()
        {
        }

        virtual ~RadialField();

    protected:
        // Use this structure to represent a directional interpolation
        // point in a field discretized on the surface of a sphere.  We
        // assume that the theta and phi angle spaces haves been discretized
        // uniformly (with a constant angular division).
        struct InterpolationDirection {
            // Angle divisions
            REAL                   _thetaDiv;
            REAL                   _phiDiv;

            // Angles of the four directions surrounding the one to be
            // interpolated
            REAL                   _theta_1;
            REAL                   _theta_2;
            REAL                   _phi_1;
            REAL                   _phi_2;

            // Indices for these 4 directions
            int                    _psi_idx_11;
            int                    _psi_idx_12;
            int                    _psi_idx_21;
            int                    _psi_idx_22;

            // Actual angles at which to evaluate the field
            REAL                   _theta;
            REAL                   _phi;
            REAL                   _r;

        };

        // Determines an interpolation direction for the given
        // point in 3D space
        InterpolationDirection interpolationDirection( const Point3d &x ) const;

        // This version of the function returns the interpolation direction
        // assuming that only the positive octant (x,y,z >= 0) is discretized
        InterpolationDirection octantInterpolationDirection( const Point3d &x ) const;

        // Interpolates the field at the given time sample
        //
        // Provide field_idx here in case one of several different fields
        // needs to be sampled
        REAL interpolatedFieldSample( int sample_idx,
                                      const InterpolationDirection &direction,
                                      int term_idx = 0, int field_idx = 0 ) const;

        // Interpolation coefficients for each of the 4 directions
        //
        //   corner_idx:   0 <--> 11
        //                 1 <--> 21
        //                 2 <--> 12
        //                 3 <--> 22
        static REAL interpolationCoefficient( const InterpolationDirection &direction,
                                              int corner_idx );

        // Gets the given direction field at the given sample
        virtual REAL sampleDirection( int direction_idx, int sample_idx,
                                      int term_idx, int field_idx ) const = 0;

        // For fields with symmetric geometry about each axis, we can consider
        // only the values of the field in the first quadrant (ie. x, y, z > 0).
        // This function returns the appropriate scaling factor for a given
        // position and field ID.
        //
        // It also modifies the listening position to put it in the first
        // quadrant
        static REAL SymmetricScalingFactor( Point3d &listeningPosition, int field_id );

        static void SymmetricScalingFactor( Point3d &listeningPosition,
                                            REAL pulseAmplitudes[ NUM_ACCEL_DIRECTIONS ] );

        static REAL SymmetricScalingFactor( unsigned int quadrant, int field_id );

        // Returns the quadrant of this position, and reflects it
        // in to the first quadrant
        static unsigned int PositionQuadrant( Point3d &listeningPosition );

    protected:
        // Look-up table for scaling factors when dealing with symmetric objects
        // (e.g., ellipsoids)
        static REAL              _SymmetricScalingFactors[][NUM_OCTANTS];

        Point3d                  _fieldCenter;
        REAL                     _fieldRadius;
        int                      _fieldResolution;

};

//////////////////////////////////////////////////////////////////////
// Parameter structure for passing around signal data (e.g., sampling
// rate, sample division size, etc.)
//////////////////////////////////////////////////////////////////////
struct SignalProperties {
    SignalProperties( REAL sampleDiv, REAL signalLength )
        : _sampleDiv( sampleDiv ),
          _signalLength( signalLength ),
          // Assume no resampling
          _rateScale( 1 ),
          _rateModificationScale( 1.0 )
    {
        _samplingRate = 1.0 / _sampleDiv;
    }

    void updateSampleDiv( int sampleDiv )
    {
        _sampleDiv = sampleDiv;
        _samplingRate = 1.0 / _sampleDiv;
    }

    //////////////////////////////////////////////////////////////////////
    // Basic signal information
    //////////////////////////////////////////////////////////////////////
    REAL                       _sampleDiv;
    REAL                       _samplingRate;
    int                        _signalLength;

    //////////////////////////////////////////////////////////////////////
    // Properties used when downsampling a signal from its original
    // sampling rate
    //////////////////////////////////////////////////////////////////////

    // Integer by which we are dividing the original sampling rate
    // (e.g., rateScale = 4 means we are sampling at 1/4th the input rate)
    int                        _rateScale;

    // In case resample requires us to re-scale the signal to
    // obtain valid results, store this rescaling factor here
    REAL                       _rateModificationScale;

};

//////////////////////////////////////////////////////////////////////
// Compacted radial approximation that handles a set of multiple
// directions all at once
//////////////////////////////////////////////////////////////////////
class CompactRadialApproximation : public RadialField
{
    public:
        CompactRadialApproximation( const Point3d &fieldCenter,
                                    REAL fieldRadius,
                                    REAL fieldResolution,
                                    REAL fieldTimeScale, REAL fieldSupport,
                                    int fieldSampleRate,
                                    int signalLength,
                                    bool useSymmetry,
                                    bool matchSamplingRate,
                                    REAL c = 343.0 );

        virtual ~CompactRadialApproximation();

        // Assuming the object undergoes a half-sine acceleration pulse
        // of the given length, add the contribution of this pulse to
        // the given listening position
        virtual bool addSinePulse( Point3d listeningPosition,
                                   REAL pulseStartTime,
                                   REAL pulseLength,
                                   REAL pulseAmplitude[ NUM_ACCEL_DIRECTIONS ],
                                   SampledFunction &outputSignal,
                                   REAL scale = 1.0 );

        // Adds the precomputed acceleration noise signal in the given
        // direction to the output signal, with the given time delay
        virtual void addPANSignal( const InterpolationDirection &direction,
                                   const int directions[],
                                   const SignalProperties &signalProps,
                                   REAL delay,
                                   // sampleScale: the amount by which we
                                   //     with to scale the signal added to
                                   ///    outputSignal
                                   // spatialScale: the amount by which we
                                   //     are spatially scaling the noise source
                                   REAL sampleScale, REAL spatialScale,
                                   REAL pulseAmplitude[NUM_ACCEL_DIRECTIONS],
                                   const InterpolationFunction &interp,
                                   SampledFunction &outputSignal ) {}

        REAL h() const
        {
            return _interp->h();
        }

    private:
        // Adjusts signal properties to (approximately) match the sampling
        // rate of outputSignal.
        //
        // For wavelet signals, this can only be done within a factor of 2
        virtual void matchSamplingRate( const SampledFunction &outputSignal,
                                        SignalProperties &signalProps ) {}

    protected:
        REAL                     _fieldTimeScale;
        REAL                     _fieldSupport;
        int                      _fieldSampleRate;

        REAL                     _sampleDiv;
        int                      _signalLength;

        // If this is true, we only store one octant of the field and
        // use the field's symmetry to evaluate the others
        bool                     _useSymmetry;

        // Whether or not to approximately match the field sampling rate
        // to that of output signals
        bool                     _matchSamplingRate;

        REAL                     _c;

        InterpolationFunction   *_interp;

};

//////////////////////////////////////////////////////////////////////
// RadialApproximation class
//
// Models data-driven approximation of far field acceleration noise
// from a single object
//////////////////////////////////////////////////////////////////////
class RadialApproximation : public RadialField {
    public:
        RadialApproximation( const Point3d &fieldCenter, REAL fieldRadius,
                             int fieldResolution,
                             REAL fieldTimeScale, REAL fieldSupport,
                             int fieldSampleRate, int field_id,
                             REAL c = 343.0 );

        // Construct by reading from a file
        RadialApproximation( REAL c = 343.0 );

        // Destructor
        virtual ~RadialApproximation();

        // Saves the field to disk
        virtual void write( const std::string &filename ) const = 0;

        // Assuming the object undergoes a half-sine acceleration pulse
        // of the given length, add the contribution of this pulse to
        // the given listening position
        virtual bool addSinePulse( Point3d listeningPosition,
                                   REAL pulseStartTime,
                                   REAL pulseLength, REAL pulseAmplitude,
                                   SampledFunction &outputSignal,
                                   REAL scale = 1.0 ) = 0;

        // Evaluates the field stored by this approximation at the given
        // listening position
        virtual void evaluatePANField( Point3d listeningPosition,
                                       SampledFunction &outputSignal ) = 0;

        // Returns the raw data associated with this field
        virtual void *data() = 0;

        virtual int signalLength() const = 0;

        struct AccelerationSet {
            RadialApproximation    *_allFields[ NUM_ACCEL_DIRECTIONS ];

            AccelerationSet()
            {
                for ( int accel_dir = 0; accel_dir < NUM_ACCEL_DIRECTIONS; accel_dir++ )
                {
                    _allFields[ accel_dir ] = NULL;
                }
            }

            ~AccelerationSet()
            {
            }

            void clear()
            {
                for ( int accel_dir = 0; accel_dir < NUM_ACCEL_DIRECTIONS; accel_dir++ )
                {
                    if ( _allFields[ accel_dir ] )
                    {
                        delete _allFields[ accel_dir ];
                    }
                }
            }
        };

        virtual CompactRadialApproximation *buildCompactApproximation( 
                                            AccelerationSet &allFields,
                                            const std::string *directionSetFile = NULL ) = 0;

    public:

        static void WriteAccelerationSet( const std::string &filePrefix,
                                          const AccelerationSet &fields );

        // Adds pulses for all acceleration sets, given a single listening
        // position and time scale (but multiple amplitudes)
        static bool AddSinePulse( AccelerationSet &fields,
                                  const Point3d &listeningPosition,
                                  REAL pulseStartTime, REAL pulseLength,
                                  REAL pulseAmplitude[ NUM_ACCEL_DIRECTIONS ],
                                  SampledFunction &outputSignal );

    protected:
        // Reads field data from disk
        virtual size_t read( const std::string &filename ) = 0;

    protected:
        REAL                     _fieldTimeScale;
        REAL                     _fieldSupport;
        int                      _fieldSampleRate;
        REAL                     _sampleDiv;

        REAL                     _c;

        InterpolationFunction   *_interp;

};

#endif
