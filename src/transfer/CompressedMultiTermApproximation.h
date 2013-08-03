//////////////////////////////////////////////////////////////////////
// CompressedMultiTermApproximation.h: Interface for the
//                                     MultiTermApproximation class
//
//////////////////////////////////////////////////////////////////////

#ifndef _COMPRESSED_MULTI_TERM_APPROXIMATION_H
#define _COMPRESSED_MULTI_TERM_APPROXIMATION_H

#include <geometry/Point3.hpp>
#include <geometry/MeshTree.h>
#include <geometry/TriangleMesh.hpp>

#include <linearalgebra/MATRIX.h>

#include <math/SampledFunction.h>

#include <utils/trace.h>

#include <TYPES.h>

#include <vector>

// Wavelet includes
#include <gsl/gsl_wavelet.h>

#include <boost/thread/tss.hpp>

typedef gsl_wavelet            Wavelet;
typedef gsl_wavelet_workspace  WaveletWorkspace;

class CompressedMultiTermApproximation;

#include "RadialApproximation.h"

//////////////////////////////////////////////////////////////////////
// A simple wrapper for the GSL Wavelet library
//////////////////////////////////////////////////////////////////////
class WaveletBasis {
    public:
        enum WaveletType {
            DAUBECHIES_BASIS = 0,
            HAAR_BASIS,
            BSPLINE_BASIS
        };

        // Provide the basis type and order here.  Daubechies wavelets of
        // order 10 seem to be a good choice, so we will use this as
        // the default.
        WaveletBasis(WaveletType type = DAUBECHIES_BASIS, int order = 10)
            : _type(type),
            _order(10),
            _wavelet(NULL),
            _workspace(NULL)
    {
    }

        WaveletBasis(const WaveletBasis &basis)
            : _type(basis._type) ,
            _order(basis._order),
            _wavelet(NULL),
            _workspace(NULL)
    {
    }

        ~WaveletBasis()
        {
            clear();
        }

        WaveletBasis &operator=(const WaveletBasis &basis)
        {
            this->_type = basis._type;
            this->_order = basis._order;

            // This wavelet is now uninitialized
            clear();
        }

        void clear() {
            if ( _wavelet ) {
                gsl_wavelet_free( _wavelet );
                _wavelet = NULL;
            }

            if ( _workspace ) {
                gsl_wavelet_workspace_free( _workspace );
                _workspace = NULL;
            }
        }

        // Initialize this wavelet to perform transforms on signals
        // of the given length
        void initBasis(int waveletSignalLength)
        {
            clear();

            switch (_type) {
                case DAUBECHIES_BASIS:
                default:
                    _wavelet = gsl_wavelet_alloc( gsl_wavelet_daubechies, _order );
                    break;
                case HAAR_BASIS:
                    _wavelet = gsl_wavelet_alloc( gsl_wavelet_haar, _order );
                    break;
                case BSPLINE_BASIS:
                    _wavelet = gsl_wavelet_alloc( gsl_wavelet_bspline, _order );
                    break;
            }

            TRACE_ASSERT( _wavelet != NULL, "GSL Wavelet construction failed" );

            _workspace = gsl_wavelet_workspace_alloc( waveletSignalLength );
        }

        // Forward and inverse transforms on the provided signal.
        // Note that this overrides the signal with wavelet basis coefficients.
        // Use the effectiveLength argument if we only wish to consider the
        // contents of signal up to a given length
        void transform(REAL *signal, int effectiveLength)
        {
            int                    status;

            status = gsl_wavelet_transform_forward( _wavelet, signal, 1,
                    effectiveLength, _workspace );

            TRACE_ASSERT( status == GSL_EINVAL, "GSL wavelet transform failed" );
        }

        void inverseTransform(REAL *signal, int effectiveLength)
        {
            int                    status;

            status = gsl_wavelet_transform_inverse( _wavelet, signal, 1,
                    effectiveLength, _workspace );

            TRACE_ASSERT( status == GSL_SUCCESS,
                    "GSL inverse wavelet trasnform failed" );
        }

    private:
        WaveletType              _type;
        int                      _order;
        Wavelet                 *_wavelet;
        WaveletWorkspace        *_workspace;

};

//////////////////////////////////////////////////////////////////////
// WaveletManager class
//
// Since we will have both a CompressedMultiTermApproximation and
// and compacted version of this class, we will collect
// all wavelet-related functionality in to a single parent class.
// This WaveletManager class provides functions for performing foward
// and inverse wavelet transforms, given some wavelet basis
//////////////////////////////////////////////////////////////////////

class WaveletManager {
    public:
        WaveletManager( int signalLength, const WaveletBasis &basis );

        // Performs a wavelet transform using our wavelet data on the given
        // data array
        void waveletTransform( REAL *data, int effectiveLength );
        void waveletTransform( FloatArray &signal, int effectiveLength );

        // Matrix version - applies transformation to all rows of the given matrix
        void waveletTransform( MATRIX &data, int effectiveLength );

        // Inverse wavelet transform
        void inverseWaveletTransform( REAL *data, int effectiveLength );
        void inverseWaveletTransform( FloatArray &signal, int effectiveLength );
        void inverseWaveletTransform( MATRIX &data, int effectiveLength );

        // Gets the transform storage workspace for this thread,
        // allocating it if necessary
        FloatArray &getTransformStorage();

        static int WaveletSignalLength( int signalLength )
        {
            // Wavelet signal length has to be a power of 2 - we can pad signals
            // shorter than this with zeros
            int waveletSignalLength = 1;
            while ( waveletSignalLength < signalLength ) {
                waveletSignalLength *= 2;
            }

            return waveletSignalLength;
        }

        void updateSignalLength( int signalLength )
        {
            _waveletSignalLength = WaveletSignalLength( signalLength );
        }

        int getWaveletSignalLength() const
        {
            return _waveletSignalLength;
        }

        void clearTransformStorage()
        {
            FloatArray            &transformStorage = getTransformStorage();

            for ( int i = 0; i < transformStorage.size(); i++ ) {
                transformStorage[ i ] = 0.0;
            }
        }

    private:
        //////////////////////////////////////////////////////////////////////
        // Wavelet data; signal length, and a buffer we can use for constructing
        // and performing wavelet transforms
        //////////////////////////////////////////////////////////////////////

        int                                    _waveletSignalLength;

        // Thread local copies of storage for accumulating wavelet
        // coefficients.
        boost::thread_specific_ptr<FloatArray>    _transformStorage;

        // Since the WaveletBasis class has a workspace that is modified 
        // during wavelet transform application, we need thread-local copies
        // so that this can be evaluated in parallel.
        WaveletBasis                              _waveletBasis;
        boost::thread_specific_ptr<WaveletBasis>  _threadBasis;


};


//////////////////////////////////////////////////////////////////////
// CompressedMultiTermApproximation class
//
// Approximates a radial field with a series of radially-decaying
// terms of the form (1/r, 1/r^2, ...).  The signals associated with
// each of these terms are represented using the provided wavelet
// compression basis.  Compression is achieved by thresholding the
// coefficients in this basis.
//////////////////////////////////////////////////////////////////////
class CompressedMultiTermApproximation : public RadialApproximation,
    private WaveletManager
{
    public:
        // Builds a wavelet approximation directly from field data by converting
        // the data at each field point in to a wavelet representation, then
        // fitting in wavelet space.
        CompressedMultiTermApproximation( const FloatArray &radii, int numTerms,
                                          const std::string &shellDataFilePrefix,
                                          const Point3d &fieldCenter,
                                          int fieldResolution,
                                          REAL fieldTimeScale,
                                          REAL fieldSupport,
                                          int fieldSampleRate, int field_id,
                                          REAL tolerance,
                                          REAL c = 343.0 );

        // Construct by reading from a file
        CompressedMultiTermApproximation( const std::string &filename,
                                          int field_id,
                                          bool useSymmetry,
                                          bool matchSamplingRate,
                                          REAL c = 343.0 );

        CompressedMultiTermApproximation()
            // For now, use Daubechies wavelets of order 10, since these
            // seem to work well.  We can expose this as a parameter later if
            // it is helpful.
            : WaveletManager( 0, WaveletBasis( WaveletBasis::DAUBECHIES_BASIS, 10 ) )
        {
        }

        // Destructor
        virtual ~CompressedMultiTermApproximation();

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

    public:

        static void ReadAccelerationSet( const std::string &filePrefix,
                                         AccelerationSet &fields,
                                         bool useSymmetry,
                                         bool matchSamplingRate,
                                         const std::string *directionSetSuffix );

    protected:
        // Reads field data from disk
        virtual void read( const std::string &filename );

        // Gets the given direction field at the given sample
        virtual REAL sampleDirection( int direction_idx, int sample_idx,
                                      int term_idx, int field_idx ) const
        {
            TRACE_ASSERT( NULL, "Should never get here" );
            return 0.0;
        }

    private:
        // Initializes the content of _fieldTerms by performing least
        // squares fits in wavelet space
        void initFieldTerms( const std::string &shellDataFilePrefix,
                             const FloatArray &radii, int numTerms,
                             int field_id, REAL tolerance );

        // Initializes a single field term
        void initSingleDirection( const FloatArray &radii, int numTerms,
                                  const MATRIX &fieldData, int direction_idx,
                                  REAL tolerance );

        // Copies all terms from one direction, adjusting for offsets due
        // to radii differences
        void copyDirectionTerms( const FloatArray &radii,
                                 const MATRIX &fieldData,
                                 int numDirections, int direction_idx,
                                 MATRIX &termStorage );

        // Figures out which coefficients to keep for a single direction
        void thresholdDirectionCoefficients( int direction_idx,
                                             const MATRIX &termStorage,
                                             REAL tolerance );

        // Solves for the coefficients in one direction
        void solveDirectionCoefficients( int direction_idx,
                                         const MATRIX &termStorage,
                                         const FloatArray &radii );

        // Builds the least squares system for a single direction
        void buildWaveletLSSystem( int direction_idx,
                                   const MATRIX &termStorage,
                                   const FloatArray &radii,
                                   MATRIX &systemMatrix, MATRIX &rhs );

        // Copies solved coefficients for a single direction
        void extractSolvedCoefficients( int direction_idx, const MATRIX &rhs );

        // Finds the base sample for this field, based on provided shell data
        void findBaseSample( const FloatArray &radii );

        // Adds one of our sparse field terms to the transform workspace,
        // with the given scaling
        inline void addScaledSignal( int term_idx, int direction_idx, REAL scale )
        {
            const CoefficientArray  &signal = _fieldTerms[ term_idx ][ direction_idx ];
            FloatArray              &transformStorage = getTransformStorage();

            for ( int coef_idx = 0; coef_idx < signal.size(); coef_idx++ ) {
                const WaveletCoefficient  &coef = signal[ coef_idx ];

                transformStorage[ coef.first ] += scale * coef.second;
            }
        }

    private:
        typedef std::pair<int, float>            WaveletCoefficient;
        typedef std::vector<WaveletCoefficient>  CoefficientArray;
        typedef std::vector<CoefficientArray>    FieldTerm;

        std::vector<FieldTerm>                 _fieldTerms;

        // The "base" sample index.  That is, sample zero is assumed to have
        // this (presumably negative) index.  This is because the pulse at the
        // origin may have had to start before time zero to reach the field
        // points by the desired time.
        int                                    _baseSample;
        int                                    _signalLength;

        // The starting sample (ie. we assume that all samples between
        // zero and this one have zero value)
        int                                    _startSample;

        int                                    _field_id;

        // If this is true, we only store one octant of the field and
        // use the field's symmetry to evaluate the others
        bool                                   _useSymmetry;
        bool                                   _matchSamplingRate;

};

//////////////////////////////////////////////////////////////////////
// Version of the above class that stores everything in one set
// so that we can save the cost of computing interpolation
// coefficients, etc.
//////////////////////////////////////////////////////////////////////
class CompactCompressedMultiTermApproximation : public CompactRadialApproximation,
    private WaveletManager
{
    private:
        typedef std::pair<int, float>            WaveletCoefficient;
        typedef std::vector<WaveletCoefficient>  CoefficientArray;
        typedef std::vector<CoefficientArray>    FieldTerm;

    public:
        CompactCompressedMultiTermApproximation(
                                    std::vector<FieldTerm> *fieldTerms[ NUM_ACCEL_DIRECTIONS ],
                                    int baseSample, int signalLength, int startSample,
                                    const Point3d &fieldCenter, REAL fieldRadius,
                                    int fieldResolution,
                                    REAL fieldTimeScale, REAL fieldSupport,
                                    int fieldSampleRate,
                                    REAL c = 343.0,
                                    bool useSymmetry = false,
                                    bool matchSamplingRate = false );

        virtual ~CompactCompressedMultiTermApproximation();

#if 0
        // Assuming the object undergoes a half-sine acceleration pulse
        // of the given length, add the contribution of this pulse to
        // the given listening position
        virtual bool addSinePulse( Point3d listeningPosition,
                REAL pulseStartTime,
                REAL pulseLength,
                REAL pulseAmplitude[ NUM_ACCEL_DIRECTIONS ],
                SampledFunction &outputSignal,
                REAL scale = 1.0 );
#endif

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
                                   SampledFunction &outputSignal );

        // Set this approximation to use a specific direction set,
        // rather than a uniform discretization of angle space
        void setDirectionSet( const std::string &directionSetFile );

        // Gets the given direction field at the given sample
        virtual REAL sampleDirection( int direction_idx, int sample_idx,
                                      int term_idx, int field_idx ) const
        {
            TRACE_ASSERT( NULL, "Should never get here" );
            return 0.0;
        }

    protected:

    private:
        // Adjusts signal properties to (approximately) match the sampling
        // rate of outputSignal.
        //
        // For wavelet signals, this can only be done within a factor of 2
        virtual void matchSamplingRate( const SampledFunction &outputSignal,
                                        SignalProperties &signalProps );

        // Adds a sine pulse assuming a uniformly discretized grid in
        // angle space
        bool addSinePulse_uniformAngle( Point3d listeningPosition,
                                        REAL pulseStartTime,
                                        REAL pulseLength,
                                        REAL pulseAmplitude[NUM_ACCEL_DIRECTIONS],
                                        SampledFunction &outputSignal,
                                        REAL scale );

        // Adds a sine pulse assuming that angular directions are specified
        // by _directionMesh
        bool addSinePulse_directionMesh( Point3d listeningPosition,
                                         REAL pulseStartTime,
                                         REAL pulseLength,
                                         REAL pulseAmplitude[NUM_ACCEL_DIRECTIONS],
                                         SampledFunction &outputSignal,
                                         REAL scale );

        // Samples a direction from our direction mesh
        void findInterpolationPoints( Vector3d direction,
                                      Tuple3i &pointIndices,
                                      Point3d &coefficients );

        // Adds one of our sparse field terms to the transform workspace,
        // with the given scaling
        inline void addScaledSignal( int field_idx, int term_idx,
                                     int direction_idx, REAL scale )
        {
            const CoefficientArray  &signal
                = _fields[ field_idx ]->at( term_idx )[ direction_idx ];

            FloatArray              &transformStorage = getTransformStorage();

            for ( int coef_idx = 0; coef_idx < signal.size(); coef_idx++ ) {
                const WaveletCoefficient  &coef = signal[ coef_idx ];

                transformStorage[ coef.first ] += scale * coef.second;
            }
        }

    private:
        std::vector<FieldTerm>  *_fields[ NUM_ACCEL_DIRECTIONS ];

        int                      _baseSample;
#if 0
        int                      _signalLength;
#endif
        int                      _startSample;

        int                      _numTerms;

        FloatArray               _outputBuffer;

#if 0
        // If this is true, we only store one octant of the field and
        // use the field's symmetry to evaluate the others
        bool                                   _useSymmetry;

        // Whether or not to approximately match the field sampling rate
        // to that of output signals
        bool                                   _matchSamplingRate;
#endif

        // Whether or not interpret the field as being represented at
        // a specific set of directions (rather than a regular grid in
        // angle space)
        bool                                   _useDirectionSet;

        TriangleMesh<REAL>                    *_directionMesh;
        MeshTree                              *_directionTree;

};

#endif // _COMPRESSED_MULTI_TERM_APPROXIMATION_H
