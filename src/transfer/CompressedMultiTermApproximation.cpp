//////////////////////////////////////////////////////////////////////
// CompressedMultiTermApproximation.cpp: Implementation of the
//                                       CompressedMultiTermApproximation
//                                       class
//
//////////////////////////////////////////////////////////////////////

#include "CompressedMultiTermApproximation.h"

#include <math/InterpolationFunction.h>

#include <linearalgebra/VECTOR.h>

#include <utils/MathUtil.h>
#include <utils/utils_IO.h>
#include <utils/STLUtil.h>

static const REAL ROOT_2_INV = 1.0 / sqrt( 2 );

//////////////////////////////////////////////////////////////////////
// WaveletManager implementation
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
WaveletManager::WaveletManager( int signalLength, const WaveletBasis &basis )
    : _waveletBasis( basis )
{
    updateSignalLength( signalLength );
}

//////////////////////////////////////////////////////////////////////
// Functions for performing wavelet transforms on provided data.
//////////////////////////////////////////////////////////////////////
void WaveletManager::waveletTransform( REAL *data, int effectiveLength )
{
    // Check to make sure this thread's basis is constructed and initialized
    if ( !_threadBasis.get() ) {
        // Just copy the "master" basis
        _threadBasis.reset( new WaveletBasis( _waveletBasis ) );
        _threadBasis->initBasis( _waveletSignalLength );
    }

    _threadBasis->transform( data, effectiveLength );
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void WaveletManager::waveletTransform( FloatArray &signal, int effectiveLength )
{
    TRACE_ASSERT( signal.size() >= effectiveLength );

    waveletTransform( signal.data(), effectiveLength );
}

//////////////////////////////////////////////////////////////////////
// Matrix version - applies transformation to all rows of the given matrix
//////////////////////////////////////////////////////////////////////
void WaveletManager::waveletTransform( MATRIX &data, int effectiveLength )
{
    TRACE_ASSERT( data.cols() >= effectiveLength );

    for ( int i = 0; i < data.rows(); i++ ) {
        waveletTransform( data.row(i), effectiveLength );
    }
}

//////////////////////////////////////////////////////////////////////
// Inverse wavelet transform
//////////////////////////////////////////////////////////////////////
void WaveletManager::inverseWaveletTransform( REAL *data, int effectiveLength )
{
    // Check to make sure this thread's basis is constructed and initialized
    if ( !_threadBasis.get() ) {
        // Just copy the "master" basis
        _threadBasis.reset( new WaveletBasis( _waveletBasis ) );
        _threadBasis->initBasis( _waveletSignalLength );
    }

    _threadBasis->inverseTransform( data, effectiveLength );
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void WaveletManager::inverseWaveletTransform( FloatArray &signal, int effectiveLength )
{
    TRACE_ASSERT( signal.size() >= effectiveLength );

    inverseWaveletTransform( signal.data(), effectiveLength );
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void WaveletManager::inverseWaveletTransform( MATRIX &data, int effectiveLength )
{
    TRACE_ASSERT( data.cols() >= effectiveLength );

    for ( int i = 0; i < data.rows(); i++ ) {
        inverseWaveletTransform( data.row(i), effectiveLength );
    }
}

//////////////////////////////////////////////////////////////////////
// Gets the transform storage workspace for this thread,
// allocating it if necessary
//////////////////////////////////////////////////////////////////////
FloatArray &WaveletManager::getTransformStorage()
{
    // Check to make sure that this thread's storage copy is constructed
    // and initialized
    if ( !_transformStorage.get() ) {
        _transformStorage.reset( new FloatArray( _waveletSignalLength ) );
    }

    FloatArray                &storage = *_transformStorage;

    if ( storage.size() != _waveletSignalLength ) {
        storage.resize( _waveletSignalLength );
    }

    return storage;
}

//////////////////////////////////////////////////////////////////////
// CompactCompressedMultiTermApproximation implementation
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
CompactCompressedMultiTermApproximation::CompactCompressedMultiTermApproximation(
                                    std::vector<FieldTerm> *fieldTerms[ NUM_ACCEL_DIRECTIONS ],
                                    int baseSample, int signalLength, int startSample,
                                    const Point3d &fieldCenter, REAL fieldRadius,
                                    int fieldResolution,
                                    REAL fieldTimeScale, REAL fieldSupport,
                                    int fieldSampleRate,
                                    REAL c,
                                    bool useSymmetry,
                                    bool matchSamplingRate )
    : CompactRadialApproximation( fieldCenter, fieldRadius, fieldResolution,
                                  fieldTimeScale, fieldSupport,
                                  fieldSampleRate, signalLength, 
                                  useSymmetry, matchSamplingRate, c ),
      // For now, use Daubechies wavelets of order 10, since these
      // seem to work well.  We can expose this as a parameter later if
      // it is helpful.
      WaveletManager( signalLength,
                      WaveletBasis( WaveletBasis::DAUBECHIES_BASIS, 10 ) ),
      _baseSample( baseSample ),
#if 0
      _signalLength( signalLength ),
#endif
      _startSample( startSample ),
#if 0
      _useSymmetry( useSymmetry ),
      _matchSamplingRate( matchSamplingRate ),
#endif
      _useDirectionSet( false ),
      _directionMesh( NULL ),
      _directionTree( NULL )
{
    for ( int direction = 0; direction < NUM_ACCEL_DIRECTIONS; direction++ )
    {
        _fields[ direction ] = fieldTerms[ direction ];

        cout << "Copying term direction " << direction << endl;

        if ( direction > 0 )
        {
            TRACE_ASSERT( _fields[ direction ]->size()
                    == _fields[ direction - 1 ]->size(),
                    "Field term number size mismatch" );

            TRACE_ASSERT( _fields[ direction ]->at( 0 ).size()
                    == _fields[ direction - 1 ]->at( 0 ).size(),
                    "Field term direction resolution mismatch" );
        }

        TRACE_ASSERT( _fields[ direction ]->at( 0 )[ 0 ].size() > 0,
                "Invalid term length" );
    }

    _numTerms = _fields[ 0 ]->size();

    _outputBuffer.reserve( signalLength );

    //REAL                       sampleDiv = 1.0 / (REAL)_fieldSampleRate;
    REAL                       sampleDiv = _sampleDiv;

    _interp = new InterpolationMitchellNetravali( sampleDiv );
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
CompactCompressedMultiTermApproximation::~CompactCompressedMultiTermApproximation()
{
    delete _directionMesh;
    delete _directionTree;
}

#if 0
//////////////////////////////////////////////////////////////////////
// Assuming the object undergoes a half-sine acceleration pulse
// of the given length, add the contribution of this pulse to
// the given listening position
//////////////////////////////////////////////////////////////////////
bool CompactCompressedMultiTermApproximation::addSinePulse(
        Point3d listeningPosition,
        REAL pulseStartTime,
        REAL pulseLength,
        REAL pulseAmplitude[ NUM_ACCEL_DIRECTIONS ],
        SampledFunction &outputSignal,
        REAL scale )
{
    if ( _useDirectionSet ) {
        return addSinePulse_directionMesh( listeningPosition,
                pulseStartTime,
                pulseLength,
                pulseAmplitude,
                outputSignal,
                scale );
    } else {
        return addSinePulse_uniformAngle( listeningPosition,
                pulseStartTime,
                pulseLength,
                pulseAmplitude,
                outputSignal,
                scale );
    }
}
#endif

//////////////////////////////////////////////////////////////////////
// Adds the precomputed acceleration noise signal in the given
// direction to the output signal, with the given time delay
//////////////////////////////////////////////////////////////////////
void CompactCompressedMultiTermApproximation::addPANSignal(
                                                const InterpolationDirection &direction,
                                                const int directions[],
                                                const SignalProperties &signalProps,
                                                REAL delay,
                                                REAL sampleScale, REAL spatialScale,
                                                REAL pulseAmplitude[NUM_ACCEL_DIRECTIONS],
                                                const InterpolationFunction &interp,
                                                SampledFunction &outputSignal )
{
    REAL                       rPower;
    REAL                       totalScaling;
    int                        waveletSignalLength;
    REAL                       transScale, rotScale;

    transScale = spatialScale;
    rotScale = spatialScale * spatialScale;

    // New wavelet signal length based on whether or not we are
    // resampling
    waveletSignalLength = getWaveletSignalLength() / signalProps._rateScale;

    sampleScale *= signalProps._rateModificationScale;

    // Clear the transform storage buffer so that we can start
    // accumulating wavelet coefficients
    clearTransformStorage();

    // Add the contribution from each direction
    for ( int direction_idx = 0; direction_idx < 4; direction_idx++ )
    {
        REAL                   interpCoef = interpolationCoefficient(
                direction, direction_idx );

        // Add the wavelet coefficients for each term from each of the
        // direction fields
        for ( int dir = 0; dir < NUM_ACCEL_DIRECTIONS; dir++ )
        {
            rPower = direction._r;

            for ( int term_idx = 0; term_idx < _numTerms; term_idx++ )
            {
                totalScaling = sampleScale * pulseAmplitude[ dir ]
                    * interpCoef / rPower;

                if ( dir < ROT_X ) {
                    totalScaling *= transScale;
                } else {
                    totalScaling *= rotScale;
                }

                addScaledSignal(
                        dir /* field ID */, term_idx /* term ID */,
                        directions[ direction_idx ] /* direction ID */,
                        totalScaling );

                rPower *= direction._r;
            }
        }
    }

    // Transform the signal back in to the time domain
    FloatArray              &transformStorage = getTransformStorage();

    inverseWaveletTransform( transformStorage, waveletSignalLength );

    // Add the transformed signal to the output
    REAL startTime = delay + signalProps._sampleDiv * (REAL)( _baseSample );
    outputSignal.addSignal( &interp, transformStorage.data(),
            (size_t)signalProps._signalLength, startTime );
}

//////////////////////////////////////////////////////////////////////
// Set this approximation to use a specific direction set,
// rather than a uniform discretization of angle space
//////////////////////////////////////////////////////////////////////
void CompactCompressedMultiTermApproximation::setDirectionSet( const string &directionSetFile )
{
#if 0
    Vector3Array               vertices;
#endif
    vector<Point3d>            vertices;
    vector<Tuple3i>            triangles;
#if 0
    IntArray                   flattenedTriangles;
#endif

    char                       buf[ 1024 ];

    // Read the vertex and triangle sets from disk
    sprintf( buf, "%s_vertices.vector", directionSetFile.c_str() );
    readVector( buf, vertices );

    sprintf( buf, "%s_triangles.vector", directionSetFile.c_str() );
    readVector( buf, triangles );

#if 0
    flattenedTriangles.reserve( triangles.size() * 3 );
    for ( int tri_idx = 0; tri_idx < triangles.size(); tri_idx++ ) {
        flattenedTriangles.push_back( triangles[ tri_idx ][ 0 ] );
        flattenedTriangles.push_back( triangles[ tri_idx ][ 1 ] );
        flattenedTriangles.push_back( triangles[ tri_idx ][ 2 ] );
    }

    _directionMesh = new TriMesh( vertices, flattenedTriangles, 1.0 );
#endif

    _directionMesh = new TriangleMesh<REAL>();

    for ( int vert_idx = 0; vert_idx < vertices.size(); vert_idx++ ) {
        _directionMesh->add_vertex( vertices[ vert_idx ] );
    }

    for ( int tri_idx = 0; tri_idx < triangles.size(); tri_idx++ ) {
        _directionMesh->add_triangle( triangles[ tri_idx ][ 0 ],
                                      triangles[ tri_idx ][ 1 ],
                                      triangles[ tri_idx ][ 2 ] );
    }

    // Initialize the mesh (not sure if this is actually necessary here, but whatever)
    _directionMesh->generate_normals();

    _directionTree = new MeshTree();
    _directionTree->addMesh( _directionMesh );
    _directionTree->buildTree();

    _useDirectionSet = true;
}

//////////////////////////////////////////////////////////////////////
// Adjusts signal properties to (approximately) match the sampling
// rate of outputSignal.
//
// For wavelet signals, this can only be done within a factor of 2
//////////////////////////////////////////////////////////////////////
void CompactCompressedMultiTermApproximation::matchSamplingRate(
                                                        const SampledFunction &outputSignal,
                                                        SignalProperties &signalProps )
{
    // If we want to approximately match the sampling rate of this PAN function
    // to that of the output function, we discard wavelet coefficients until
    // our sampling rate is no more than twice as large as the output sampling
    // rate.
    REAL                     outputSamplingRate;

    outputSamplingRate = (REAL)outputSignal.samplingRate();

    while ( signalProps._samplingRate > 2.0 * outputSamplingRate ) {
        signalProps._samplingRate /= 2.0;
        signalProps._sampleDiv *= 2.0;
        signalProps._rateModificationScale *= ROOT_2_INV;
        signalProps._signalLength /= 2;
        signalProps._rateScale *= 2;
    }
}

//////////////////////////////////////////////////////////////////////
// Adds a sine pulse assuming a uniformly discretized grid in
// angle space
//////////////////////////////////////////////////////////////////////
bool CompactCompressedMultiTermApproximation::addSinePulse_uniformAngle(
                                                        Point3d listeningPosition,
                                                        REAL pulseStartTime,
                                                        REAL pulseLength,
                                                        REAL pulseAmplitude[NUM_ACCEL_DIRECTIONS],
                                                        SampledFunction &outputSignal,
                                                        REAL scale )
{
    REAL                       pulseHalfLength = pulseLength / 2.0;
    REAL                       tPulse;
    REAL                       tDelay;
    REAL                       pulseMiddle = pulseStartTime + pulseHalfLength;
    REAL                       pulseScale = pulseLength / M_PI;
    int                        contributionWidth;

    // Sample division size must be scaled by the spatial scaling
    // applied to the field
    SignalProperties           signalProps( scale * _sampleDiv, _signalLength );

    REAL                       delay;

    REAL                       fieldTimeScale = _fieldTimeScale;
    REAL                       fieldSupport = _fieldSupport;

    REAL                       totalScaling;

    REAL                       sampleScale;
    REAL                       rPower;

    InterpolationDirection     direction;

    // Apply any scaling that might be necessary here due to
    // spatial scaling of the field
    fieldTimeScale *= scale;
    fieldSupport *= scale;

    if ( _matchSamplingRate ) {
#if 0
        // If we want to approximately match the sampling rate of this PAN function
        // to that of the output function, we discard wavelet coefficients until
        // our sampling rate is no more than twice as large as the output sampling
        // rate.
        REAL                     outputSamplingRate;

        outputSamplingRate = (REAL)outputSignal.samplingRate();

        while ( signalProps._samplingRate > 2.0 * outputSamplingRate ) {
            signalProps._samplingRate /= 2.0;
            signalProps._sampleDiv *= 2.0;
            signalProps._rateModificationScale *= ROOT_2_INV;
            signalProps._signalLength /= 2;
            signalProps._rateScale *= 2;
        }
#endif
        matchSamplingRate( outputSignal, signalProps );
    }

    InterpolationMitchellNetravali interp( signalProps._sampleDiv );

    if ( _useSymmetry ) {
        SymmetricScalingFactor( listeningPosition, pulseAmplitude );

        // Rescale according to symmetries
        direction = octantInterpolationDirection( listeningPosition );
    } else {
        direction = interpolationDirection( listeningPosition );
    }

    // Delay based on our proximity to the receiver
    delay = direction._r / _c;

    // Apply spatial scaling
    direction._r /= scale;

    // Figure out how many contributions we will add
    contributionWidth = (int)ceil( pulseHalfLength / fieldTimeScale );

    for ( int pulse_idx = -contributionWidth; pulse_idx <= contributionWidth;
            pulse_idx++ )
    {
        tPulse = pulseMiddle + fieldTimeScale * (REAL)pulse_idx;

        sampleScale = sin( ( tPulse - pulseStartTime ) / pulseScale );
        sampleScale = max( sampleScale, 0.0 );

        // Our pulses are centered in time at _fieldSupport, so take this
        // in to account when computing the delay
        tDelay = tPulse - fieldSupport + delay;

        int                      directions[] = { direction._psi_idx_11,
            direction._psi_idx_21,
            direction._psi_idx_12,
            direction._psi_idx_22 };

        addPANSignal( direction, directions, signalProps, tDelay,
                sampleScale, scale,
                pulseAmplitude, interp, outputSignal );
#if 0
        tPulse = pulseMiddle + fieldTimeScale * (REAL)pulse_idx;

        sampleScale = sin( ( tPulse - pulseStartTime ) / pulseScale );
        sampleScale = max( sampleScale, 0.0 );

        sampleScale *= rateModificationScale;

        // Our pulses are centered in time at _fieldSupport, so take this
        // in to account when computing the delay
        tDelay = tPulse - fieldSupport;

        // Clear the transform storage buffer so that we can start
        // accumulating wavelet coefficients
        clearTransformStorage();

        int                      directions[] = { direction._psi_idx_11,
            direction._psi_idx_21,
            direction._psi_idx_12,
            direction._psi_idx_22 };

        // Add the contribution from each direction
        for ( int direction_idx = 0; direction_idx < 4; direction_idx++ )
        {
            REAL                   interpCoef = interpolationCoefficient(
                    direction, direction_idx );

            // Add the wavelet coefficients for each term from each of the
            // direction fields
            for ( int dir = 0; dir < NUM_ACCEL_DIRECTIONS; dir++ )
            {
                rPower = direction._r;

                for ( int term_idx = 0; term_idx < _numTerms; term_idx++ )
                {
                    totalScaling = sampleScale * pulseAmplitude[ dir ]
                        * interpCoef / rPower;

                    if ( dir < ROT_X ) {
                        totalScaling *= transScale;
                    } else {
                        totalScaling *= rotScale;
                    }

                    addScaledSignal(
                            dir /* field ID */, term_idx /* term ID */,
                            directions[ direction_idx ] /* direction ID */,
                            totalScaling );

                    rPower *= direction._r;
                }
            }
        }

        // Transform the signal back in to the time domain
        FloatArray              &transformStorage = getTransformStorage();

        inverseWaveletTransform( transformStorage, waveletSignalLength );

        // Add the transformed signal to the output
        REAL startTime = tDelay + delay + sampleDiv * (REAL)( _baseSample );
        outputSignal.addSignal( &interp, transformStorage.data(),
                (size_t)signalLength, startTime );
#endif
    }

    return true;
}

//////////////////////////////////////////////////////////////////////
// Adds a sine pulse assuming that angular directions are specified
// by _directionMesh
//////////////////////////////////////////////////////////////////////
bool CompactCompressedMultiTermApproximation::addSinePulse_directionMesh(
                                                    Point3d listeningPosition,
                                                    REAL pulseStartTime,
                                                    REAL pulseLength,
                                                    REAL pulseAmplitude[NUM_ACCEL_DIRECTIONS],
                                                    SampledFunction &outputSignal,
                                                    REAL scale )
{
    REAL                       pulseHalfLength = pulseLength / 2.0;
    REAL                       t;
    REAL                       tPulse;
    REAL                       tDelay;
    REAL                       pulseMiddle = pulseStartTime + pulseHalfLength;
    REAL                       pulseScale = pulseLength / M_PI;
    int                        contributionWidth;

    REAL                       delay;
    REAL                       sampleDiv = _sampleDiv;

    REAL                       fieldTimeScale = _fieldTimeScale;
    REAL                       fieldSupport = _fieldSupport;

    REAL                       transScale = scale;
    REAL                       rotScale = scale * scale;

    REAL                       totalScaling;

    REAL                       sampleValue;
    REAL                       sampleScale;
    REAL                       rPower;
    REAL                       R;

    int                        waveletSignalLength = getWaveletSignalLength();
    int                        signalLength = _signalLength;
    REAL                       samplingRate;
    REAL                       rateModificationScale = 1.0;

    Tuple3i                    directionIndices;
    Point3d                    directionCoefficients;

    Vector3d                   listeningDirection;

    // Apply any scaling that might be necessary here
    sampleDiv *= scale;
    fieldTimeScale *= scale;
    fieldSupport *= scale;

    if ( _matchSamplingRate ) {
        // If we want to approximately match the sampling rate of this PAN function
        // to that of the output function, we discard wavelet coefficients until
        // our sampling rate is no more than twice as large as the output sampling
        // rate.
        REAL                     outputSamplingRate;

        outputSamplingRate = (REAL)outputSignal.samplingRate();
        samplingRate = 1.0 / sampleDiv;

        while ( samplingRate > 2.0 * outputSamplingRate ) {
            samplingRate /= 2.0;
            sampleDiv *= 2.0;
            rateModificationScale *= ROOT_2_INV;
            waveletSignalLength /= 2;
            signalLength /= 2;
        }
    }

    InterpolationMitchellNetravali interp( sampleDiv );

    listeningPosition -= _fieldCenter;

    if ( _useSymmetry ) {
        SymmetricScalingFactor( listeningPosition, pulseAmplitude );
    }

    listeningDirection = listeningPosition;

#if 0
    R = norm( listeningPosition );
#endif
    R = listeningDirection.length();

    findInterpolationPoints( listeningDirection, directionIndices, directionCoefficients );

    // Time delay based on our distance from the receiver
    delay = R / _c;

    // Apply spatial scaling
    R /= scale;

    // Figure out how many contributions we will add
    contributionWidth = (int)ceil( pulseHalfLength / fieldTimeScale );

    for ( int pulse_idx = -contributionWidth; pulse_idx <= contributionWidth;
            pulse_idx++ )
    {
        tPulse = pulseMiddle + fieldTimeScale * (REAL)pulse_idx;

        sampleScale = sin( ( tPulse - pulseStartTime ) / pulseScale );
        sampleScale = max( sampleScale, 0.0 );

        sampleScale *= rateModificationScale;

        // Our pulses are centered in time at _fieldSupport, so take this
        // in to account when computing the delay
        tDelay = tPulse - fieldSupport;

        // Clear the transform storage buffer so that we can start
        // accumulating wavelet coefficients
        clearTransformStorage();

        // Add the contribution from each direction
        for ( int direction_idx = 0; direction_idx < 3; direction_idx++ )
        {
            REAL                   interpCoef = directionCoefficients[ direction_idx ];

            // Add the wavelet coefficients for each term from each of the
            // direction fields
            for ( int dir = 0; dir < NUM_ACCEL_DIRECTIONS; dir++ )
            {
                rPower = R;

                for ( int term_idx = 0; term_idx < _numTerms; term_idx++ )
                {
                    totalScaling = sampleScale * pulseAmplitude[ dir ]
                        * interpCoef / rPower;

                    if ( dir < ROT_X ) {
                        totalScaling *= transScale;
                    } else {
                        totalScaling *= rotScale;
                    }

                    addScaledSignal(
                            dir /* field ID */, term_idx /* term ID */,
                            directionIndices[ direction_idx ] /* direction ID */,
                            totalScaling );

                    rPower *= R;
                }
            }
        }

        // Transform the signal back in to the time domain
        FloatArray              &transformStorage = getTransformStorage();

        inverseWaveletTransform( transformStorage, waveletSignalLength );

        // Add the transformed signal to the output
        REAL startTime = tDelay + delay + sampleDiv * (REAL)( _baseSample );
        outputSignal.addSignal( &interp, transformStorage.data(),
                (size_t)signalLength, startTime );
    }

    return true;
}

//////////////////////////////////////////////////////////////////////
// Samples a direction from our direction mesh
//////////////////////////////////////////////////////////////////////
void CompactCompressedMultiTermApproximation::findInterpolationPoints( Vector3d direction,
                                                                       Tuple3i &pointIndices,
                                                                       Point3d &coefficients )
{
    MeshTree::RayHit           hitData;
    int                        tri_idx;

    direction.normalize();

    // Shoot a ray from the origin
    hitData = _directionTree->intersectRay( Point3d( 0.0, 0.0, 0.0 ), direction );

    TRACE_ASSERT( hitData._hitReference.first == 0,
            "No mesh intersection found" );

    tri_idx = hitData._hitReference.second;
    coefficients = hitData._barycentricPosition;

    pointIndices[ 0 ] = _directionMesh->triangles()[ tri_idx ][ 0 ];
    pointIndices[ 1 ] = _directionMesh->triangles()[ tri_idx ][ 1 ];
    pointIndices[ 2 ] = _directionMesh->triangles()[ tri_idx ][ 2 ];
}

//////////////////////////////////////////////////////////////////////
// CompressedMultiTermApproximation implementation
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Builds a wavelet approximation directly from field data by converting
// the data at each field point in to a wavelet representation, then
// fitting in wavelet space.
//////////////////////////////////////////////////////////////////////
CompressedMultiTermApproximation::CompressedMultiTermApproximation(
                                                        const FloatArray &radii, int numTerms,
                                                        const std::string &shellDataFilePrefix,
                                                        const Point3d &fieldCenter,
                                                        int fieldResolution,
                                                        REAL fieldTimeScale,
                                                        REAL fieldSupport,
                                                        int fieldSampleRate, int field_id,
                                                        REAL tolerance,
                                                        REAL c )
    : RadialApproximation( fieldCenter, 0.0, fieldResolution,
                           fieldTimeScale, fieldSupport, fieldSampleRate,
                           field_id, c ),
      // For now, use Daubechies wavelets of order 10, since these
      // seem to work well.  We can expose this as a parameter later if
      // it is helpful.
      WaveletManager( 0 /* empty signal */, WaveletBasis( WaveletBasis::DAUBECHIES_BASIS, 10 ) ),
      _baseSample( 0 ),
      _signalLength( 0 ),
      _startSample( 0 )
{
    initFieldTerms( shellDataFilePrefix, radii, numTerms, field_id, tolerance );
}

//////////////////////////////////////////////////////////////////////
// Construct by reading from a file
//////////////////////////////////////////////////////////////////////
CompressedMultiTermApproximation::CompressedMultiTermApproximation(
                                                        const std::string &filename,
                                                        int field_id, bool useSymmetry,
                                                        bool matchSamplingRate, REAL c )
    : RadialApproximation( c ),
    // For now, use Daubechies wavelets of order 10, since these
    // seem to work well.  We can expose this as a parameter later if
    // it is helpful.
    WaveletManager( 0 /* empty signal */,
            WaveletBasis( WaveletBasis::DAUBECHIES_BASIS, 10 ) ),
    _baseSample( 0 ),
    _signalLength( 0 ),
    _startSample( 0 ),
    _field_id( field_id ),
    _useSymmetry( useSymmetry ),
    _matchSamplingRate( matchSamplingRate )
{
    read( filename );
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
CompressedMultiTermApproximation::~CompressedMultiTermApproximation()
{
}

//////////////////////////////////////////////////////////////////////
// Saves the field to disk
//////////////////////////////////////////////////////////////////////
void CompressedMultiTermApproximation::write( const std::string &filename ) const
{
    FILE                       *file;
    int                         size;
    int                         numTerms;
    char                        buf[ 1024 ];

    size_t                      bytes_written;

    sprintf( buf, "%s.dat", filename.c_str() );

    file = fopen( buf, "wb" );
    if ( !file )
    {
        cerr << "ERROR: Could not open " << filename << " for writing" << endl;
        abort();
    }

    int                        waveletSignalLength = getWaveletSignalLength();

    // Write basic info about the field first
    bytes_written = fwrite( (void *)&_fieldCenter, sizeof( Point3d ), 1, file );
    bytes_written = fwrite( (void *)&_fieldRadius, sizeof( REAL ), 1, file );
    bytes_written = fwrite( (void *)&_fieldResolution, sizeof( int ), 1, file );
    bytes_written = fwrite( (void *)&_fieldTimeScale, sizeof( REAL ), 1, file );
    bytes_written = fwrite( (void *)&_fieldSupport, sizeof( REAL ), 1, file );
    bytes_written = fwrite( (void *)&_fieldSampleRate, sizeof( int ), 1, file );
    bytes_written = fwrite( (void *)&_baseSample, sizeof( int ), 1, file );
    bytes_written = fwrite( (void *)&_signalLength, sizeof( int ), 1, file );
    bytes_written = fwrite( (void *)&_startSample, sizeof( int ), 1, file );
    bytes_written = fwrite( (void *)&waveletSignalLength, sizeof( int ), 1,
            file );

    numTerms = _fieldTerms.size();
    bytes_written = fwrite( (void *)&numTerms, sizeof( int ), 1, file );

    // Write out each term
    //
    // We will just store everything in the same file here, since we can
    // no longer store things individually as matrices
    for ( int term_idx = 0; term_idx < _fieldTerms.size(); term_idx++ )
    {
        size = _fieldTerms[ term_idx ].size();
        bytes_written = fwrite( (void *)&size, sizeof( int ), 1, file );

        // Write each direction
        for ( int direction_idx = 0; direction_idx < _fieldTerms[ term_idx ].size();
                direction_idx++ )
        {
            const CoefficientArray &term = _fieldTerms[ term_idx ][ direction_idx ];

            size = term.size();
            bytes_written = fwrite( (void *)&size, sizeof( int ), 1, file );

            // Write the data itself
            bytes_written = fwrite( (void *)term.data(), sizeof( WaveletCoefficient ),
                    size, file );
        }
    }

    fclose( file );
}

//////////////////////////////////////////////////////////////////////
// Assuming the object undergoes a half-sine acceleration pulse
// of the given length, add the contribution of this pulse to
// the given listening position
//////////////////////////////////////////////////////////////////////
bool CompressedMultiTermApproximation::addSinePulse( Point3d listeningPosition,
                                                     REAL pulseStartTime,
                                                     REAL pulseLength, REAL pulseAmplitude,
                                                     SampledFunction &outputSignal,
                                                     REAL scale )
{
    REAL                       pulseHalfLength = pulseLength / 2.0;
    REAL                       t;
    REAL                       tPulse;
    REAL                       tDelay;
    REAL                       pulseMiddle = pulseStartTime + pulseHalfLength;
    REAL                       pulseScale = pulseLength / M_PI;
    int                        contributionWidth;

    REAL                       delay;
    REAL                       sampleDiv = _sampleDiv;

    REAL                       fieldTimeScale = _fieldTimeScale;
    REAL                       fieldSupport = _fieldSupport;

    REAL                       sampleValue;
    REAL                       sampleScale;
    REAL                       rPower;

    InterpolationDirection     direction;

    // Apply any scaling that might be necessary here
    sampleDiv *= scale;
    fieldTimeScale *= scale;
    fieldSupport *= scale;

    InterpolationMitchellNetravali interp( sampleDiv );

    if ( _useSymmetry ) {
        pulseAmplitude *= SymmetricScalingFactor( listeningPosition, _field_id );

        // Rescale according to symmetries
        direction = octantInterpolationDirection( listeningPosition );
    } else {
        direction = interpolationDirection( listeningPosition );
    }

    // Time delay based on our proximity to the receiver
    delay = direction._r / _c;

    // Apply spatial scaling
    direction._r /= scale;

    // Figure out how many contributions we will add
    contributionWidth = (int)ceil( pulseHalfLength / fieldTimeScale );

    for ( int pulse_idx = -contributionWidth; pulse_idx <= contributionWidth;
            pulse_idx++ )
    {
        tPulse = pulseMiddle + fieldTimeScale * (REAL)pulse_idx;

        sampleScale = sin( ( tPulse - pulseStartTime ) / pulseScale );
        sampleScale = max( sampleScale, 0.0 );

        // Our pulses are centered in time at _fieldSupport, so take this
        // in to account when computing the delay
        tDelay = tPulse - fieldSupport;

        // Clear the transform storage buffer so that we can start
        // accumulating wavelet coefficients
        clearTransformStorage();

        int                      directions[] = { direction._psi_idx_11,
            direction._psi_idx_21,
            direction._psi_idx_12,
            direction._psi_idx_22 };

        // Add the contribution from each direction
        for ( int direction_idx = 0; direction_idx < 4; direction_idx++ )
        {
            REAL                   interpCoef = interpolationCoefficient(
                    direction, direction_idx );

            // Add the wavelet coefficients for each term.
            // Scale by the interpolation coefficient, and the radial attenuation
            // factor for each term.
            rPower = direction._r;
            for ( int term_idx = 0; term_idx < _fieldTerms.size(); term_idx++ ) {
                addScaledSignal( term_idx, directions[ direction_idx ],
                        sampleScale * pulseAmplitude * interpCoef / rPower );

                rPower *= direction._r;
            }
        }

        // Transform the signal back in to the time domain
        FloatArray              &transformStorage = getTransformStorage();

        inverseWaveletTransform( transformStorage, getWaveletSignalLength() );

        // Add the transformed signal to the output
        for ( int sample_idx = 0; sample_idx < signalLength(); sample_idx++ ) {
            t = sampleDiv * (REAL)( _baseSample + sample_idx );
            t += tDelay + delay;

            outputSignal.addSample( &interp, t, transformStorage[ sample_idx ] );
        }
    }

    return true;
}

//////////////////////////////////////////////////////////////////////
// Evaluates the field stored by this approximation at the given
// listening position
//////////////////////////////////////////////////////////////////////
void CompressedMultiTermApproximation::evaluatePANField( Point3d listeningPosition,
                                                         SampledFunction &outputSignal )
{
    REAL                       delay;
    REAL                       sampleDiv = _sampleDiv;
    REAL                       sampleValue;
    REAL                       t;
    REAL                       pulseAmplitude = 1.0;

    REAL                       rPower;

    InterpolationDirection     direction;

    InterpolationMitchellNetravali interp( sampleDiv );

    // Get the angles to interpolate between
    if ( _useSymmetry ) {
        pulseAmplitude *= SymmetricScalingFactor( listeningPosition, _field_id );

        // Rescale according to symmetries
        direction = octantInterpolationDirection( listeningPosition );
    } else {
        direction = interpolationDirection( listeningPosition );
    }

    delay = direction._r / _c;

    // Clear the transform storage buffer so that we can start
    // accumulating wavelet coefficients
    clearTransformStorage();

    int                        directions[] = { direction._psi_idx_11,
        direction._psi_idx_21,
        direction._psi_idx_12,
        direction._psi_idx_22 };
    // Add the contribution from each direction
    for ( int direction_idx = 0; direction_idx < 4; direction_idx++ )
    {
        REAL                     interpCoef = interpolationCoefficient(
                direction, direction_idx );

        // Add the wavelet coefficients for each term.
        // Scale by the interpolation coefficient, and the radial attenuation
        // factor for each term.
        rPower = direction._r;
        for ( int term_idx = 0; term_idx < _fieldTerms.size(); term_idx++ ) {
            addScaledSignal( term_idx, directions[ direction_idx ],
                    pulseAmplitude * interpCoef / rPower );

            rPower *= direction._r;
        }
    }

    // Transform the signal back in to the time domain
    FloatArray                &transformStorage = getTransformStorage();

    inverseWaveletTransform( transformStorage, getWaveletSignalLength() );

    // Add to the output signal
    for ( int sample_idx = 0; sample_idx < signalLength(); sample_idx++ ) {
        t = sampleDiv * (REAL)( _baseSample + sample_idx );
        t += delay;

        outputSignal.addSample( &interp, t, transformStorage[ sample_idx ] );
    }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
CompactRadialApproximation *
CompressedMultiTermApproximation::buildCompactApproximation( AccelerationSet &allFields,
                                                             const string *directionSetFile )
{
    vector<FieldTerm>         *fields[ NUM_ACCEL_DIRECTIONS ];

    for ( int dir = 0; dir < NUM_ACCEL_DIRECTIONS; dir++ )
    {
        // Yes, this is a bit ugly
        fields[ dir ]
            = ( vector<FieldTerm> * )allFields._allFields[ dir ]->data();
    }

    CompactCompressedMultiTermApproximation *newField
        = new CompactCompressedMultiTermApproximation( fields,
                _baseSample,
                _signalLength,
                _startSample,
                _fieldCenter,
                _fieldRadius,
                _fieldResolution,
                _fieldTimeScale,
                _fieldSupport,
                _fieldSampleRate,
                _c,
                _useSymmetry,
                _matchSamplingRate );

    if ( directionSetFile != NULL ) {
        newField->setDirectionSet( *directionSetFile );
    }

    return newField;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void CompressedMultiTermApproximation::ReadAccelerationSet( const std::string &filePrefix,
                                                            AccelerationSet &fields,
                                                            bool useSymmetry,
                                                            bool matchSamplingRate,
                                                            const string *directionSetSuffix )
{
    char                       buf[ 1024 ];

    for ( int dir = TRANS_X; dir < NUM_ACCEL_DIRECTIONS; dir++ )
    {
        if ( directionSetSuffix != NULL ) {
            // Append the given direction set suffix to the input file name
            sprintf( buf, "%s%s_%d", filePrefix.c_str(),
                    directionSetSuffix->c_str(), dir );
        } else if ( useSymmetry ) {
            sprintf( buf, "%s_octant_%d", filePrefix.c_str(), dir );
        } else {
            sprintf( buf, "%s_%d", filePrefix.c_str(), dir );
        }

        fields._allFields[ dir ] = new CompressedMultiTermApproximation(
                buf, dir, useSymmetry, matchSamplingRate );
    }
}

//////////////////////////////////////////////////////////////////////
// Reads field data from disk
//////////////////////////////////////////////////////////////////////
size_t CompressedMultiTermApproximation::read( const std::string &filename )
{
    FILE                       *file;
    int                         size;
    int                         numTerms;
    char                        buf[ 1024 ];

    size_t                      bytes_read;

    sprintf( buf, "%s.dat", filename.c_str() );

    printf( "Opening CompressedMultiTermApproximation in %s\n", buf );

    file = fopen( buf, "rb" );
    if ( !file )
    {
        cerr << "ERROR: Could not open " << filename << " for reading" << endl;
        abort();
    }

    int                        waveletSignalLength;

    // Read basic info about the field first
    bytes_read = fread( (void *)&_fieldCenter, sizeof( Point3d ), 1, file );
    bytes_read = fread( (void *)&_fieldRadius, sizeof( REAL ), 1, file );
    bytes_read = fread( (void *)&_fieldResolution, sizeof( int ), 1, file );
    bytes_read = fread( (void *)&_fieldTimeScale, sizeof( REAL ), 1, file );
    bytes_read = fread( (void *)&_fieldSupport, sizeof( REAL ), 1, file );
    bytes_read = fread( (void *)&_fieldSampleRate, sizeof( int ), 1, file );
    bytes_read = fread( (void *)&_baseSample, sizeof( int ), 1, file );
    bytes_read = fread( (void *)&_signalLength, sizeof( int ), 1, file );
    bytes_read = fread( (void *)&_startSample, sizeof( int ), 1, file );

    bytes_read = fread( (void *)&waveletSignalLength, sizeof( int ), 1, file );
    updateSignalLength( _signalLength );

    bytes_read = fread( (void *)&numTerms, sizeof( int ), 1, file );
    _fieldTerms.resize( numTerms );

    // Read in each term
    for ( int term_idx = 0; term_idx < numTerms; term_idx++ )
    {
        bytes_read = fread( (void *)&size, sizeof( int ), 1, file );
        _fieldTerms[ term_idx ].resize( size );

        // Read each direction
        for ( int direction_idx = 0; direction_idx < _fieldTerms[ term_idx ].size();
                direction_idx++ )
        {
            CoefficientArray &term = _fieldTerms[ term_idx ][ direction_idx ];

            bytes_read = fread( (void *)&size, sizeof( int ), 1, file );
            term.resize( size );

            // Read in the data itself
            bytes_read = fread( (void *)term.data(), sizeof( WaveletCoefficient ),
                    size, file );
        }
    }

    _sampleDiv = 1.0 / (REAL)_fieldSampleRate;

    fclose( file );

    return bytes_read; 
}

//////////////////////////////////////////////////////////////////////
// Initializes the content of _fieldTerms by performing least
// squares fits in wavelet space
//////////////////////////////////////////////////////////////////////
void CompressedMultiTermApproximation::initFieldTerms( const std::string &shellDataFilePrefix,
                                                       const FloatArray &radii, int numTerms,
                                                       int field_id, REAL tolerance )
{
    MATRIX                     fieldData;
    int                        numDirections;
    char                       filename[ 1024 ];

    // Read in the data at all field points.
    sprintf( filename, "%s_%d.matrix", shellDataFilePrefix.c_str(), field_id );
    fieldData.read( filename );

    numDirections = fieldData.rows() / radii.size();

    TRACE_ASSERT( fieldData.rows() % radii.size() == 0 );

    _signalLength = fieldData.cols();
    updateSignalLength( _signalLength );

    // Set up our field term vectors
    _fieldTerms.resize( numTerms );
    for ( int term_idx = 0; term_idx < numTerms; term_idx++ ) {
        _fieldTerms[ term_idx ].resize( numDirections );
    }

    findBaseSample( radii );

    for ( int direction_idx = 0; direction_idx < numDirections; direction_idx++ )
    {
        initSingleDirection( radii, numTerms, fieldData, direction_idx, tolerance );
    }
}

//////////////////////////////////////////////////////////////////////
// Initializes a single field term
//////////////////////////////////////////////////////////////////////
void CompressedMultiTermApproximation::initSingleDirection( const FloatArray &radii, int numTerms,
                                                            const MATRIX &fieldData,
                                                            int direction_idx,
                                                            REAL tolerance )
{
    int                        numDirections;

    // Space for storing the signals we want to transform
    MATRIX                     signalWorkspace( radii.size(),
            getWaveletSignalLength() );

    numDirections = fieldData.rows() / radii.size();

    copyDirectionTerms( radii, fieldData, numDirections, direction_idx,
            signalWorkspace );

    // Get wavelet transforms of the field at each radii
    waveletTransform( signalWorkspace, getWaveletSignalLength() );

    // Figure out which wavelet coefficients we want to use to represent
    // this part of the signal
    thresholdDirectionCoefficients( direction_idx, signalWorkspace, tolerance );

    solveDirectionCoefficients( direction_idx, signalWorkspace, radii );
}

//////////////////////////////////////////////////////////////////////
// Copies all terms from one direction, adjusting for offsets due
// to radii differences
//////////////////////////////////////////////////////////////////////
void CompressedMultiTermApproximation::copyDirectionTerms( const FloatArray &radii,
                                                           const MATRIX &fieldData,
                                                           int numDirections, int direction_idx,
                                                           MATRIX &termStorage )
{
    REAL                       radius;
    int                        sampleOffset;
    int                        signalEnd;

    for ( int radius_idx = 0; radius_idx < radii.size(); radius_idx++ ) {
        radius = radii[ radius_idx ];

        sampleOffset = (int)( radius / _c / _sampleDiv );

        sampleOffset = -sampleOffset - _baseSample;

        TRACE_ASSERT( sampleOffset >= 0 );

        signalEnd = getWaveletSignalLength() - sampleOffset;
        signalEnd = min( fieldData.cols(), signalEnd );

        TRACE_ASSERT( sampleOffset + fieldData.cols() < termStorage.cols() );

        // Copy the data
        for ( int i = 0; i < signalEnd; i++ ) {
            termStorage( radius_idx, sampleOffset + i )
                = fieldData( direction_idx + radius_idx * numDirections, i );
        }
    }
}

//////////////////////////////////////////////////////////////////////
// Figures out which coefficients to keep for a single direction
//////////////////////////////////////////////////////////////////////
void CompressedMultiTermApproximation::thresholdDirectionCoefficients( int direction_idx,
                                                                       const MATRIX &termStorage,
                                                                       REAL tolerance )
{
    FloatArray                 thresholds( termStorage.rows() );

    for ( int radius_idx = 0; radius_idx < termStorage.rows(); radius_idx++ ) {
        REAL                     maxValue = 0.0;

        for ( int coef_idx = 0; coef_idx < termStorage.cols(); coef_idx++ ) {
            maxValue = max( maxValue, abs( termStorage( radius_idx, coef_idx ) ) );
        }

        thresholds[ radius_idx ] = maxValue * tolerance;
    }

    // Threshold coefficients based on a tolerance
    for ( int coef_idx = 0; coef_idx < termStorage.cols(); coef_idx++ ) {
        bool                     useCoef = false;

        for ( int radius_idx = 0; radius_idx < termStorage.rows(); radius_idx++ ) {
            if ( abs(termStorage(radius_idx, coef_idx)) > thresholds[radius_idx] ) {
                useCoef = true;
                break;
            }
        }

        if ( useCoef ) {
            for ( size_t term_idx = 0; term_idx < _fieldTerms.size(); term_idx++ ) {
                // Leave the coefficient value empty for now
                _fieldTerms[ term_idx ][ direction_idx ].push_back(
                        WaveletCoefficient( coef_idx, 0.0 ) );
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////
// Solves for the coefficients in one direction
//////////////////////////////////////////////////////////////////////
void CompressedMultiTermApproximation::solveDirectionCoefficients( int direction_idx,
                                                                   const MATRIX &termStorage,
                                                                   const FloatArray &radii )
{
    MATRIX                     systemMatrix, rhs;
    VECTOR                     singularValues;
    int                        rank;

    buildWaveletLSSystem( direction_idx, termStorage, radii, systemMatrix, rhs );

    // Solve and copy coefficients
    MATRIX::LeastSquares_TSVD( systemMatrix, rhs, singularValues, rank, 1e-10 );

    extractSolvedCoefficients( direction_idx, rhs );
}

//////////////////////////////////////////////////////////////////////
// Builds the least squares system for a single direction
//////////////////////////////////////////////////////////////////////
void CompressedMultiTermApproximation::buildWaveletLSSystem( int direction_idx,
                                                             const MATRIX &termStorage,
                                                             const FloatArray &radii,
                                                             MATRIX &systemMatrix, MATRIX &rhs )
{
    int                        numTerms = _fieldTerms.size();
    const CoefficientArray    &directionCoefs = _fieldTerms[ 0 ][ direction_idx ];

    systemMatrix.resizeAndWipe( radii.size(), _fieldTerms.size() );
    rhs.resizeAndWipe( radii.size(), _fieldTerms[ 0 ][ direction_idx ].size() );

    // Build the system matrix
    for ( int radius_idx = 0; radius_idx < radii.size(); radius_idx++ ) {
        REAL                     r = radii[ radius_idx ];

        for ( int term_idx = 0; term_idx < numTerms; term_idx++ ) {
            systemMatrix( radius_idx, term_idx ) = 1.0 / r;

            r *= radii[ radius_idx ];
        }
    }

    // Build the right hand side
    for ( size_t radius_idx = 0; radius_idx < radii.size(); radius_idx++ ) {
        for ( size_t comp_coef_idx = 0; comp_coef_idx < directionCoefs.size();
                comp_coef_idx++ )
        {
            // The actual wavelet coefficient we are storing
            int                    coef_idx = directionCoefs[ comp_coef_idx ].first;

            rhs( radius_idx, comp_coef_idx ) = termStorage( radius_idx, coef_idx );
        }
    }
}

//////////////////////////////////////////////////////////////////////
// Copies solved coefficients for a single direction
//////////////////////////////////////////////////////////////////////
void CompressedMultiTermApproximation::extractSolvedCoefficients( int direction_idx,
                                                                  const MATRIX &rhs )
{
    //int                        numTerms = _fieldTerms.size();

    // Copy terms
    for ( int term_idx = 0; term_idx < _fieldTerms.size(); term_idx++ ) {
        CoefficientArray        &directionCoefs = _fieldTerms[ term_idx ]
            [ direction_idx ];

        TRACE_ASSERT( directionCoefs.size() == rhs.cols() );

        for ( size_t coef_idx = 0; coef_idx < directionCoefs.size(); coef_idx++ ) {
            directionCoefs[ coef_idx ].second = rhs( term_idx, coef_idx );
        }
    }
}

//////////////////////////////////////////////////////////////////////
// Finds the base sample for this field, based on provided shell data
//////////////////////////////////////////////////////////////////////
void CompressedMultiTermApproximation::findBaseSample( const FloatArray &radii )
{
    REAL                       maxRadius = radii.back();

    int                        sampleOffset;

    // Figure out how far back in time the signal at maxRadius could have
    // come from, assuming that it started at the origin
    sampleOffset = (int)ceil( maxRadius / _c / _sampleDiv );

    _baseSample = -sampleOffset;
}
