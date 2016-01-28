//////////////////////////////////////////////////////////////////////
// PulseApproximation.cpp: Implementation of the PulseApproximation
//                         class
//
//////////////////////////////////////////////////////////////////////

#include "PulseApproximation.h"

#include <math/InterpolationFunction.h>

#include <utils/MathUtil.h>
#include <utils/IO.h>

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
PulseApproximation::PulseApproximation( const vector<vector<FloatArray> > &fieldData,
                                        const Point3d &fieldCenter,
                                        REAL fieldRadius,
                                        int fieldResolution,
                                        REAL fieldTimeScale,
                                        REAL fieldSupport,
                                        int fieldSampleRate,
                                        int field_id,
                                        REAL c )
    : RadialApproximation( fieldCenter, fieldRadius, fieldResolution,
                           fieldTimeScale, fieldSupport, fieldSampleRate,
                           field_id, c )
{
    _fieldData.resize( fieldData.size() );

    // Copy field field_id from fieldData
    for ( size_t field_position = 0; field_position < fieldData.size();
            field_position++ )
    {
        _fieldData[ field_position ] = fieldData[ field_position ][ field_id ];
    }
}

//////////////////////////////////////////////////////////////////////
// Construct by reading from a file
//////////////////////////////////////////////////////////////////////
PulseApproximation::PulseApproximation( const string &filename, REAL c )
    : RadialApproximation( c )
{
    read( filename );
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
PulseApproximation::~PulseApproximation()
{
}

//////////////////////////////////////////////////////////////////////
// Saves the field to disk
//////////////////////////////////////////////////////////////////////
void PulseApproximation::write( const string &filename ) const
{
    FILE                      *file;
    //int                        size;

    //size_t                     bytes_written;

    char                       buf[ 1024 ];

    MATRIX                     dataMatrix( _fieldData.size(), signalLength() );

    sprintf( buf, "%s.dat", filename.c_str() );

    file = fopen( buf, "wb" );
    if ( !file )
    {
        cerr << "ERROR: Could not open " << filename << " for writing" << endl;
        abort();
    }

    // Copy field data to the matrix
    for ( size_t field_point_idx = 0; field_point_idx < _fieldData.size();
            field_point_idx++ )
        for ( int sample_idx = 0; sample_idx < signalLength(); sample_idx++ )
        {
            dataMatrix( field_point_idx, sample_idx )
                = _fieldData[ field_point_idx ][ sample_idx ];
        }

    sprintf( buf, "%s.matrix", filename.c_str() );
    dataMatrix.write( buf );

#if 0
    // Write the field data
    size = _fieldData.size();

    bytes_written = fwrite( (void *)&size, sizeof( int ), 1, file );

    for ( int field_point_idx = 0; field_point_idx < _fieldData.size();
            field_point_idx++ )
    {
        size = _fieldData[ field_point_idx ].size();

        // Write field point size
        bytes_written = fwrite( (void *)&size, sizeof( int ), 1, file );

        // Write field point data
        bytes_written = fwrite( (void *)_fieldData[ field_point_idx ].data(),
                sizeof( REAL ), size, file );
    }
#endif

    // Write the remaining info needed to reconstruct the field
    fwrite( (void *)&_fieldCenter, sizeof( Point3d ), 1, file );
    fwrite( (void *)&_fieldRadius, sizeof( REAL ), 1, file );
    fwrite( (void *)&_fieldResolution, sizeof( int ), 1, file );
    fwrite( (void *)&_fieldTimeScale, sizeof( REAL ), 1, file );
    fwrite( (void *)&_fieldSupport, sizeof( REAL ), 1, file );
    fwrite( (void *)&_fieldSampleRate, sizeof( int ), 1, file );

    fclose( file );
}

//////////////////////////////////////////////////////////////////////
// Assuming the object undergoes a half-sine acceleration pulse
// of the given length, add the contribution of this pulse to
// the given listening position
//////////////////////////////////////////////////////////////////////
bool PulseApproximation::addSinePulse( Point3d listeningPosition,
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
    REAL                       attenuation;
    REAL                       sampleDiv = 1.0 / (REAL)_fieldSampleRate;

    REAL                       sampleValue;
    REAL                       sampleScale;

    InterpolationDirection     direction;

    if ( !_interp )
    {
        _interp = new InterpolationMitchellNetravali( sampleDiv );
    }

    // Apply any scaling that might be necessary here
    sampleDiv *= scale;
    _interp->setTimeScale( sampleDiv );

    direction = interpolationDirection( listeningPosition );

    delay = ( direction._r - _fieldRadius ) / _c;
    attenuation = _fieldRadius / direction._r;

    // Figure out how many contributions we will add
    contributionWidth = (int)ceil( pulseHalfLength / _fieldTimeScale );

    for ( int pulse_idx = -contributionWidth; pulse_idx <= contributionWidth;
            pulse_idx++ )
    {
        tPulse = pulseMiddle + _fieldTimeScale * (REAL)pulse_idx;

        sampleScale = sin( ( tPulse - pulseStartTime ) / pulseScale );
        sampleScale = max( sampleScale, 0.0 );
#if 0
        cout << SDUMP( sampleScale ) << endl;
#endif
        sampleScale *= attenuation;

        // Our pulses are centered in time at _fieldSupport, so take this
        // in to account when computing the delay
        tDelay = tPulse - _fieldSupport;

        // Add contributions from this pulse to the output signal
        for ( int sample_idx = 0; sample_idx < signalLength(); sample_idx++ )
        {
            t = sampleDiv * (REAL)sample_idx;
            t += tDelay + delay;

            sampleValue = interpolatedFieldSample( sample_idx, direction );

            //outputSignal.addSample( t, sampleValue * sampleScale * pulseAmplitude );
            outputSignal.addSample( _interp, t,
                    sampleValue * sampleScale * pulseAmplitude );
        }
    }

    return true;
}

//////////////////////////////////////////////////////////////////////
// Evaluates the field stored by this approximation at the given
// listening position
//////////////////////////////////////////////////////////////////////
void PulseApproximation::evaluatePANField( Point3d listeningPosition,
                                           SampledFunction &outputSignal )
{
    TRACE_ASSERT( NULL, "Not implemented" );
}

//////////////////////////////////////////////////////////////////////
// Resamples this field on only one octant of angle space
//////////////////////////////////////////////////////////////////////
PulseApproximation *PulseApproximation::resampleOctant( int octantResolution ) const
{
    PulseApproximation        *newField = new PulseApproximation();
    Vector3Array               newDirections;

    int                        numRadii;
    int                        numDirections;
    int                        new_numDirections;

    numDirections = 2 * _fieldResolution * _fieldResolution + 2;

    newField->_fieldCenter = _fieldCenter;
    newField->_fieldRadius = _fieldRadius;
    newField->_fieldResolution = octantResolution;
    newField->_fieldTimeScale = _fieldTimeScale;
    newField->_fieldSupport = _fieldSupport;
    newField->_fieldSampleRate = _fieldSampleRate;
    newField->_c = _c;
    newField->_interp = NULL;

    // Generate the octant direction set
    new_numDirections = MathUtil::GenerateOctantPoints( newField->_fieldCenter,
            1.0, octantResolution,
            newDirections );

    // FIXME: debugging
    for ( size_t i = 0; i < newDirections.size(); i++ ) {
        printf( "Direction %lu: ", i );
        cout << newDirections[ i ] << endl;
    }

    numRadii = _fieldData.size() / numDirections;

    TRACE_ASSERT( _fieldData.size() % numDirections == 0 );

    newField->_fieldData.resize( numRadii * new_numDirections );
    for ( size_t dir_idx = 0; dir_idx < newField->_fieldData.size(); dir_idx++ ) {
        newField->_fieldData[ dir_idx ].resize( signalLength() );
    }

    for ( int dir_idx = 0; dir_idx < new_numDirections; dir_idx++ ) {
        InterpolationDirection   direction;

        direction = interpolationDirection( newDirections[ dir_idx ] );

        for ( int rad_idx = 0; rad_idx < numRadii; rad_idx++ ) {
            int                    dof_idx = dir_idx + rad_idx * new_numDirections;

            for ( int sample_idx = 0; sample_idx < signalLength(); sample_idx++ ) {
                // Use the term_idx argument as an offset
                newField->_fieldData[ dof_idx ][ sample_idx ]
                    = interpolatedFieldSample( sample_idx, direction,
                            rad_idx * numDirections );
            }
        }
    }

    return newField;
}

//////////////////////////////////////////////////////////////////////
// Resamples the field to a fixed set of points
//////////////////////////////////////////////////////////////////////
PulseApproximation *PulseApproximation::resamplePoints( const Vector3Array &newDirections ) const
{
    PulseApproximation        *newField = new PulseApproximation();

    int                        numRadii;
    int                        numDirections;
    int                        new_numDirections;

    numDirections = 2 * _fieldResolution * _fieldResolution + 2;

    newField->_fieldCenter = _fieldCenter;
    newField->_fieldRadius = _fieldRadius;
    // No longer explicitly store a resolution
    newField->_fieldResolution = -1;
    newField->_fieldTimeScale = _fieldTimeScale;
    newField->_fieldSupport = _fieldSupport;
    newField->_fieldSampleRate = _fieldSampleRate;
    newField->_c = _c;
    newField->_interp = NULL;

    // FIXME: debugging
    for ( size_t i = 0; i < newDirections.size(); i++ ) {
        printf( "Direction %lu: ", i );
        cout << newDirections[ i ] << endl;
    }

    numRadii = _fieldData.size() / numDirections;
    new_numDirections = newDirections.size();

    TRACE_ASSERT( _fieldData.size() % numDirections == 0 );

    newField->_fieldData.resize( numRadii * new_numDirections );
    for ( size_t dir_idx = 0; dir_idx < newField->_fieldData.size(); dir_idx++ ) {
        newField->_fieldData[ dir_idx ].resize( signalLength() );
    }

    for ( int dir_idx = 0; dir_idx < new_numDirections; dir_idx++ ) {
        InterpolationDirection   direction;

        direction = interpolationDirection( newDirections[ dir_idx ] );

        for ( int rad_idx = 0; rad_idx < numRadii; rad_idx++ ) {
            int                    dof_idx = dir_idx + rad_idx * new_numDirections;

            for ( int sample_idx = 0; sample_idx < signalLength(); sample_idx++ ) {
                // Use the term_idx argument as an offset
                newField->_fieldData[ dof_idx ][ sample_idx ]
                    = interpolatedFieldSample( sample_idx, direction,
                            rad_idx * numDirections );
            }
        }
    }

    return newField;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void PulseApproximation::ReadAccelerationSet( const std::string &filePrefix,
                                              AccelerationSet &fields )
{
    char                       buf[ 1024 ];

    for ( int dir = TRANS_X; dir < NUM_ACCEL_DIRECTIONS; dir++ )
    {
        sprintf( buf, "%s_%d", filePrefix.c_str(), dir );

        fields._allFields[ dir ] = new PulseApproximation( buf );
    }
}

//////////////////////////////////////////////////////////////////////
// Reads field data from disk
//////////////////////////////////////////////////////////////////////
size_t PulseApproximation::read( const string &filename )
{
    FILE                      *file;
    //int                        size;

    size_t                     bytes_read;

    char                       buf[ 1024 ];

    MATRIX                     dataMatrix;

    sprintf( buf, "%s.dat", filename.c_str() );

    _fieldData.clear();

    file = fopen( buf, "rb" );
    if ( !file )
    {
        cerr << "ERROR: Could not open " << filename << " for reading" << endl;
        abort();
    }

    sprintf( buf, "%s.matrix", filename.c_str() );
    dataMatrix.read( buf );

    _fieldData.resize( dataMatrix.rows() );

    for ( int field_point_idx = 0; field_point_idx < dataMatrix.rows();
            field_point_idx++ )
    {
        _fieldData[ field_point_idx ].resize( dataMatrix.cols() );

        for ( int sample_idx = 0; sample_idx < dataMatrix.cols(); sample_idx++ )
        {
            _fieldData[ field_point_idx ][ sample_idx ]
                = dataMatrix( field_point_idx, sample_idx );
        }
    }

    // Write the remaining info needed to reconstruct the field
    bytes_read = fread( (void *)&_fieldCenter, sizeof( Point3d ), 1, file );
    bytes_read = fread( (void *)&_fieldRadius, sizeof( REAL ), 1, file );
    bytes_read = fread( (void *)&_fieldResolution, sizeof( int ), 1, file );
    bytes_read = fread( (void *)&_fieldTimeScale, sizeof( REAL ), 1, file );
    bytes_read = fread( (void *)&_fieldSupport, sizeof( REAL ), 1, file );
    bytes_read = fread( (void *)&_fieldSampleRate, sizeof( int ), 1, file );

    fclose( file );

    return bytes_read; 
}
