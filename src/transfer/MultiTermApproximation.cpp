//////////////////////////////////////////////////////////////////////
// MultiTermApproximation.cpp: Implementation of the
//                             MultiTermApproximation class
//
//////////////////////////////////////////////////////////////////////

#include "MultiTermApproximation.h"

#include <math/InterpolationFunction.h>

#include <linearalgebra/VECTOR.h>

#include <utils/MathUtil.h>
#include <utils/utils_IO.h>
#include <utils/STLUtil.h>

//////////////////////////////////////////////////////////////////////
// CompactMultiTermApproximation implementation
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
CompactMultiTermApproximation::CompactMultiTermApproximation(
                                        vector<FieldTerm> *fieldTerms[ NUM_ACCEL_DIRECTIONS ],
                                        int baseSample, int signalLength, int startSample,
                                        const Point3d &fieldCenter, REAL fieldRadius,
                                        int fieldResolution,
                                        REAL fieldTimeScale, REAL fieldSupport,
                                        int fieldSampleRate,
                                        REAL c )
    : CompactRadialApproximation( fieldCenter, fieldRadius, fieldResolution,
                                  fieldTimeScale, fieldSupport,
                                  fieldSampleRate, signalLength,
                                  false, false, c ),
      _baseSample( baseSample ),
#if 0
      _signalLength( signalLength ),
#endif
      _startSample( startSample )
{
    printf( "Building CompactMultiTermApproximation with baseSample = %d, "
            "signalLength = %d, startSample = %d, "
            "fieldTimeScale = %f, fieldSupport = %f, "
            "fieldSampleRate = %d\n",
            _baseSample, _signalLength, _startSample,
            fieldTimeScale, fieldSupport, fieldSampleRate );

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

        TRACE_ASSERT( _fields[ direction ]->at( 0 )[ 0 ].size() == _signalLength,
                "Invalid term length" );
    }

    //REAL                       sampleDiv = 1.0 / (REAL)_fieldSampleRate;
    REAL                       sampleDiv = _sampleDiv;

    _interp = new InterpolationMitchellNetravali( sampleDiv );
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
CompactMultiTermApproximation::~CompactMultiTermApproximation()
{
}

//////////////////////////////////////////////////////////////////////
// Assuming the object undergoes a half-sine acceleration pulse
// of the given length, add the contribution of this pulse to
// the given listening position
//////////////////////////////////////////////////////////////////////
bool CompactMultiTermApproximation::addSinePulse( Point3d listeningPosition,
                                                  REAL pulseStartTime,
                                                  REAL pulseLength,
                                                  REAL pulseAmplitude[ NUM_ACCEL_DIRECTIONS ],
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

    REAL                       sampleValue;
    REAL                       sampleScale;
    REAL                       rPower;

    InterpolationDirection     direction;

    VECTOR                     signalStorage( _signalLength );

    int                        numTerms = _fields[ 0 ]->size();
    REAL                       totalScaling;

    if ( !_interp )
    {
        _interp = new InterpolationMitchellNetravali( sampleDiv );
    }

    // Apply any scaling that might be necessary here
    sampleDiv *= scale;
    fieldTimeScale *= scale;
    fieldSupport *= scale;

    InterpolationMitchellNetravali interp( sampleDiv );

    direction = interpolationDirection( listeningPosition );

    // Additional time delay depending on our distance from the receiver
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

        int                      directions[] = { direction._psi_idx_11,
            direction._psi_idx_21,
            direction._psi_idx_12,
            direction._psi_idx_22 };

        signalStorage.clear();

        // Add the contribution from each direction
        for ( int direction_idx = 0; direction_idx < 4; direction_idx++ )
        {
            REAL                   interpCoef = interpolationCoefficient(
                    direction, direction_idx );

            // Add all NUM_ACCEL_DIRECTIONS fields in this direction
            for ( int field_id = 0; field_id < NUM_ACCEL_DIRECTIONS; field_id++ )
            {
                rPower = direction._r;

                // Add each term in the radial series
                for ( int term_idx = 0; term_idx < numTerms; term_idx++ )
                {
                    const FloatArray  &signal
                        = _fields[ field_id ]->at( term_idx )[ directions[ direction_idx ] ];

                    // This gives the scaling for this particular angular
                    totalScaling
                        = sampleScale * pulseAmplitude[ field_id ] * interpCoef / rPower;

                    // The scaling varies depending on whether or not this is a
                    // translational field or a rotational field
                    if ( field_id < ROT_X ) {
                        totalScaling *= transScale;
                    } else {
                        totalScaling *= rotScale;
                    }

                    // Add to the output signal
                    for ( int i = 0; i < _signalLength; i++ ) {
                        signalStorage( i ) += totalScaling * signal[ i ];
                    }

                    rPower *= direction._r;
                }
            }
        }

        // Add to the output
        REAL startTime = tDelay + delay + sampleDiv * (REAL)( _baseSample );
        outputSignal.addSignal( &interp, signalStorage.data(), signalStorage.size(),
                startTime );
    }

    return true;
}

//////////////////////////////////////////////////////////////////////
// MultiTermApproximation implementation
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
MultiTermApproximation::MultiTermApproximation( const FloatArray &radii, int numTerms,
                                                const string &shellDataFilePrefix,
                                                const Point3d &fieldCenter,
                                                int fieldResolution,
                                                REAL fieldTimeScale,
                                                REAL fieldSupport,
                                                int fieldSampleRate,
                                                int field_id,
                                                REAL svdTolerance,
                                                REAL c )
    : RadialApproximation( fieldCenter, 0.0, fieldResolution,
                           fieldTimeScale, fieldSupport, fieldSampleRate,
                           field_id, c ),
      _baseSample( 0 ),
      _signalLength( 0 ),
      _startSample( 0 )
{
    initFieldTerms( shellDataFilePrefix, radii, numTerms, field_id, svdTolerance );
}

//////////////////////////////////////////////////////////////////////
// Construct by reading from a file
//////////////////////////////////////////////////////////////////////
MultiTermApproximation::MultiTermApproximation( const string &filename, REAL c )
    : RadialApproximation( c ),
      _baseSample( 0 ),
      _signalLength( 0 ),
      _startSample( 0 )
{
    read( filename );
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
MultiTermApproximation::~MultiTermApproximation()
{
}

//////////////////////////////////////////////////////////////////////
// Saves the field to disk
//////////////////////////////////////////////////////////////////////
void MultiTermApproximation::write( const string &filename ) const
{
    FILE                      *file;
    int                        size;
    int                        numTerms;

    size_t                     bytes_written;

    char                       buf[ 1024 ];

    sprintf( buf, "%s.dat", filename.c_str() );

    file = fopen( buf, "wb" );
    if ( !file )
    {
        cerr << "ERROR: Could not open " << filename << " for writing" << endl;
        abort();
    }

    numTerms = _fieldTerms.size();
    bytes_written = fwrite( (void *)&numTerms, sizeof( int ), 1, file );

    for ( int term_idx = 0; term_idx < _fieldTerms.size(); term_idx++ )
    {
        MATRIX                   dataMatrix( _fieldTerms[ term_idx ].size(),
                signalLength() );

        // Copy field data to the matrix
        for ( int field_point_idx = 0;
                field_point_idx < _fieldTerms[ term_idx ].size();
                field_point_idx++ )
            for ( int sample_idx = 0; sample_idx < signalLength(); sample_idx++ )
            {
                dataMatrix( field_point_idx, sample_idx )
                    = _fieldTerms[ term_idx ][ field_point_idx ][ sample_idx ];
            }

        sprintf( buf, "%s_term_%d.matrix", filename.c_str(), term_idx );
        dataMatrix.write( buf );
    }

    // Write the remaining info needed to reconstruct the field
    bytes_written = fwrite( (void *)&_fieldCenter, sizeof( Point3d ), 1, file );
    bytes_written = fwrite( (void *)&_fieldRadius, sizeof( REAL ), 1, file );
    bytes_written = fwrite( (void *)&_fieldResolution, sizeof( int ), 1, file );
    bytes_written = fwrite( (void *)&_fieldTimeScale, sizeof( REAL ), 1, file );
    bytes_written = fwrite( (void *)&_fieldSupport, sizeof( REAL ), 1, file );
    bytes_written = fwrite( (void *)&_fieldSampleRate, sizeof( int ), 1, file );
    bytes_written = fwrite( (void *)&_baseSample, sizeof( int ), 1, file );
    bytes_written = fwrite( (void *)&_signalLength, sizeof( int ), 1, file );
    bytes_written = fwrite( (void *)&_startSample, sizeof( int ), 1, file );

    fclose( file );
}

//////////////////////////////////////////////////////////////////////
// Assuming the object undergoes a half-sine acceleration pulse
// of the given length, add the contribution of this pulse to
// the given listening position
//////////////////////////////////////////////////////////////////////
bool MultiTermApproximation::addSinePulse( Point3d listeningPosition,
                                           REAL pulseStartTime,
                                           REAL pulseLength,
                                           REAL pulseAmplitude,
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

    if ( !_interp )
    {
        _interp = new InterpolationMitchellNetravali( sampleDiv );
    }

    // Apply any scaling that might be necessary here
    sampleDiv *= scale;
    fieldTimeScale *= scale;
    fieldSupport *= scale;

    InterpolationMitchellNetravali interp( sampleDiv );

    direction = interpolationDirection( listeningPosition );

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
        cout << SDUMP( fieldSupport ) << endl;
        cout << SDUMP( _baseSample ) << endl;
        cout << SDUMP( _fieldCenter ) << endl;
        cout << SDUMP( _c ) << endl;

        // Add contributions from this pulse to the output signal
        for ( int sample_idx = 0; sample_idx < signalLength(); sample_idx++ )
        {
            t = sampleDiv * (REAL)( _baseSample + sample_idx );
            t += tDelay + delay;

            // Add each term
            rPower = direction._r;
            for ( int term_idx = 0; term_idx < _fieldTerms.size(); term_idx++ )
            {
                sampleValue = interpolatedFieldSample( sample_idx, direction,
                        term_idx );
                sampleValue /= rPower;

                //outputSignal.addSample( t, sampleValue * sampleScale * pulseAmplitude );
                outputSignal.addSample( &interp, t,
                        sampleValue * sampleScale * pulseAmplitude );

                rPower *= direction._r;
            }
        }
    }

    return true;
}

//////////////////////////////////////////////////////////////////////
// Evaluates the field stored by this approximation at the given
// listening position
//////////////////////////////////////////////////////////////////////
void MultiTermApproximation::evaluatePANField( Point3d listeningPosition,
                                               SampledFunction &outputSignal )
{
    REAL                       delay;
    REAL                       sampleDiv = _sampleDiv;
    REAL                       sampleValue;
    REAL                       t;

    REAL                       rPower;

    InterpolationDirection     direction;

    InterpolationMitchellNetravali interp( sampleDiv );

    // Get the angles to interpolate between
    direction = interpolationDirection( listeningPosition );

    delay = direction._r / _c;

    rPower = direction._r;

    for ( int term_idx = 0; term_idx < _fieldTerms.size(); term_idx++ ) {
        for ( int sample_idx = 0; sample_idx < signalLength(); sample_idx++ ) {
            t = sampleDiv * (REAL)( _baseSample + sample_idx );
            t += delay;

            sampleValue = interpolatedFieldSample( sample_idx, direction, term_idx );
            sampleValue /= rPower;

            outputSignal.addSample( &interp, t, sampleValue );
        }

        rPower *= direction._r;
    }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
CompactRadialApproximation *
MultiTermApproximation::buildCompactApproximation( AccelerationSet &allFields,
                                                   const string *directionSetFile )
{
    vector<FieldTerm>         *fields[ NUM_ACCEL_DIRECTIONS ];

    for ( int dir = 0; dir < NUM_ACCEL_DIRECTIONS; dir++ )
    {
        // Yes - this is a bit ugly
        fields[ dir ]
            = ( vector<FieldTerm> * )allFields._allFields[ dir ]->data();
    }

    return new CompactMultiTermApproximation( fields,
            _baseSample, _signalLength,
            _startSample,
            _fieldCenter, _fieldRadius,
            _fieldResolution,
            _fieldTimeScale, _fieldSupport,
            _fieldSampleRate, _c );
}

//////////////////////////////////////////////////////////////////////
// Builds a new multi-term approximation by resampling this one to a new
// angular resolution
//////////////////////////////////////////////////////////////////////
MultiTermApproximation *
MultiTermApproximation::resampleApproximation( int newResolution ) const
{
    MultiTermApproximation    *newField = new MultiTermApproximation();
    Vector3Array               newDirections;

    newField->_fieldCenter = _fieldCenter;
    newField->_fieldRadius = _fieldRadius;
    newField->_fieldResolution = newResolution;
    newField->_fieldTimeScale = _fieldTimeScale;
    newField->_fieldSupport = _fieldSupport;
    newField->_fieldSampleRate = _fieldSampleRate;
    newField->_c = _c;
    newField->_interp = NULL;
    newField->_baseSample = _baseSample;
    newField->_signalLength = _signalLength;
    newField->_startSample = _startSample;

    // Get the new direction set
    MathUtil::GenerateSpherePoints( _fieldCenter, 1.0, newResolution,
            newDirections );

    newField->_fieldTerms.resize( _fieldTerms.size() );
    for ( int term_idx = 0; term_idx < _fieldTerms.size(); term_idx++ )
    {
        newField->_fieldTerms[ term_idx ].resize( newDirections.size() );

        for ( int direction_idx = 0; direction_idx < newDirections.size();
                direction_idx++ )
        {
            newField->_fieldTerms[ term_idx ][ direction_idx ].resize(
                    signalLength() );
        }
    }

    for ( int dof_idx = 0; dof_idx < newDirections.size(); dof_idx++ )
    {
        InterpolationDirection   direction;

        direction = interpolationDirection( newDirections[ dof_idx ] );

        for ( int term_idx = 0; term_idx < _fieldTerms.size(); term_idx++ )
            for ( int sample_idx = 0; sample_idx < signalLength(); sample_idx++ )
            {
                newField->_fieldTerms[ term_idx ][ dof_idx ][ sample_idx ]
                    = interpolatedFieldSample( sample_idx, direction, term_idx );
            }
    }

    return newField;
}

//////////////////////////////////////////////////////////////////////
// Builds the least squares system matrix for a multi-term expansion
//////////////////////////////////////////////////////////////////////
IndexRange MultiTermApproximation::BuildMultiTermSystem( const FloatArray &inputRadii,
                                                         int numSamples, int numTerms,
                                                         int fieldSampleRate,
                                                         MATRIX &systemMatrix,
                                                         REAL c )
{
    int                        startIndex_full = 0;
    int                        endIndex_full = 0;
    int                        numRows, numCols, numTermCols;
    REAL                       listeningRadius;
    IndexRange                 columnRange;

    // The number of rows is just the number of samples in the input
    // fields, multiplied by the number of radii for which we have
    // input values
    numRows = numSamples * inputRadii.size();

    // Determine the number of columns that we need
    for ( int shell_idx = 0; shell_idx < inputRadii.size(); shell_idx++ )
    {
        int                      startIndex, endIndex;

        listeningRadius = inputRadii[ shell_idx ];

        MultiTermApproximation::GetRadiusIndexRange( listeningRadius,
                numSamples,
                fieldSampleRate,
                startIndex, endIndex, c );

        startIndex_full = min( startIndex, startIndex_full );
        endIndex_full = max( endIndex, endIndex_full );
    }

    numTermCols = endIndex_full - startIndex_full + 1;
    numCols = numTermCols * numTerms;

    systemMatrix.resizeAndWipe( numRows, numCols );

    // Add each block row
    for ( int shell_idx = 0; shell_idx < inputRadii.size(); shell_idx++ )
    {
        int                      start_row = shell_idx * numSamples;

        // Add each column block
        for ( int term_idx = 0; term_idx < numTerms; term_idx++ )
        {
            int                    start_col = term_idx * numTermCols;

            MultiTermApproximation::AddRadiusTerms( inputRadii[ shell_idx ],
                    numSamples,
                    term_idx + 1, /* Radius exp. */
                    fieldSampleRate,
                    systemMatrix,
                    startIndex_full, /* Base column */
                    start_row, start_col, c );
        }
    }

    columnRange.first = startIndex_full;
    columnRange.second = endIndex_full;

    return columnRange;
}

//////////////////////////////////////////////////////////////////////
// Reorders field data so that all radii for a single shell
// are put together in the same column
//////////////////////////////////////////////////////////////////////
void MultiTermApproximation::BuildMultiTermRHS( const MATRIX &fieldData,
                                                MATRIX &rhs,
                                                int numRadii )
{
    int                        pointsPerShell;
    int                        numSamples;

    TRACE_ASSERT( fieldData.rows() % numRadii == 0, "Invalid field data" );

    pointsPerShell = fieldData.rows() / numRadii;
    numSamples = fieldData.cols();

    rhs.resizeAndWipe( numSamples * numRadii, pointsPerShell );

    for ( int row_idx = 0; row_idx < fieldData.rows(); row_idx++ )
    {
        int                      shell_idx = row_idx / pointsPerShell;
        int                      direction_idx = row_idx % pointsPerShell;
        int                      base_row = shell_idx * numSamples;

        // Place row data from the input matrix along a row of the
        // new RHS matrix
        for ( int col_idx = 0; col_idx < fieldData.cols(); col_idx++ )
        {
            rhs( base_row + col_idx, direction_idx ) = fieldData( row_idx, col_idx );
        }
    }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void MultiTermApproximation::ReadAccelerationSet( const std::string &filePrefix,
                                                  AccelerationSet &fields )
{
    char                       buf[ 1024 ];

    for ( int dir = TRANS_X; dir < NUM_ACCEL_DIRECTIONS; dir++ )
    {
        sprintf( buf, "%s_%d", filePrefix.c_str(), dir );

        fields._allFields[ dir ] = new MultiTermApproximation( buf );
    }
}

//////////////////////////////////////////////////////////////////////
// When building a multi-term expansion, adds matrix entries for
// a single radius value
//////////////////////////////////////////////////////////////////////
void MultiTermApproximation::AddRadiusTerms( REAL radius, int numSamples,
                                             int radiusExponent,
                                             int fieldSampleRate,
                                             SPARSE_MATRIX &systemMatrix,
                                             int base_column,
                                             int start_row,
                                             REAL c )
{
    REAL                             h = 1.0 / (REAL)fieldSampleRate;
    REAL                             radiusPower = radius;

    // The interpolation function to use
    InterpolationMitchellNetravali   f( h );

    // Make sure we have enough space for the first and last sample
    int                              startIndex, endIndex;

    GetRadiusIndexRange( radius, numSamples, fieldSampleRate,
            startIndex, endIndex, c );

    TRACE_ASSERT( startIndex >= base_column );
    TRACE_ASSERT( endIndex - base_column < systemMatrix.cols() );
    TRACE_ASSERT( start_row + numSamples <= systemMatrix.rows() );

    for ( int i = 1; i < radiusExponent; i++ )
    {
        radiusPower *= radius;
    }

    // Fill in each system row
    for ( int sample_idx = 0; sample_idx < numSamples; sample_idx++ )
    {
        REAL             sampleTime = h * (REAL)sample_idx - radius / c;
        int              sampleStart;
        int              sampleEnd;

        sampleStart = (int)floor( sampleTime / h ) - f.support();
        sampleEnd = (int)floor( sampleTime / h ) + f.support();

        for ( int col_idx = sampleStart; col_idx <= sampleEnd; col_idx++ )
        {
            REAL           interpSampleTime = h * (REAL)col_idx;
            REAL           interpValue = f.evaluate( sampleTime, interpSampleTime );

            interpValue /= radiusPower;

            systemMatrix( start_row + sample_idx, col_idx - base_column )
                = interpValue;
        }
    }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void MultiTermApproximation::AddRadiusTerms( REAL radius, int numSamples,
                                             int radiusExponent,
                                             int fieldSampleRate,
                                             MATRIX &systemMatrix,
                                             int base_column,
                                             int start_row, int start_col,
                                             REAL c )
{
    REAL                             h = 1.0 / (REAL)fieldSampleRate;
    REAL                             radiusPower = radius;

    // The interpolation function to use
    InterpolationMitchellNetravali   f( h );

    // Make sure we have enough space for the first and last sample
    int                              startIndex, endIndex;

    GetRadiusIndexRange( radius, numSamples, fieldSampleRate,
            startIndex, endIndex, c );

    TRACE_ASSERT( startIndex >= base_column );
    TRACE_ASSERT( endIndex - base_column + start_col  < systemMatrix.cols() );
    TRACE_ASSERT( start_row + numSamples <= systemMatrix.rows() );

    for ( int i = 1; i < radiusExponent; i++ )
    {
        radiusPower *= radius;
    }

    // Fill in each system row
    for ( int sample_idx = 0; sample_idx < numSamples; sample_idx++ )
    {
        REAL             sampleTime = h * (REAL)sample_idx - radius / c;
        int              sampleStart;
        int              sampleEnd;

        sampleStart = (int)floor( sampleTime / h ) - f.support();
        sampleEnd = (int)floor( sampleTime / h ) + f.support();

        for ( int col_idx = sampleStart; col_idx <= sampleEnd; col_idx++ )
        {
            REAL           interpSampleTime = h * (REAL)col_idx;
            REAL           interpValue = f.evaluate( sampleTime, interpSampleTime );

            interpValue /= radiusPower;

            systemMatrix( start_row + sample_idx, col_idx - base_column + start_col )
                = interpValue;
        }
    }
}

//////////////////////////////////////////////////////////////////////
// Helper function for the above
//////////////////////////////////////////////////////////////////////
void MultiTermApproximation::GetRadiusIndexRange( REAL radius,
                                                  int numSamples,
                                                  int fieldSampleRate,
                                                  int &startIndex,
                                                  int &endIndex,
                                                  REAL c )
{
    REAL                             h = 1.0 / (REAL)fieldSampleRate;

    // The interpolation function to use
    InterpolationMitchellNetravali   f( h );

    // Make sure we have enough space for the first and last sample
    REAL                             tStart, tEnd;

    tStart = -radius / c;
    tEnd = h * (REAL)( numSamples - 1 ) - radius / c;

    startIndex = (int)floor( tStart / h ) - f.support();
    endIndex = (int)floor( tEnd / h ) + f.support();
}

//////////////////////////////////////////////////////////////////////
// Reads field data from disk
//////////////////////////////////////////////////////////////////////
size_t MultiTermApproximation::read( const string &filename )
{
    FILE                      *file;
    int                        size;
    int                        numTerms;

    size_t                     bytes_read;

    char                       buf[ 1024 ];

    MATRIX                     dataMatrix;

    sprintf( buf, "%s.dat", filename.c_str() );

    _fieldTerms.clear();

    file = fopen( buf, "rb" );
    if ( !file )
    {
        cerr << "ERROR: Could not open " << filename << " for reading" << endl;
        abort();
    }

    bytes_read = fread( (void *)&numTerms, sizeof( int ), 1, file );
    _fieldTerms.resize( numTerms );

    for ( int term_idx = 0; term_idx < _fieldTerms.size(); term_idx++ )
    {
        sprintf( buf, "%s_term_%d.matrix", filename.c_str(), term_idx );
        dataMatrix.read( buf );

        printf( "Read term %d with %d points and %d samples\n",
                term_idx, dataMatrix.rows(), dataMatrix.cols() );

        _fieldTerms[ term_idx ].resize( dataMatrix.rows() );

        for ( int field_point_idx = 0; field_point_idx < dataMatrix.rows();
                field_point_idx++ )
        {
            _fieldTerms[ term_idx ][ field_point_idx ].resize( dataMatrix.cols() );

            for ( int sample_idx = 0; sample_idx < dataMatrix.cols(); sample_idx++ )
            {
                _fieldTerms[ term_idx ][ field_point_idx ][ sample_idx ]
                    = dataMatrix( field_point_idx, sample_idx );
            }
        }
    }

    // Write the remaining info needed to reconstruct the field
    bytes_read = fread( (void *)&_fieldCenter, sizeof( Point3d ), 1, file );
    bytes_read = fread( (void *)&_fieldRadius, sizeof( REAL ), 1, file );
    bytes_read = fread( (void *)&_fieldResolution, sizeof( int ), 1, file );
    bytes_read = fread( (void *)&_fieldTimeScale, sizeof( REAL ), 1, file );
    bytes_read = fread( (void *)&_fieldSupport, sizeof( REAL ), 1, file );
    bytes_read = fread( (void *)&_fieldSampleRate, sizeof( int ), 1, file );
    bytes_read = fread( (void *)&_baseSample, sizeof( int ), 1, file );
    bytes_read = fread( (void *)&_signalLength, sizeof( int ), 1, file );
    bytes_read = fread( (void *)&_startSample, sizeof( int ), 1, file );

    fclose( file );

    return bytes_read;
}

//////////////////////////////////////////////////////////////////////
// Initializes the content of _fieldTerms
//
// Performs a least squares solve to find the best solution for
// the given number of terms
//////////////////////////////////////////////////////////////////////
void MultiTermApproximation::initFieldTerms( const string &shellDataFilePrefix,
                                             const FloatArray &radii, int numTerms,
                                             int field_id, REAL svdTolerance )
{
    MATRIX                     rhs;
    MATRIX                     systemMatrix;
    VECTOR                     singularValues;
    int                        rank;
    int                        numSamples;
    IndexRange                 columnRange;
    char                       filename[ 1024 ];

    // Read in the right hand side matrix, then reorder it's entries to put
    // them in the format we need
    {
        MATRIX                     fieldData;

        sprintf( filename, "%s_%d.matrix", shellDataFilePrefix.c_str(), field_id );
        fieldData.read( filename );

        cout << "Building RHS" << endl;

        BuildMultiTermRHS( fieldData, rhs, radii.size() );

        numSamples = rhs.rows() / radii.size();

        TRACE_ASSERT( rhs.rows() % radii.size() == 0 );
    }

    cout << "Building system" << endl;
    columnRange = BuildMultiTermSystem( radii, numSamples, numTerms,
            _fieldSampleRate, systemMatrix, _c );

    // Get the solution
    systemMatrix.write( "test.matrix" );
    rhs.write( "testRHS.matrix" );
    cout << "Least squares solve " << endl;
    cout << SDUMP( systemMatrix.rows() ) << "   "
        << SDUMP( systemMatrix.cols() ) << endl;
    MATRIX::LeastSquares_TSVD( systemMatrix, rhs, singularValues, rank,
            svdTolerance );
    cout << "Done" << endl;

    _baseSample = columnRange.first;

    extractFieldTerms( rhs, systemMatrix, columnRange, numTerms );
}

//////////////////////////////////////////////////////////////////////
// Helper function for init field terms.  Once the least squares
// solve is done, this packs the matrix data from the solution in
// to vectors in _fieldTerms
//////////////////////////////////////////////////////////////////////
void MultiTermApproximation::extractFieldTerms( const MATRIX &solvedData,
                                                const MATRIX &systemMatrix,
                                                const IndexRange &columnRange,
                                                int numTerms )
{
    int                        termLength = systemMatrix.cols() / numTerms;

    TRACE_ASSERT( systemMatrix.cols() % termLength == 0 );
    TRACE_ASSERT( solvedData.rows() >= systemMatrix.cols() );

    // TODO: At some point we should prune out zeros in the signal,
    // but for now we'll just take the whole thing
    _baseSample = columnRange.first;
    _startSample = 0;

    _signalLength = termLength;

    _fieldTerms.resize( numTerms, FieldTerm( solvedData.cols() ) );

    for ( int col_idx = 0; col_idx < solvedData.cols(); col_idx++ ) {
        for ( int term_idx = 0; term_idx < numTerms; term_idx++ )
        {
            int                    base_idx = term_idx * termLength;

            _fieldTerms[ term_idx ][ col_idx ].resize( termLength );

            for ( int sample_idx = 0; sample_idx < termLength; sample_idx++ )
            {
                // The column index is the index of the direction we
                // are filling in here
                _fieldTerms[ term_idx ][ col_idx ][ sample_idx ]
                    = solvedData( base_idx + sample_idx, col_idx );
            }
        }
    }
}


