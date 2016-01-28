//////////////////////////////////////////////////////////////////////
// RadialApproximation.cpp: Implementation of the RadialApproximation
//                          class
//
//////////////////////////////////////////////////////////////////////

#include "RadialApproximation.h"

#include <math/InterpolationFunction.h>

#include <utils/MathUtil.h>
#include <utils/IO.h>

//////////////////////////////////////////////////////////////////////
// RadialField implementation
//////////////////////////////////////////////////////////////////////

RadialField::RadialField( const Point3d &fieldCenter, REAL fieldRadius, int fieldResolution )
    : _fieldCenter( fieldCenter ),
      _fieldRadius( fieldRadius ),
      _fieldResolution( fieldResolution )
{
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
RadialField::~RadialField()
{
}

//////////////////////////////////////////////////////////////////////
// Determines an interpolation direction for the given
// point in 3D space
//////////////////////////////////////////////////////////////////////
RadialField::InterpolationDirection
RadialField::interpolationDirection( const Point3d &x ) const
{
    InterpolationDirection     result;

    // Get the spherical coordinates of this point
    Point3d                    R;

    R = MathUtil::spherical_coordinates( x - _fieldCenter );

    // Find the bounding theta and phi coordinates
    REAL                       r = R[0];
    REAL                       theta = R[1];
    REAL                       phi = R[2];

    REAL                       thetaDiv;
    REAL                       phiDiv;

    int                        fieldDOFS;

    thetaDiv = M_PI / (REAL)( _fieldResolution + 1 );
    phiDiv = 2.0 * M_PI / (REAL)( 2 * _fieldResolution );

    fieldDOFS = 2 * _fieldResolution * _fieldResolution + 2;

    if ( theta < 0 )
    {
        theta = fabs( theta );
    }

    if ( phi < 0 )
    {
        phi += 2.0 * M_PI;
    }

    int                        theta_idx_1 = (int)floor( theta / thetaDiv );
    int                        phi_idx_1 = (int)floor( phi / phiDiv );
    int                        theta_idx_2, phi_idx_2;
    int                        psi_idx_11, psi_idx_12, psi_idx_21, psi_idx_22;
    REAL                       theta_1, theta_2, phi_1, phi_2;

    // Clamp these to the appropriate range
    if ( theta_idx_1 < 0 )
    {
        theta_idx_1 = 0;
    }
    else if ( theta_idx_1 > _fieldResolution )
    {
        theta_idx_1 = _fieldResolution;
    }

    if ( phi_idx_1 < 0 )
    {
        phi_idx_2 = 0;
    }
    else if ( phi_idx_1 > 2 * _fieldResolution - 1 )
    {
        phi_idx_1 = 2 * _fieldResolution - 1;
    }

    theta_idx_2 = theta_idx_1 + 1;
    phi_idx_2 = phi_idx_1 + 1;

    theta_1 = (REAL)theta_idx_1 * thetaDiv;
    theta_2 = (REAL)theta_idx_2 * thetaDiv;
    phi_1 = (REAL)phi_idx_1 * phiDiv;
    phi_2 = (REAL)phi_idx_2 * phiDiv;

    if ( phi_idx_2 == 2 * _fieldResolution )
    {
        phi_idx_2 = 0;
    }

    // Find the starting index for each of the four needed
    // angular coordinates.
    psi_idx_11 = ( theta_idx_1 - 1 ) * ( 2 * _fieldResolution ) + 1 + phi_idx_1;
    psi_idx_12 = ( theta_idx_1 - 1 ) * ( 2 * _fieldResolution ) + 1 + phi_idx_2;
    psi_idx_21 = ( theta_idx_2 - 1 ) * ( 2 * _fieldResolution ) + 1 + phi_idx_1;
    psi_idx_22 = ( theta_idx_2 - 1 ) * ( 2 * _fieldResolution ) + 1 + phi_idx_2;

    // Special cases:
    if ( theta_idx_1 == 0 )
    {
        // If we're near the top of the sphere (theta_idx_1 = 0)
        // then we need the data from the first cell.
        psi_idx_11 = 0;
        psi_idx_12 = 0;
    }
    else if ( theta_idx_1 == _fieldResolution )
    {
        // If we're near the bottom of the sphere (theta_idx_1 = resolution)
        // then we need the data from the last cell.
        psi_idx_21 = fieldDOFS - 1;
        psi_idx_22 = psi_idx_21;
    }

    result._thetaDiv = thetaDiv;
    result._phiDiv = phiDiv;
    result._theta_1 = theta_1;
    result._theta_2 = theta_2;
    result._phi_1 = phi_1;
    result._phi_2 = phi_2;
    result._psi_idx_11 = psi_idx_11;
    result._psi_idx_12 = psi_idx_12;
    result._psi_idx_21 = psi_idx_21;
    result._psi_idx_22 = psi_idx_22;
    result._theta = theta;
    result._phi = phi;
    result._r = r;

    return result;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
RadialField::InterpolationDirection
RadialField::octantInterpolationDirection( const Point3d &x ) const
{
    InterpolationDirection     result;

    // Get the spherical coordinates of this point
    Point3d                    R;

    R = MathUtil::spherical_coordinates( x - _fieldCenter );

    // Find the bounding theta and phi coordinates
    REAL                       r = R[0];
    REAL                       theta = R[1];
    REAL                       phi = R[2];

    REAL                       thetaDiv;
    REAL                       phiDiv;

    //int                        fieldDOFS;

    int                        azimuthResolution;

    azimuthResolution = _fieldResolution + 1;

    thetaDiv = M_PI / 2.0 / (REAL)( _fieldResolution );
    phiDiv = M_PI / 2.0 / (REAL)( azimuthResolution - 1 );

    //fieldDOFS = _fieldResolution * ( azimuthResolution ) + 1;

    if ( theta < 0 )
    {
        theta = fabs( theta );
    }

    if ( phi < 0 )
    {
        phi += 2.0 * M_PI;
    }

    TRACE_ASSERT( theta >= 0.0 );
    TRACE_ASSERT( theta <= M_PI / 2.0 );
    TRACE_ASSERT( phi >= 0.0 );
    TRACE_ASSERT( phi <= M_PI / 2.0 );

    int                        theta_idx_1 = (int)floor( theta / thetaDiv );
    int                        phi_idx_1 = (int)floor( phi / phiDiv );
    int                        theta_idx_2, phi_idx_2;
    int                        psi_idx_11, psi_idx_12, psi_idx_21, psi_idx_22;
    REAL                       theta_1, theta_2, phi_1, phi_2;

    // Clamp these to the appropriate range
    if ( theta_idx_1 < 0 ) {
        theta_idx_1 = 0;
    } else if ( theta_idx_1 > _fieldResolution ) {
        theta_idx_1 = _fieldResolution;
    }

    if ( phi_idx_1 < 0 ) {
        phi_idx_1 = 0;
    } else if ( phi_idx_1 > azimuthResolution - 1 ) {
        phi_idx_1 = azimuthResolution - 1;
    }

    theta_idx_2 = theta_idx_1 + 1;
    phi_idx_2 = phi_idx_1 + 1;

    // Clamp
    if ( theta_idx_2 > _fieldResolution ) {
        theta_idx_2 = _fieldResolution;
    }
    if ( phi_idx_2 > azimuthResolution - 1 ) {
        phi_idx_2 = azimuthResolution - 1;
    }

    theta_1 = (REAL)theta_idx_1 * thetaDiv;
    theta_2 = (REAL)theta_idx_2 * thetaDiv;
    phi_1 = (REAL)phi_idx_1 * phiDiv;
    phi_2 = (REAL)phi_idx_2 * phiDiv;

    // Find the starting index for each of the four needed
    // angular coordinates.
    psi_idx_11 = ( theta_idx_1 - 1 ) * ( azimuthResolution ) + 1 + phi_idx_1;
    psi_idx_12 = ( theta_idx_1 - 1 ) * ( azimuthResolution ) + 1 + phi_idx_2;
    psi_idx_21 = ( theta_idx_2 - 1 ) * ( azimuthResolution ) + 1 + phi_idx_1;
    psi_idx_22 = ( theta_idx_2 - 1 ) * ( azimuthResolution ) + 1 + phi_idx_2;

    // Special cases:
    if ( theta_idx_1 == 0 ) {
        // If we're near the top of the sphere (theta_idx_1 = 0)
        // then we need the data from the first cell.
        psi_idx_11 = 0;
        psi_idx_12 = 0;
    }

    result._thetaDiv = thetaDiv;
    result._phiDiv = phiDiv;
    result._theta_1 = theta_1;
    result._theta_2 = theta_2;
    result._phi_1 = phi_1;
    result._phi_2 = phi_2;
    result._psi_idx_11 = psi_idx_11;
    result._psi_idx_12 = psi_idx_12;
    result._psi_idx_21 = psi_idx_21;
    result._psi_idx_22 = psi_idx_22;
    result._theta = theta;
    result._phi = phi;
    result._r = r;

    return result;
}

//////////////////////////////////////////////////////////////////////
// Interpolates the field at the given time sample
//////////////////////////////////////////////////////////////////////
REAL RadialField::interpolatedFieldSample( int sample_idx,
                                           const InterpolationDirection &direction,
                                           int term_idx, int field_idx ) const
{
    REAL                       psi_11;
    REAL                       psi_21;
    REAL                       psi_12;
    REAL                       psi_22;

    REAL                       theta_1 = direction._theta_1;
    REAL                       theta_2 = direction._theta_2;
    REAL                       phi_1 = direction._phi_1;
    REAL                       phi_2 = direction._phi_2;

    REAL                       theta = direction._theta;
    REAL                       phi = direction._phi;

    REAL                       thetaDiv = direction._thetaDiv;
    REAL                       phiDiv = direction._phiDiv;

    psi_11 = sampleDirection( direction._psi_idx_11, sample_idx,
            term_idx, field_idx );
    psi_21 = sampleDirection( direction._psi_idx_21, sample_idx,
            term_idx, field_idx );
    psi_12 = sampleDirection( direction._psi_idx_12, sample_idx,
            term_idx, field_idx );
    psi_22 = sampleDirection( direction._psi_idx_22, sample_idx,
            term_idx, field_idx );

    return ( psi_11 * ( theta_2 - theta ) * ( phi_2 - phi )
            + psi_21 * ( theta - theta_1 ) * ( phi_2 - phi )
            + psi_12 * ( theta_2 - theta ) * ( phi - phi_1 )
            + psi_22 * ( theta - theta_1 ) * ( phi - phi_1 ) )
        / ( thetaDiv * phiDiv );
}

//////////////////////////////////////////////////////////////////////
// Interpolation coefficients for each of the 4 directions
//
//   corner_idx:   0 <--> 11
//                 1 <--> 21
//                 2 <--> 12
//                 3 <--> 22
//////////////////////////////////////////////////////////////////////
REAL RadialField::interpolationCoefficient( const InterpolationDirection &direction,
                                            int corner_idx )
{

    REAL                       theta_1 = direction._theta_1;
    REAL                       theta_2 = direction._theta_2;
    REAL                       phi_1 = direction._phi_1;
    REAL                       phi_2 = direction._phi_2;

    REAL                       theta = direction._theta;
    REAL                       phi = direction._phi;

    REAL                       thetaDiv = direction._thetaDiv;
    REAL                       phiDiv = direction._phiDiv;

    switch ( corner_idx )
    {
        case 0:
        default:
            {
                return ( theta_2 - theta ) * ( phi_2 - phi ) / ( thetaDiv * phiDiv );
            }
        case 1:
            {
                return ( theta - theta_1 ) * ( phi_2 - phi ) / ( thetaDiv * phiDiv );
            }
        case 2:
            {
                return ( theta_2 - theta ) * ( phi - phi_1 ) / ( thetaDiv * phiDiv );
            }
        case 3:
            {
                return ( theta - theta_1 ) * ( phi - phi_1 ) / ( thetaDiv * phiDiv );
            }
    }
}

//////////////////////////////////////////////////////////////////////
// For fields with symmetric geometry about each axis, we can consider
// only the values of the field in the first quadrant (ie. x, y, z > 0).
// This function returns the appropriate scaling factor for a given
// position and field ID.
//
// It also modifies the listening position to put it in the first
// quadrant
//////////////////////////////////////////////////////////////////////
REAL RadialField::SymmetricScalingFactor( Point3d &listeningPosition,
                                          int field_id )
{
    unsigned int               quadrant = PositionQuadrant( listeningPosition );

    return SymmetricScalingFactor( quadrant, field_id );
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void RadialField::SymmetricScalingFactor( Point3d &listeningPosition,
                                          REAL pulseAmplitudes[ NUM_ACCEL_DIRECTIONS ] )
{
    unsigned int               quadrant = PositionQuadrant( listeningPosition );

    for ( int dir = 0; dir < NUM_ACCEL_DIRECTIONS; dir++ ) {
        pulseAmplitudes[ dir ] *= SymmetricScalingFactor( quadrant, dir );
    }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
REAL RadialField::_SymmetricScalingFactors[][NUM_OCTANTS] =
{
    { 1, 1, -1, -1, 1, 1, -1, -1 },
    { 1, 1, 1, 1, -1, -1, -1, -1 },
    { 1, -1, 1, -1, 1, -1, 1, -1 },
    { 1, -1, 1, -1, -1, 1, -1, 1 },
    { 1, -1, -1, 1, 1, -1, -1, 1 },
    { 1, 1, -1, -1, -1, -1, 1, 1 }
};

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
REAL RadialField::SymmetricScalingFactor( unsigned int quadrant, int field_id )
{
    return _SymmetricScalingFactors[ field_id ][ quadrant ];
#if 0
    switch ( field_id )
    {
        case TRANS_X: {
                          if ( quadrant == 0 || quadrant == 1 || quadrant == 4 || quadrant == 5 ) {
                              return 1.0;
                          } else {
                              return -1.0;
                          }
                      }
        case TRANS_Y: {
                          if ( quadrant == 0 || quadrant == 1 || quadrant == 2 || quadrant == 3 ) {
                              return 1.0;
                          } else {
                              return -1.0;
                          }
                      }
        case TRANS_Z: {
                          if ( quadrant == 0 || quadrant == 2 || quadrant == 4 || quadrant == 6 ) {
                              return 1.0;
                          } else {
                              return -1.0;
                          }
                      }
        case ROT_X: {
                        if ( quadrant == 0 || quadrant == 2 || quadrant == 5 || quadrant == 7 ) {
                            return 1.0;
                        } else {
                            return -1.0;
                        }
                    }
        case ROT_Y: {
                        if ( quadrant == 0 || quadrant == 3 || quadrant == 4 || quadrant == 7 ) {
                            return 1.0;
                        } else {
                            return -1.0;
                        }
                    }
        case ROT_Z: {
                        if ( quadrant == 0 || quadrant == 1 || quadrant == 6 || quadrant == 7 ) {
                            return 1.0;
                        } else {
                            return -1.0;
                        }
                    }
    }

    return 1.0;
#endif
}

//////////////////////////////////////////////////////////////////////
// Returns the quadrant of this position, and reflects it
// in to the first quadrant
//////////////////////////////////////////////////////////////////////
unsigned int RadialField::PositionQuadrant( Point3d &listeningPosition )
{
    unsigned const int         zNegative = 1;
    unsigned const int         xNegative = 2;
    unsigned const int         yNegative = 4;
    unsigned int               quadrant = 0;

    // Construct quadrant bitmask
    if ( listeningPosition[ 0 ] < 0.0 ) {
        quadrant |= xNegative;

        listeningPosition[ 0 ] *= -1.0;
    }

    if ( listeningPosition[ 1 ] < 0.0 ) {
        quadrant |= yNegative;

        listeningPosition[ 1 ] *= -1.0;
    }

    if ( listeningPosition[ 2 ] < 0.0 ) {
        quadrant |= zNegative;

        listeningPosition[ 2 ] *= -1.0;
    }

    return quadrant;
}

//////////////////////////////////////////////////////////////////////
// CompactRadialApproximation implementation
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
CompactRadialApproximation::CompactRadialApproximation( const Point3d &fieldCenter,
                                                        REAL fieldRadius,
                                                        REAL fieldResolution,
                                                        REAL fieldTimeScale, REAL fieldSupport,
                                                        int fieldSampleRate,
                                                        int signalLength,
                                                        bool useSymmetry,
                                                        bool matchSamplingRate,
                                                        REAL c )
: RadialField( fieldCenter, fieldRadius, fieldResolution ),
    _fieldTimeScale( fieldTimeScale ),
    _fieldSupport( fieldSupport ),
    _fieldSampleRate( fieldSampleRate ),
    _signalLength( signalLength ),
    _useSymmetry( useSymmetry ),
    _matchSamplingRate( matchSamplingRate ),
    _c( c ),
    _interp( NULL )
{
    _sampleDiv = 1.0 / (REAL)_fieldSampleRate;
}

//////////////////////////////////////////////////////////////////////
// Assuming the object undergoes a half-sine acceleration pulse
// of the given length, add the contribution of this pulse to
// the given listening position
//////////////////////////////////////////////////////////////////////
bool CompactRadialApproximation::addSinePulse( Point3d listeningPosition,
                                               REAL pulseStartTime,
                                               REAL pulseLength,
                                               REAL pulseAmplitude[ NUM_ACCEL_DIRECTIONS ],
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

    //REAL                       totalScaling;

    REAL                       sampleScale;
    //REAL                       rPower;

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
//////////////////////////////////////////////////////////////////////
CompactRadialApproximation::~CompactRadialApproximation()
{
}

//////////////////////////////////////////////////////////////////////
// RadialApproximation implementation
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
RadialApproximation::RadialApproximation( const Point3d &fieldCenter,
                                          REAL fieldRadius,
                                          int fieldResolution,
                                          REAL fieldTimeScale,
                                          REAL fieldSupport,
                                          int fieldSampleRate,
                                          int field_id,
                                          REAL c )
: RadialField( fieldCenter, fieldRadius, fieldResolution ),
    _fieldTimeScale( fieldTimeScale ),
    _fieldSupport( fieldSupport ),
    _fieldSampleRate( fieldSampleRate ),
    _c( c ),
    _interp( NULL )
{
    _sampleDiv = 1.0 / (REAL)_fieldSampleRate;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
    RadialApproximation::RadialApproximation( REAL c )
: _c( c ),
    _interp( NULL )
{
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
RadialApproximation::~RadialApproximation()
{
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void RadialApproximation::WriteAccelerationSet( const std::string &filePrefix,
                                                const AccelerationSet &fields )
{
    char                       buf[ 1024 ];

    for ( int dir = TRANS_X; dir < NUM_ACCEL_DIRECTIONS; dir++ )
    {
        sprintf( buf, "%s_%d", filePrefix.c_str(), dir );

        fields._allFields[ dir ]->write( buf );
    }
}

//////////////////////////////////////////////////////////////////////
// Adds pulses for all acceleration sets, given a single listening
// position and time scale (but multiple amplitudes)
//////////////////////////////////////////////////////////////////////
bool RadialApproximation::AddSinePulse( AccelerationSet &fields,
                                        const Point3d &listeningPosition,
                                        REAL pulseStartTime, REAL pulseLength,
                                        REAL pulseAmplitude[ NUM_ACCEL_DIRECTIONS ],
                                        SampledFunction &outputSignal )
{
    bool                       addingSound = false;

    // Add the pulse for each acceleration direction
    for ( int dir = TRANS_X; dir < NUM_ACCEL_DIRECTIONS; dir++ )
    {
        bool                     termAdded;

        termAdded = fields._allFields[ dir ]->addSinePulse(
                listeningPosition,
                pulseStartTime, pulseLength,
                pulseAmplitude[ dir ],
                outputSignal );

        addingSound = addingSound || termAdded;
    }

    return addingSound;
}
