//////////////////////////////////////////////////////////////////////
// ScalarField.cpp: Implementation of the ScalarField class
//
//////////////////////////////////////////////////////////////////////

#include "ScalarField.h"

#include <utils/IO.h>

const int ScalarField::NUM_NEIGHBOURS = 6;

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
ScalarField::ScalarField()
    : _bbox( Vector3d( 0, 0, 0 ), Vector3d( 0, 0, 0 ) ),
      _divisions( 1, 1, 1 ),
      _cellSize( 0.0 ),
      _fieldValues( NULL )
{
}

//////////////////////////////////////////////////////////////////////
// Constructs scalar field within the given bounding box,
// with the given division size.  The box may be enlarged
// to ensure square cells.
//////////////////////////////////////////////////////////////////////
ScalarField::ScalarField( const BoundingBox &bbox, REAL cellSize )
    : _cellSize( cellSize ),
      _fieldValues( NULL )
{
    Vector3d                   minBound;
    Vector3d                   maxBound;

    for ( int dimension = 0; dimension < 3; dimension++ )
    {
        _divisions[ dimension ]
            = (int)ceil( bbox.axislength( dimension ) / cellSize );

        minBound[ dimension ] = bbox.minBound()[ dimension ];

        maxBound[ dimension ] = minBound[ dimension ];
        maxBound[ dimension ] += cellSize * (REAL)_divisions[ dimension ];
    }

    _bbox = BoundingBox( minBound, maxBound );
}

//////////////////////////////////////////////////////////////////////
// Builds a scalar field starting at the given minimum bound,
// and with the specified number of cell divisions and cell size
//////////////////////////////////////////////////////////////////////
ScalarField::ScalarField( const Vector3d &minBound, const Tuple3i &divisions,
        REAL cellSize )
: _divisions( divisions ),
    _cellSize( cellSize ),
    _fieldValues( NULL )
{
    Vector3d                 maxBound;

    for ( int dimension = 0; dimension < 3; dimension++ )
    {
        maxBound[ dimension ] = minBound[ dimension ];
        maxBound[ dimension ] += cellSize * (REAL)( _divisions[ dimension ] - 1 );
    }

    _bbox = BoundingBox( minBound, maxBound );
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
ScalarField::~ScalarField()
{
    if ( _fieldValues )
    {
        delete _fieldValues;
    }
}

//////////////////////////////////////////////////////////////////////
// Returns the flat cell index for the given x, y, z index
//////////////////////////////////////////////////////////////////////
int ScalarField::cellIndex( const Tuple3i &index ) const
{
    int                        flatIndex;

    flatIndex = index[ 0 ] * ( _divisions[ 1 ] * _divisions[ 2 ] );
    flatIndex += index[ 1 ] * _divisions[ 2 ];
    flatIndex += index[ 2 ];

    return flatIndex;
}

//////////////////////////////////////////////////////////////////////
// Returns the x, y and z index for the given flat cell index
//////////////////////////////////////////////////////////////////////
Tuple3i ScalarField::cellIndex( int flatIndex ) const
{
    Tuple3i                    index;

    index[ 0 ] = flatIndex / ( _divisions[ 1 ] * _divisions[ 2 ] );

    flatIndex -= index[ 0 ] * ( _divisions[ 1 ] * _divisions[ 2 ] );

    index[ 1 ] = flatIndex / _divisions[ 2 ];

    flatIndex -= index[ 1 ] * _divisions[ 2 ];

    index[ 2 ] = flatIndex;

    return index;
}

//////////////////////////////////////////////////////////////////////
// Cell center for the given index
//////////////////////////////////////////////////////////////////////
Vector3d ScalarField::cellPosition( const Tuple3i &index ) const
{
    Vector3d                 pos;

    for ( int dimension = 0; dimension < 3; dimension++ )
    {
        pos[ dimension ] = _bbox.minBound()[ dimension ];
        pos[ dimension ] += ( 0.5 + (REAL)index[ dimension ] ) * _cellSize;
    }

    return pos;
}

//////////////////////////////////////////////////////////////////////
// Returns this cell's neighbours
//////////////////////////////////////////////////////////////////////
void ScalarField::cellNeighbours( const Tuple3i &index, IntArray &neighbours ) const
{
    neighbours.clear();

    for ( int dimension = 0; dimension < 3; dimension++ )
    {
        Tuple3i                  neighbourIndex = index;

        for ( int idx = -1; idx <= 1; idx += 2 )
        {
            neighbourIndex[ dimension ] = index[ dimension ] + idx;

            if ( neighbourIndex[ dimension ] >= 0
                    && neighbourIndex[ dimension ] < _divisions[ dimension ] )
            {
                neighbours.push_back( cellIndex( neighbourIndex ) );
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////
// Returns the indices of all cells whose values contribute
// to the given position, as well as trilinear interpolation
// coefficients for each cell.
//
// NOTE: This will fail if it lies outside of the bounds
//////////////////////////////////////////////////////////////////////
void ScalarField::interpolationCoefficients( const Vector3d &x,
                                             KeyValueArray &points ) const
{
    REAL                       normalization;
    REAL                       coefficient;
    Vector3d                   minPoint = x - _bbox.minBound();
    Vector3d                   maxPoint;
    Tuple3i                    minIndex;
    Tuple3i                    index;

    points.clear();

    normalization = _cellSize * _cellSize * _cellSize;

    for ( int dimension = 0; dimension < 3; dimension++ )
    {
        minPoint[ dimension ] -= 0.5 * _cellSize;

        minIndex[ dimension ] = (int)( minPoint[ dimension ] / _cellSize );

        TRACE_ASSERT( minIndex[ dimension ] >= 0
                && minIndex[ dimension ] < _divisions[ dimension ] - 1,
                "Position out of range" );

        minPoint[ dimension ] = _bbox.minBound()[ dimension ];
        minPoint[ dimension ] += _cellSize * ( 0.5 + (REAL)minIndex[ dimension ] );
        maxPoint[ dimension ] = minPoint[ dimension ] + _cellSize;
    }

    maxPoint -= x;
    minPoint = x - minPoint;

    // (0,0,0) point
    coefficient = maxPoint[ 0 ] * maxPoint[ 1 ] * maxPoint[ 2 ];
    coefficient /= normalization;

    index = minIndex;

    points.push_back( KeyValuePair( cellIndex( index ), coefficient ) );

    // (0,0,1) point
    coefficient = maxPoint[ 0 ] * maxPoint[ 1 ] * minPoint[ 2 ];
    coefficient /= normalization;

    index = minIndex;
    index[ 2 ] += 1;

    points.push_back( KeyValuePair( cellIndex( index ), coefficient ) );

    // (0,1,0) point
    coefficient = maxPoint[ 0 ] * minPoint[ 1 ] * maxPoint[ 2 ];
    coefficient /= normalization;

    index = minIndex;
    index[ 1 ] += 1;

    points.push_back( KeyValuePair( cellIndex( index ), coefficient ) );

    // (0,1,1) point
    coefficient = maxPoint[ 0 ] * minPoint[ 1 ] * minPoint[ 2 ];
    coefficient /= normalization;

    index = minIndex;
    index[ 1 ] += 1;
    index[ 2 ] += 1;

    points.push_back( KeyValuePair( cellIndex( index ), coefficient ) );

    // (1,0,0) point
    coefficient = minPoint[ 0 ] * maxPoint[ 1 ] * maxPoint[ 2 ];
    coefficient /= normalization;

    index = minIndex;
    index[ 0 ] += 1;

    points.push_back( KeyValuePair( cellIndex( index ), coefficient ) );

    // (1,0,1) point
    coefficient = minPoint[ 0 ] * maxPoint[ 1 ] * minPoint[ 2 ];
    coefficient /= normalization;

    index = minIndex;
    index[ 0 ] += 1;
    index[ 2 ] += 1;

    points.push_back( KeyValuePair( cellIndex( index ), coefficient ) );

    // (1,1,0) point
    coefficient = minPoint[ 0 ] * minPoint[ 1 ] * maxPoint[ 2 ];
    coefficient /= normalization;

    index = minIndex;
    index[ 0 ] += 1;
    index[ 1 ] += 1;

    points.push_back( KeyValuePair( cellIndex( index ), coefficient ) );

    // (1,1,1) point
    coefficient = minPoint[ 0 ] * minPoint[ 1 ] * minPoint[ 2 ];
    coefficient /= normalization;

    index = minIndex;
    index[ 0 ] += 1;
    index[ 1 ] += 1;
    index[ 2 ] += 1;

    points.push_back( KeyValuePair( cellIndex( index ), coefficient ) );
}

//////////////////////////////////////////////////////////////////////
// Interpolates a scalar field defined on this grid
//////////////////////////////////////////////////////////////////////
REAL ScalarField::interpolateScalarField( const Vector3d &x,
                                          const VECTOR &data ) const
{
    KeyValueArray              coefficients;
    REAL                       interpResult = 0.0;

    interpolationCoefficients( x, coefficients );

    for ( int i = 0; i < coefficients.size(); i++ )
    {
        const KeyValuePair      &coef = coefficients[ i ];

        TRACE_ASSERT( coef.first >= 0 && coef.first < data.size(),
                "Invalid interpolation coefficient" );

        interpResult += coef.second * data( coef.first );
    }

    return interpResult;
}

//////////////////////////////////////////////////////////////////////
// Interpolates a vector field in which each vector is the row of
// a matrix
//////////////////////////////////////////////////////////////////////
void ScalarField::interpolateVectorField( const Vector3d &x,
                                          const MATRIX &data,
                                          VECTOR &output ) const
{
    KeyValueArray              coefficients;

    interpolationCoefficients( x, coefficients );

    if ( output.size() != data.cols() )
    {
        output.resizeAndWipe( data.cols() );
    }
    else
    {
        output.clear();
    }

    for ( int i = 0; i < coefficients.size(); i++ )
    {
        const KeyValuePair      &coef = coefficients[ i ];

        TRACE_ASSERT( coef.first >= 0 && coef.first < data.rows(),
                "Invalid interpolation coefficient" );

        MATRIX::axpy( output.data(),
                // The desired row in data
                data.data() + coef.first * data.cols(),
                1, data.cols(), /* Copy size */
                coef.second /* multiplier */ );
    }
}
