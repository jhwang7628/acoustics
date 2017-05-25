//////////////////////////////////////////////////////////////////////
// ScalarField.cpp: Implementation of the ScalarField class
//
//////////////////////////////////////////////////////////////////////

#include "ScalarField.h"

#include <utils/IO.h>
#include <config.h> 

const int ScalarField::NUM_NEIGHBOURS = 6;

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
ScalarField::ScalarField()
    : _bbox( Vector3d( 0, 0, 0 ), Vector3d( 0, 0, 0 ) ),
      _divisions( 1, 1, 1 ),
      _indexOffset(0,0,0),
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
    : _indexOffset(0,0,0),
      _cellSize( cellSize ),
      _fieldValues( NULL )
{
    Vector3d                   minBound;
    Vector3d                   maxBound;

    for ( int dimension = 0; dimension < 3; dimension++ )
    {
        _divisions[ dimension ]
            = (int)ceil( (bbox.axislength( dimension )-EPS) / cellSize ); // minus EPS to prevent rounding error in case division gives integer

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
ScalarField::ScalarField( const Vector3d &minBound, const Tuple3i &divisions, REAL cellSize )
    : _divisions( divisions ),
      _indexOffset(0,0,0),
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
    // with offset
    const int newIndex[3] = {(index[0] + _indexOffset[0]) % _divisions[0], 
                             (index[1] + _indexOffset[1]) % _divisions[1],
                             (index[2] + _indexOffset[2]) % _divisions[2]}; 
    flatIndex = newIndex[0] * ( _divisions[1] * _divisions[2] );
    flatIndex += newIndex[ 1 ] * _divisions[ 2 ];
    flatIndex += newIndex[ 2 ];
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

    // with offset
    for (int d=0; d<3; ++d)
        index[d] = (index[d] - _indexOffset[d] + _divisions[d]) % _divisions[d]; 
    return index;
}

//////////////////////////////////////////////////////////////////////
// Cell center for the given index
//////////////////////////////////////////////////////////////////////
Vector3d ScalarField::cellPosition( const Tuple3i &index ) const
{
    Vector3d                 pos;
    cellPosition(index,pos);
    return pos;
}

void ScalarField::cellPosition( const Tuple3i &index, Vector3d &pos ) const
{
    //Vector3d                 pos;

    for ( int dimension = 0; dimension < 3; ++dimension )
    {
        pos[dimension] = _bbox.minBound()[ dimension ];
        pos[dimension] += ( 0.5 + (REAL)index[ dimension ] ) * _cellSize;
    }
}

//////////////////////////////////////////////////////////////////////
// Returns the cell that encloses the position
//////////////////////////////////////////////////////////////////////
Tuple3i ScalarField::enclosingCell(const Vector3d &position) const
{
    const Vector3d diff = position - _bbox.minBound(); 
    Tuple3i indices; 
    indices.x = (int) (diff.x/_cellSize); 
    indices.y = (int) (diff.y/_cellSize); 
    indices.z = (int) (diff.z/_cellSize); 
    if (indices.x < 0 || indices.x >= _divisions[0]) 
        throw std::runtime_error("**ERROR** enclosingCell out of bounds"); 
    if (indices.y < 0 || indices.y >= _divisions[1]) 
        throw std::runtime_error("**ERROR** enclosingCell out of bounds"); 
    if (indices.z < 0 || indices.z >= _divisions[2]) 
        throw std::runtime_error("**ERROR** enclosingCell out of bounds"); 
    return indices; 
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

void ScalarField::AddCornerNeighbours( const int &flatIndex, IntArray &neighbours ) const
{
    const Tuple3i index = cellIndex(flatIndex); 

    neighbours.push_back(cellIndex(index[0]-1, index[1]-1, index[1]-1)); 
    neighbours.push_back(cellIndex(index[0]-1, index[1]-1, index[1]+1)); 
    neighbours.push_back(cellIndex(index[0]-1, index[1]+1, index[1]-1)); 
    neighbours.push_back(cellIndex(index[0]-1, index[1]+1, index[1]+1)); 
    neighbours.push_back(cellIndex(index[0]+1, index[1]+1, index[1]-1)); 
    neighbours.push_back(cellIndex(index[0]+1, index[1]+1, index[1]+1)); 
}


void ScalarField::cell26Neighbours( const int &flatIndex, IntArray &neighbours) const 
{
    neighbours.clear();

    const Tuple3i index = cellIndex(flatIndex); 

    Tuple3i neighbourIndex = index;
    for (int kk=-1; kk<=1; kk++) 
        for (int jj=-1; jj<=1; jj++) 
            for (int ii=-1; ii<=1; ii++) 
            {
                neighbourIndex[0] = index[0]+ii;
                neighbourIndex[1] = index[1]+jj;
                neighbourIndex[2] = index[2]+kk;

                if (neighbourIndex[0] < 0 || neighbourIndex[0] >= _divisions[0] || 
                    neighbourIndex[1] < 0 || neighbourIndex[1] >= _divisions[1] ||
                    neighbourIndex[2] < 0 || neighbourIndex[2] >= _divisions[2] ||
                    (ii==0 && jj==0 && kk==0)                                     )
                    continue; 

                neighbours.push_back(cellIndex(neighbourIndex)); 
            }
}

void ScalarField::enclosingNeighbours(const Vector3d &position, IntArray &neighbours) const 
{

    neighbours.clear(); 

    Vector3d diff = position - _bbox.minBound(); 
    diff -= _cellSize/2.0; 

    // compute the lower corner and clamp it to the (boundary - 1) cell
    const int x = min( max( (int) (diff.x/_cellSize), 1), _divisions[0]-2 ); 
    const int y = min( max( (int) (diff.y/_cellSize), 1), _divisions[1]-2 ); 
    const int z = min( max( (int) (diff.z/_cellSize), 1), _divisions[2]-2 ); 

    neighbours.push_back( cellIndex(x+0,y+0,z+0) ); 
    neighbours.push_back( cellIndex(x+0,y+0,z+1) ); 
    neighbours.push_back( cellIndex(x+0,y+1,z+0) ); 
    neighbours.push_back( cellIndex(x+0,y+1,z+1) ); 
    neighbours.push_back( cellIndex(x+1,y+0,z+0) ); 
    neighbours.push_back( cellIndex(x+1,y+0,z+1) ); 
    neighbours.push_back( cellIndex(x+1,y+1,z+0) ); 
    neighbours.push_back( cellIndex(x+1,y+1,z+1) ); 

}

void ScalarField::GetIterationBox(const Vector3d &in_minBound, const Vector3d &in_maxBound, RangeIndices &indices)
{
    Vector3d minBound = in_minBound; 
    Vector3d maxBound = in_maxBound; 

    if (!_bbox.isInside(in_minBound))
    {
        minBound.x = std::max<REAL>(_bbox.minBound().x, minBound.x);
        minBound.y = std::max<REAL>(_bbox.minBound().y, minBound.y);
        minBound.z = std::max<REAL>(_bbox.minBound().z, minBound.z);
    }
    if (!_bbox.isInside(in_maxBound))
    {
        maxBound.x = std::min<REAL>(_bbox.maxBound().x, maxBound.x);
        maxBound.y = std::min<REAL>(_bbox.maxBound().y, maxBound.y);
        maxBound.z = std::min<REAL>(_bbox.maxBound().z, maxBound.z);
    }
    //assert(_bbox.isInside(minBound) && _bbox.isInside(maxBound)); 
    const Vector3d &fieldMinBound = _bbox.minBound();
    for (int d=0; d<3; ++d) 
    {
        indices.startIndex[d] = (int)std::floor((minBound[d]-fieldMinBound[d])/_cellSize); 
        indices.endIndex[d]   = (int) std::ceil((maxBound[d]-fieldMinBound[d])/_cellSize);
        indices.dimensionIteration[d] = indices.endIndex[d] - indices.startIndex[d]; 
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

    for ( size_t i = 0; i < coefficients.size(); i++ )
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

    for ( size_t i = 0; i < coefficients.size(); i++ )
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

////////////////////////////////////////////////////////////////////////////////
// Move the center towards direction and amount indicated by amount.
// This will rotate the way cellIndex is computed, effectively shift
// all data. For example, if amount=(1,0,0), the scene center moves 
// one unit towards positive x, therefore the data should rotate 
// one unit towards negative x (or move offset towards positive x).
////////////////////////////////////////////////////////////////////////////////
void ScalarField::MoveCenter(const Tuple3i &amount)
{
    for (int d=0; d<3; ++d) 
        _indexOffset[d] = (_indexOffset[d] + amount[d] + _divisions[d]) % _divisions[d]; 
    Vector3d offset(((double)amount[0])*_cellSize,
                    ((double)amount[1])*_cellSize,
                    ((double)amount[2])*_cellSize); 
    _bbox.setMinBound(_bbox.minBound() + offset); 
    _bbox.setMaxBound(_bbox.maxBound() + offset); 
}

void ScalarField::TestSubindices()
{
    std::cout << *this << std::endl;
    const Vector3d testMin(0.39,0.40,0.41); 
    const Vector3d testMax(0.79,0.80,1.0000); 
    RangeIndices indices; 
    GetIterationBox(testMin, testMax, indices); 
    COUT_SDUMP(indices.startIndex); 
    COUT_SDUMP(indices.dimensionIteration); 
    COUT_SDUMP(testMin); 
    COUT_SDUMP(testMax);

    for (int ii=0; ii<indices.dimensionIteration.x; ++ii) 
        for (int jj=0; jj<indices.dimensionIteration.y; ++jj) 
            for (int kk=0; kk<indices.dimensionIteration.z; ++kk) 
            {
                const int ind_x = indices.startIndex.x + ii; 
                const int ind_y = indices.startIndex.y + jj; 
                const int ind_z = indices.startIndex.z + kk; 
                const Tuple3i index(ind_x, ind_y, ind_z); 
                const Vector3d position = cellPosition(index); 
                COUT_SDUMP(position);
            }
}

std::ostream &operator <<(std::ostream &os, const ScalarField &field)
{
    os << "--------------------------------------------------------------------------------\n" 
        << "Class ScalarField\n" 
        << "--------------------------------------------------------------------------------\n";
    os << "min bound: " << field._bbox.minBound() << "\n"
        << "max bound: " << field._bbox.maxBound() << "\n"
        << "divisions: " << field._divisions << "\n"
        << "cell size: " << field._cellSize << "\n"; 
    os << "--------------------------------------------------------------------------------" 
        << std::flush; 
    return os; 
}

