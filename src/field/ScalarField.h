//////////////////////////////////////////////////////////////////////
// ScalarField.h: Interface for the ScalarField class
//
//////////////////////////////////////////////////////////////////////

#ifndef SCALAR_FIELD_H
#define SCALAR_FIELD_H

#include <linearalgebra/MATRIX.h>
#include <linearalgebra/VECTOR.h>

#include <geometry/BoundingBox.h>

#include <TYPES.h>

//////////////////////////////////////////////////////////////////////
// ScalarField class
//
// Models a scalar field on a regular grid
//////////////////////////////////////////////////////////////////////
class ScalarField {
    public:
        ScalarField();

        // Constructs scalar field within the given bounding box,
        // with the given division size.  The box may be enlarged
        // to ensure square cells.
        ScalarField( const BoundingBox &bbox, REAL cellSize );

        // Builds a scalar field starting at the given minimum bound,
        // and with the specified number of cell divisions and cell size
        ScalarField( const Vector3d &minBound, const Tuple3i &divisions,
                     REAL cellSize );

        // Destructor
        virtual ~ScalarField();

        // Returns the flat cell index for the given x, y, z index
        int                      cellIndex( const Tuple3i &index ) const;

        // Returns the x, y and z index for the given flat cell index
        Tuple3i                  cellIndex( int flatIndex ) const;

        // Cell center for the given index
        Vector3d                 cellPosition( const Tuple3i &index ) const;

        inline Vector3d          cellPosition( int flatIndex ) const
        {
            return cellPosition( cellIndex( flatIndex ) );
        }

        // Returns this cell's neighbours
        void                     cellNeighbours( const Tuple3i &index,
                                                 IntArray &neighbours ) const;

        inline void              cellNeighbours( int flatIndex, 
                                                 IntArray &neighbours ) const
        {
            cellNeighbours( cellIndex( flatIndex ), neighbours );
        }

        // Returns the indices of all cells whose values contribute
        // to the given position, as well as trilinear interpolation
        // coefficients for each cell.
        //
        // NOTE: This will fail if it lies outside of the bounds
        void                     interpolationCoefficients( const Vector3d &x,
                                                            KeyValueArray &points ) const;

        // Interpolates a scalar field defined on this grid
        REAL                     interpolateScalarField( const Vector3d &x,
                                                         const VECTOR &data ) const;

        // Interpolates a vector field in which each vector is the row of
        // a matrix
        void                     interpolateVectorField( const Vector3d &x,
                                                         const MATRIX &data,
                                                         VECTOR &output ) const;

        const BoundingBox       &bbox() const
        {
            return _bbox;
        }

        inline int               numCells() const
        {
            return _divisions[ 0 ] * _divisions[ 1 ] * _divisions[ 2 ];
        }

        inline REAL              cellSize() const
        {
            return _cellSize;
        }

        inline const Tuple3i    &cellDivisions() const
        {
            return _divisions;
        }

        inline bool              isBoundaryCell( const Tuple3i &index ) const
        {
            return (    index[ 0 ] == 0
                     || index[ 0 ] == _divisions[ 0 ] - 1
                     || index[ 1 ] == 0
                     || index[ 1 ] == _divisions[ 1 ] - 1 
                     || index[ 2 ] == 0
                     || index[ 2 ] == _divisions[ 2 ] - 1 );
        }

        inline bool              isBoundaryCell( int index ) const
        {
            return isBoundaryCell( cellIndex( index ) );
        }

    public:
        static const int         NUM_NEIGHBOURS;

    protected:

    private:
        BoundingBox              _bbox;
        Tuple3i                  _divisions;
        REAL                     _cellSize;

        ScalarArray3D           *_fieldValues;

};


#endif
