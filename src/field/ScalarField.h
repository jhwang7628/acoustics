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
// 
// If wants to iterate the cell indices, the order corresponds to the 
// memory layout is x -> y -> z. 
// FOr more details, please see function cellIndex
//
// Example (get cell position according to memory location): 
//
//  for (x=0; x<Nx; ++x) 
//    for (y=0; y<Ny; ++y) 
//      for (z=0; z<Nz; ++z) 
//      {
//         cellPosition(x,y,z)
//      }
//////////////////////////////////////////////////////////////////////
class ScalarField {
    public:
        static const int         NUM_NEIGHBOURS;
        struct RangeIndices
        {
            Vector3i startIndex; 
            Vector3i endIndex; 
            Vector3i dimensionIteration;
        }; 
        enum BoundaryType 
        {
            Not_Boundary = 0x00,
            Positive_X_Boundary = 0x01,
            Negative_X_Boundary = 0x02, 
            Positive_Y_Boundary = 0x04, 
            Negative_Y_Boundary = 0x08, 
            Positive_Z_Boundary = 0x10, 
            Negative_Z_Boundary = 0x20
        };

    private:
        BoundingBox              _bbox;
        Tuple3i                  _divisions;
        Tuple3i                  _indexOffset; 
        REAL                     _cellSize;

        ScalarArray3D           *_fieldValues;

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
        inline int               cellIndex( const int &x, const int &y, const int &z) const 
        {
            Tuple3i index(x,y,z); 
            return cellIndex(index); 
        }
        // Returns the x, y and z index for the given flat cell index
        Tuple3i                  cellIndex( int flatIndex ) const;

        inline int               cellIndex(const int &d0, const int &d1, const int &d2,
                                           const int &i0, const int &i1, const int &i2) const
        {
            Tuple3i indices; 
            indices[d0] = i0; 
            indices[d1] = i1; 
            indices[d2] = i2; 
            return cellIndex(indices); 
        }

        // Cell center for the given index
        Vector3d                 cellPosition( const Tuple3i &index ) const;
        void                     cellPosition( const Tuple3i &index, Vector3d &pos ) const;

        inline Vector3d          cellPosition( int flatIndex ) const
        {
            return cellPosition( cellIndex( flatIndex ) );
        }
        Tuple3i                  enclosingCell(const Vector3d &position) const; 

        // Returns this cell's neighbours
        void                     cellNeighbours( const Tuple3i &index,
                                                 IntArray &neighbours ) const;
        void                     cellNeighbours( const int &flatIndex, 
                                                 IntArray &neighbours,
                                                 IntArray &neighbourTopology) const; 

        inline void              cellNeighbours( int flatIndex, 
                                                 IntArray &neighbours ) const
        {
            cellNeighbours( cellIndex( flatIndex ), neighbours );
        }

        // all 26 of cell's neighbours
        void                     cell26Neighbours( const int &flatIndex, IntArray &neighbours ) const; 

        // add the neighbours on the eight orthogonal corners to the current cell. 
        void                     AddCornerNeighbours( const int &flatIndex, IntArray &neighbours ) const; 

        // find the enclosing (8) neighbours for a given point in space
        void                     enclosingNeighbours(const Vector3d &position, IntArray &neighbours) const; 

        // given query box, compute lower index and number of indices in each dimension. 
        // For example, if (0,0,0) and (1,1,1) is given, then this function should return 
        // the first index that is within this boundingbox, and (Nx,Ny,Nz) for
        // iteration in x,y,z directions. 
        void                     GetIterationBox(const Vector3d &minBound, const Vector3d &maxBound, RangeIndices &indices); 

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

        inline bool              isBoundaryCell( const Tuple3i &index, const int offset = 0 ) const
        {
            return (    index[ 0 ] == offset
                     || index[ 0 ] == _divisions[ 0 ] - 1 - offset
                     || index[ 1 ] == offset
                     || index[ 1 ] == _divisions[ 1 ] - 1 - offset
                     || index[ 2 ] == offset
                     || index[ 2 ] == _divisions[ 2 ] - 1 - offset);
        }

        inline bool              isBoundaryCell( int index, const int offset = 0 ) const
        {
            return isBoundaryCell( cellIndex( index ), offset );
        }

        inline unsigned char     boundaryCellType(int index) const 
        {
            const Tuple3i indices = cellIndex(index); 
            unsigned char type = 0;
            if      (indices[0] == 0              ) type = (type|Negative_X_Boundary); 
            else if (indices[0] == _divisions[0]-1) type = (type|Positive_X_Boundary); 
            if      (indices[1] == 0              ) type = (type|Negative_Y_Boundary); 
            else if (indices[1] == _divisions[1]-1) type = (type|Positive_Y_Boundary); 
            if      (indices[2] == 0              ) type = (type|Negative_Z_Boundary); 
            else if (indices[2] == _divisions[2]-1) type = (type|Positive_Z_Boundary); 
            return type; 
        }

        inline Tuple3i          boundaryCellNormal(int index, const int offset) const 
        {
            const Tuple3i indices = cellIndex(index); 
            Tuple3i type; 
            if      (indices[0] == offset                ) type[0] = -1; 
            else if (indices[0] == _divisions[0]-1-offset) type[0] = +1; 
            if      (indices[1] == offset                ) type[1] = -1; 
            else if (indices[1] == _divisions[1]-1-offset) type[1] = +1; 
            if      (indices[2] == offset                ) type[2] = -1; 
            else if (indices[2] == _divisions[2]-1-offset) type[2] = +1; 
            return type; 
        }

        inline Tuple3i          boundaryCellNormal(int index) const 
        {
            return boundaryCellNormal(index,0); 
        }

        inline Vector3d minBound() const {return cellPosition(0);}
        inline Vector3d maxBound() const {return cellPosition(numCells()-1);}
        inline const Tuple3i &indexOffset() const {return _indexOffset;}

        void MoveCenter(const Tuple3i &amount); 

        //// debug methods ////
        void TestSubindices();

    friend std::ostream &operator <<(std::ostream &os, const ScalarField &field);

};


#endif
