//////////////////////////////////////////////////////////////////////
// Laplacian.h: Interface for the Laplacian class
//
//////////////////////////////////////////////////////////////////////

#ifndef LAPLACIAN_H
#define LAPLACIAN_H

#include <distancefield/distanceField.h>

#include <field/ScalarField.h>

#include <geometry/BoundingBox.h>
#include <geometry/TriangleMesh.hpp>

#include <linearalgebra/MATRIX.h>
#include <linearalgebra/VECTOR.h>
#include <linearalgebra/Vector3.hpp>

#include <utils/Evaluator.h>

#include <TYPES.h>

//////////////////////////////////////////////////////////////////////
// Laplacian class
//
// Builds a finite difference discretization of a Laplacian, handling
// an irregular internal boundary defined by an SDF
//
// This construction is based on [Marella, S et al.] "Sharp interface
// Cartesian grid method I: An easily implemented technique for 3D
// moving boundary computations"
//////////////////////////////////////////////////////////////////////
class Laplacian {
    private:
        typedef TriangleMesh<REAL>  TriMesh;

    public:
        // Provide the size of the domain, finite difference division
        // size, and a signed distance function for the interior boundary
        Laplacian( const BoundingBox &bbox, REAL cellSize,
                   const TriMesh &mesh,
                   const DistanceField &distanceField,
                   REAL distanceTolerance = 0.0,
                   int N = 1 );

        Laplacian( const BoundingBox &bbox, REAL cellSize,
                   std::vector<const TriMesh *> &meshes,
                   std::vector<const DistanceField *> &boundaryFields,
                   REAL distanceTolerance = 0.0,
                   int N = 1 );

        // Destructor
        virtual ~Laplacian();

        // Initializes the field boundary condition, and produces
        // a Laplacian discretization
        void initField( bool useBoundary );

        // Initializes field using a rasterized version of the boundary
        void initFieldRasterized( bool useBoundary );

        // Applies the laplacian to field p, and places the result in Lp.
        //
        // alpha is an optional multiplier to apply to the Laplacian result
        void apply( const MATRIX &p, MATRIX &Lp, REAL alpha = 1.0 ) const;

        // Forms the boundary condition contribution to Lp resulting from
        // the application of the Laplacian
        //
        // alpha is an optional multiplier to apply to the Laplacian result
        void applyBoundary( MATRIX &Lp, const BoundaryEvaluator &bc,
                REAL t, REAL alpha = 1.0 );

#if 0
        void constructLaplacianMatrix( SPARSE_MATRIX &L );
#endif

        inline int               numCells() const
        {
            return _field.numCells();
        }

        inline const Tuple3i    &fieldDivisions() const
        {
            return _field.cellDivisions();
        }

        inline Vector3d          fieldPosition( const Tuple3i &index ) const
        {
            return _field.cellPosition( index );
        }

        inline Vector3d          fieldPosition( int index ) const
        {
            return _field.cellPosition( index );
        }

        inline int               fieldVertexIndex( const Tuple3i &index ) const
        {
            return _field.cellIndex( index );
        }

        inline Tuple3i           fieldVertexIndex( int index ) const
        {
            return _field.cellIndex( index );
        }

#if 0
        inline const TriMesh    &mesh() const
        {
            return *_boundaryMeshes[ 0 ];
        }
#endif
        inline const vector<const TriMesh *> &meshes() const
        {
            return _boundaryMeshes;
        }

        inline const IntArray   &ghostCells() const
        {
            return _ghostCells;
        }

        inline const IntArray   &interfacialCells() const
        {
            return _interfacialCells;
        }

        const ScalarField       &field() const
        {
            return _field;
        }

        inline REAL              fieldDiameter() const
        {
            return _field.bbox().maxlength();
        }

        inline int               N() const
        {
            return _N;
        }

    protected:

    private:
        // Model a point on the boundary, which will be represented
        // in terms of some other points in the domain, as well as
        // a boundary contribution (just given by a coefficient).
        struct BoundaryPoint {
            BoundaryPoint()
                : _boundaryCoefficient( 0.0 )
            {
            }

            BoundaryPoint( const Vector3d &x, const Vector3d &n, REAL coefficient )
                : _x( x ),
                _n( n ),
                _boundaryCoefficient( coefficient )
            {
            }

            // Position and normal (pointing in to the solution domain)
            Vector3d               _x;
            Vector3d               _n;

            REAL                   _boundaryCoefficient;
        };

        typedef std::vector<BoundaryPoint> BoundaryPointArray;

        // Classifies cells as either a bulk cell, ghost cell, or
        // interfacial cell
        void classifyCells( bool useBoundary );

        // Checks the validity of ghost cells.  A ghost cell is only
        // valid if we can extrapolate along its gradient to the boundary
        // and then out a distance of 2 * h, and lie in a region governed
        // only by ghost, bulk or interfacial cells
        //
        // If we find any invalid cells, we will add new ghost cells
        void checkGhostCellValidity();

        // Builds all ghost cell coefficients
        void buildGhostCellCoefficients();

        // Builds the coefficient list for a single ghost cell.
        //
        // NOTE: Currently this just uses the distance field and
        // not the mesh itself.
        void buildGhostCellCoefficients( int cell_idx );

        // Given an interfacial point, builds coefficients for the second
        // derivative in one direction
        void buildDirectionalCoefficients(
                int cell_idx, GridDirection direction,
                KeyValueArray &coefficients,
                BoundaryPointArray &boundaryCoefficients );

        // Builds the coefficients resulting from the line joining cell_idx
        // and neighbour_idx and the boundary sdf.
        //
        // Returns the ratio of the distance between the intersection point
        // and cell_idx and the grid cell spacing
        REAL buildIntersectionCoefficients(
                int cell_idx, int neighbour_idx,
                KeyValueArray &coefficients,
                BoundaryPointArray &boundaryCoefficients );

        // Given an intersection point with the boundary, build coefficients
        // by extrapolating a normal line out from the boundary and constraining
        // the point at the boundary to obey a Neumann BC
        void buildIntersectionCoefficients(
                const Vector3d &intersection,
                KeyValueArray &coefficients,
                BoundaryPointArray &boundaryCoefficients );

        bool insideDomain( int cell_idx )
        {
            return ( _isBulkCell[ cell_idx ] || _isGhostCell[ cell_idx ]
                    || _isInterfacialCell[ cell_idx ] );
        }

        // Returns the distance along vector d, originating at x at which
        // an intersection with the boundary field occurs.  The
        // result is assumed to lie within (0,1).
        //
        // x is assumed to lie outside of the boundary field, and x + d is
        // assumed to lie inside.
        //
        // Intersection is determined via bisection.
        REAL boundaryIntersection( Vector3d x, Vector3d d, REAL tolerance = 1e-4 );

        // Builds coefficients for all interface points
        void buildInterfacialCoefficients();

        // Builds coefficients for interface points assuming a rasterized
        // version of the input shape
        void buildInterfacialCoefficientsRasterized();

    private:
        // Ghost cell coefficient calculations.  Loosely based on
        // equation (41) of [Marelli et al.]
        //
        // TODO: Look in to these calculations a bit more
        static REAL ghostBoundaryCoefficient( REAL d1, REAL d2, REAL h );
        static REAL ghost_pG_Coefficient( REAL d1, REAL d2 );
        static REAL ghost_pG2_Coefficient( REAL d1, REAL d2 );
        static REAL ghost_pG3_Coefficient( REAL d1, REAL d2 );

    private:
        ScalarField              _field;
        REAL                     _distanceTolerance;

        std::vector<const TriMesh *>         _boundaryMeshes;
        std::vector<const DistanceField *>   _boundaryFields;


#if 0
        const DistanceField     *_boundarySDF;
        const TriMesh           *_mesh;
#endif


        BoolArray                _isBulkCell;
        BoolArray                _isInterfacialCell;
        BoolArray                _isGhostCell;

        IntArray                 _bulkCells;
        IntArray                 _interfacialCells;
        IntArray                 _interfacialBoundaryIDs;
        IntArray                 _ghostCells;

        // Dimensionality of the data we are working with
        int                      _N;

        // Coefficient lists for ghost points
        std::vector<KeyValueArray>         _ghostCoefficients;
        std::vector<BoundaryPointArray>    _ghostBoundaryCoefficients;

        // For each cell, store a pointer in to the array above, if
        // the cell is a ghost cell
        IntArray                           _ghostCoefficientPointers;

        // Coefficient lists for interfacial cells, and a pointer array
        // in to that list
        std::vector<KeyValueArray>         _interfacialCoefficients;
        std::vector<BoundaryPointArray>    _interfacialBoundaryCoefficients;

#if 0
        IntArray                           _interfactialCellPointers;
#endif

};

#endif
