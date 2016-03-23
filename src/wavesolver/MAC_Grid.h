//////////////////////////////////////////////////////////////////////
// MAC_Grid.h: Interface for the MAC_Grid class
//
//////////////////////////////////////////////////////////////////////

#ifndef MAC_GRID_H
#define MAC_GRID_H

#include <distancefield/distanceField.h>
#include <distancefield/closestPointField.h>

#include <field/ScalarField.h>

#include <linearalgebra/MATRIX.h>
#include <linearalgebra/VECTOR.h>
#include <linearalgebra/Vector3.hpp>

#include <geometry/BoundingBox.h>
#include <geometry/TriangleMesh.hpp>

#include <utils/Evaluator.h>

#include <TYPES.h>

//////////////////////////////////////////////////////////////////////
// MAC_Grid class
//
// Builds a finite difference discretization of a MAC_Grid, handling
// an irregular internal boundary defined by an SDF
//
// This construction is based on [Marella, S et al.] "Sharp interface
// Cartesian grid method I: An easily implemented technique for 3D
// moving boundary computations"
//////////////////////////////////////////////////////////////////////
class MAC_Grid {
    private:
        typedef TriangleMesh<REAL>  TriMesh;

    public:
        // Provide the size of the domain, finite difference division
        // size, and a signed distance function for the interior boundary
        MAC_Grid( const BoundingBox &bbox, REAL cellSize,
                  const TriMesh &mesh,
                  const DistanceField &distanceField,
                  REAL distanceTolerance = 0.0,
                  int N = 1 );

        MAC_Grid( const BoundingBox &bbox, REAL cellSize,
                  std::vector<const TriMesh *> &meshes,
                  std::vector<const DistanceField *> &boundaryFields,
                  REAL distanceTolerance = 0.0,
                  int N = 1 );

        // Destructor
        virtual ~MAC_Grid();

        // Initializes field using a rasterized version of the boundary
        void initFieldRasterized( bool useBoundary );

        // Width is the number of cells we wish to absorb in
        void setPMLBoundaryWidth( REAL width, REAL strength );

        // Gets the derivative for the given velocity field, based on the
        // input pressure value
        void velocityDerivative( const MATRIX &p, const BoundaryEvaluator &bc,
                                 MATRIX &dV_dt, int dimension,
                                 REAL t, REAL alpha = 1.0 ) const;

        // Gets the derivative of the pressure field using the current
        // velocity field values in each direction
        void pressureDerivative( const MATRIX *v[ 3 ], MATRIX &p,
                                 REAL alpha = 1.0 ) const;

        // Performs a velocity update in the given direction, as detailed
        // by Liu et al. (equation (14))
        void PML_velocityUpdate( const MATRIX &p, const BoundaryEvaluator &bc,
                                 MATRIX &v, int dimension,
                                 REAL t, REAL timeStep, REAL density = 1.0 );

        // Performs a pressure update for the given pressure direction,
        // as detailed by Liu et al. (equation (16))
        void PML_pressureUpdate( const MATRIX &v, MATRIX &p, int dimension,
                                 REAL timeStep, REAL c, 
                                 const ExternalSourceEvaluator *sourceEvaluator, const REAL simulationTime,
                                 REAL density = 1.0
                                 );

        // Samples data from a z slice of the finite difference grid and
        // puts it in to a matrix
        void sampleZSlice( int slice, const MATRIX &p, MATRIX &sliceData );

        void SmoothFieldInplace(MATRIX &p1, MATRIX &p2, MATRIX &p3, REAL w1, REAL w2, REAL w3);

        // set the PML effect region flag
        inline void SetCornellBoxBoundaryCondition(const bool &isOn)
        { 
            _cornellBoxBoundaryCondition = isOn; 
        } 

        inline int               numPressureCells() const
        {
            return _pressureField.numCells();
        }

        inline int               numVelocityCellsX() const
        {
            return _velocityField[ 0 ].numCells();
        }

        inline int               numVelocityCellsY() const
        {
            return _velocityField[ 1 ].numCells();
        }

        inline int               numVelocityCellsZ() const
        {
            return _velocityField[ 2 ].numCells();
        }

        inline const Tuple3i    &pressureFieldDivisions() const
        {
            return _pressureField.cellDivisions();
        }

        inline Vector3d          pressureFieldPosition( const Tuple3i &index ) const
        {
            return _pressureField.cellPosition( index );
        }

        inline Vector3d          pressureFieldPosition( int index ) const
        {
            return _pressureField.cellPosition( index );
        }

        inline int               pressureFieldVertexIndex( const Tuple3i &index ) const
        {
            return _pressureField.cellIndex( index );
        }

        inline Tuple3i           pressureFieldVertexIndex( int index ) const
        {
            return _pressureField.cellIndex( index );
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

#if 0
        inline const IntArray   &interfacialCells() const
        {
            return _interfacialCells;
        }
#endif

        const ScalarField       &pressureField() const
        {
            return _pressureField;
        }

        inline REAL              fieldDiameter() const
        {
            return _pressureField.bbox().maxlength();
        }

        inline int               N() const
        {
            return _N;
        }

    protected:

    private:
        // Classifies cells as either a bulk cell, ghost cell, or
        // interfacial cell
        void classifyCells( bool useBoundary );


        // Returns the absorption coefficient along a certain
        // dimension for a point in space.
        //
        // FIXME: For now, we will use a quadratic profile here, though
        // we may need to try something more complex later on.
        inline REAL PML_absorptionCoefficient( const Vector3d &x, REAL absorptionWidth,
                                               int dimension );

        // Scaling factor for the initial velocity update in each time step
        //
        // We need the absorption coefficient in this dimension
        //
        // Equation (15) in Liu et al. (f_{1,\nu})
        static inline REAL PML_velocityUpdateCoefficient( REAL absorptionCoefficient,
                REAL timeStep )
        {
            REAL                   coefficient;

            coefficient = 1.0 / timeStep - absorptionCoefficient / 2.0;
            coefficient /= 1.0 / timeStep + absorptionCoefficient / 2.0;

            return coefficient;
        }

        // Scaling factor for the pressure gradient update applied to the
        // velocity in each time step
        //
        // Equation (15) in Liu et al. (f_{2,\nu})
        static inline REAL PML_pressureGradientCoefficient(
                REAL absorptionCoefficient,
                REAL timeStep,
                REAL cellSize,
                REAL density )
        {
            REAL                   coefficient;

            coefficient = 1.0 / timeStep + absorptionCoefficient / 2.0;
            coefficient *= density * cellSize;

            return ( -1.0 / coefficient );
        }

        // "Directional coefficient" D_{\nu,v} used in equations (17), (18)
        // of Liu et al
        //
        // Note, we have no damping factor (ie. gamma = 0)
        static inline REAL PML_directionalCoefficient( REAL absorptionCoefficient,
                REAL timeStep )
        {
            return ( 1.0 / timeStep + absorptionCoefficient / 2.0 );
        }

        // Pressure scaling for pressure update
        //
        // f_{3,\nu} in Liu et al.
        static inline REAL PML_pressureUpdateCoefficient( REAL absorptionCoefficient,
                REAL timeStep,
                REAL directionalCoefficient )
        {
            return ( 1.0 / timeStep - absorptionCoefficient / 2.0 )
                / directionalCoefficient;
        }

        // Scaling for the velocity divergence term in the pressure update
        //
        // f_{5,\nu} in Liu et al.
        static inline REAL PML_divergenceCoefficient( REAL density, REAL c,
                REAL directionalCoefficient )
        {
            return ( -1.0 * density * c * c / directionalCoefficient );
        }

        

    private:
        std::vector<const DistanceField *>   _boundaryFields;
        std::vector<const TriMesh *>         _boundaryMeshes;

        REAL                     _distanceTolerance;

        ScalarField              _pressureField;
        ScalarField              _velocityField[ 3 ];

        // isBulkCell and isGhostCell refer to cells in the pressure grid
        BoolArray                _isBulkCell;
        BoolArray                _isGhostCell;

        // isInterfacialCell refers to cells in the three velocity grid
        BoolArray                _isVelocityBulkCell[ 3 ];

        BoolArray                _isVelocityInterfacialCell[ 3 ];

        IntArray                 _bulkCells;
        IntArray                 _ghostCells;

        IntArray                 _velocityBulkCells[ 3 ];

        IntArray                 _velocityInterfacialCells[ 3 ];

        IntArray                 _interfacialBoundaryIDs[ 3 ];
        FloatArray               _interfacialBoundaryDirections[ 3 ];
        FloatArray               _interfacialBoundaryCoefficients[ 3 ];

        // Dimensionality of the data we are working with
        int                      _N;

        REAL                     _PML_absorptionWidth;
        REAL                     _PML_absorptionStrength;

        bool                     _cornellBoxBoundaryCondition; // see readme in the solver

};

#endif
