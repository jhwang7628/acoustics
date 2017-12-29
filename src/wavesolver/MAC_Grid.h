//////////////////////////////////////////////////////////////////////
// MAC_Grid.h: Interface for the MAC_Grid class
//
//////////////////////////////////////////////////////////////////////

#ifndef MAC_GRID_H
#define MAC_GRID_H

#include <unordered_map>
#include <set>
#include <map>
#include <distancefield/distanceField.h>
#include <distancefield/closestPointField.h>

#include <field/ScalarField.h>

#include <linearalgebra/MATRIX.h>
#include <linearalgebra/VECTOR.h>
#include <linearalgebra/Vector3.hpp>

#include <geometry/BoundingBox.h>
#include <geometry/TriangleMesh.hpp>

#include <wavesolver/PML_WaveSolver_Settings.h> 
#include <wavesolver/FDTD_Objects.h>

#include <utils/Evaluator.h>
#include <utils/timer.hpp>

#include <Eigen/Dense> 
#include <Eigen/Sparse>

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
class MAC_Grid 
{
    public: 
        typedef int GhostCellId; 

        // only used for classification with multiple threads
        struct GhostCellInfo
        {
            int dim; 
            int cellIndex; 
            int boundaryObject;
            REAL boundaryDirection; 
            REAL boundaryCoefficient; 
            int childArrayPosition; 
            int ghostCellParent; 
            Vector3d ghostCellPosition; 
            int ghostCellBoundaryID; 
        }; 

        class FVMetaData
        {
            public:
            std::unordered_map<int, std::shared_ptr<std::vector<TriangleIdentifier> > > cellMap; 
            void Clear(){cellMap.clear();}
        };

        //######################################################################
        // Struct MAC_Grid::GhostCell
        //  ownerCell flat cell index that owns this ghost cell, the owner cell
        //            is always a solid cell (but has fluid neighbour)
        //  neighbourCell flat cell index that share this ghost cell with owner,
        //                the neighbourCell is always a fluid cell
        //  topology topology (edge) between owner and neighbour. can only
        //           take values {+1, -1, +2, -2, +3, -3}
        //           +1: neighbour is at positive x direction relative to owner
        //           -1: neighbour is at negative x direction relative to owner
        //           +2: neighbour is at positive y direction relative to owner
        //           -2: neighbour is at negative y direction relative to owner
        //           +3: neighbour is at positive z direction relative to owner
        //           -3: neighbour is at negative z direction relative to owner
        //  pressure pressure value stored at this ghost cell
        //  tid identify last-seen closest triangle
        //  cache associated cached field
        //######################################################################
        struct GhostCell_Cache
        {
            TriangleIdentifier tid; 
        }; 
        //using GhostCell_Cache_UPtr = std::unique_ptr<GhostCell_Cache>; 
        using GhostCell_Cache_Ptr = std::shared_ptr<GhostCell_Cache>; 
        using GhostCell_Key  = unsigned long long int; 
        using BoundaryGhostCell_Key  = std::string; 
        struct GhostCell
        {
            int ownerCell; 
            int neighbourCell; 
            int topology;
            int type;  // 0: bulky, 1: constraint, 2: shell 
            int boundaryObject = -1; 
            REAL distance = std::numeric_limits<REAL>::max(); 
            Vector3d position; 
            REAL pressure; 
            GhostCell_Cache_Ptr cache; 
            static GhostCell_Key MakeKey(const int &_a, const int &_b)
            {
                // Szudzik's function
                const unsigned long long a = std::min(_a, _b); 
                const unsigned long long b = std::max(_a, _b); 
                return (a>=b ? a*a+a+b : a+b*b); 
            }
            inline GhostCell_Key MakeKey() const 
            {return GhostCell::MakeKey(ownerCell, neighbourCell);}

            bool boundary = false; 
            std::string ownerSolverId; 
            std::string neighbourSolverId;
            BoundaryInterface_Ptr interface; 
            static BoundaryGhostCell_Key MakeKey_Boundary(const std::string &oId, 
                                                          const std::string &nId,
                                                          const int &oCell, const int &nCell)
            {
                return (oId < nId ?
                             oId + ":" + std::to_string(oCell) + ";"
                           + nId + ":" + std::to_string(nCell)
                         : 
                             nId + ":" + std::to_string(nCell) + ";"
                           + oId + ":" + std::to_string(oCell)
                       );
            }
            BoundaryGhostCell_Key MakeKey_Boundary()
            {
                return GhostCell::MakeKey_Boundary(ownerSolverId, neighbourSolverId, 
                                                   ownerCell, neighbourCell); 
            }
            Vector3d CellNormal()
            {
                Vector3d n; 
                n[abs(topology)-1] = (topology>0 ? 1.0 : -1.0); 
                return n;
            }
        };
        using  GhostCell_UPtr = std::unique_ptr<GhostCell>;
        using  GhostCell_Ptr = std::shared_ptr<GhostCell>; 

        class GhostCell_Deprecated
        {
            public:
                static std::vector<Timer<false> > ghostCellTimers; // NOTE not thread safe
                struct BoundarySamples
                {
                    bool isBulk; 
                    int neighbour_idx; 
                    int gc_value_idx; // which subdivision child is this sample
                    Vector3d position; 
                    Vector3d normal;
                    BoundarySamples(){}
                    BoundarySamples(const Vector3d &pos, const Vector3d &nor, const int &nei, const int &gcind) 
                        : neighbour_idx(nei), gc_value_idx(gcind), position(pos), normal(nor){}
                };
                struct VolumeSamples
                {
                    bool isBulk; 
                    Vector3d position; 
                    VolumeSamples(){}
                    VolumeSamples(const Vector3d &pos)
                        : position(pos){}
                }; 
                static int valuePointer; // points to current value
                bool validSample; 
                int parent_idx; 
                REAL volume; 
                REAL dp_dn_dot_S; 
                FloatArray positions;
                std::vector<FloatArray> values; // last, current, future
                std::vector<BoundarySamples> boundarySamples; 
                std::vector<VolumeSamples> volumeSamples; 
                std::shared_ptr<std::vector<TriangleIdentifier> > hashedTriangles; 
                GhostCell_Deprecated(const int &parent, std::shared_ptr<std::vector<TriangleIdentifier> > &triangles)
                    : validSample(true), parent_idx(parent), positions(FloatArray(6)), values(6, FloatArray(6, 0.0)), hashedTriangles(triangles)
                {}
        };

        class Cell
        {
            public: 
                int         index = -1; 
                REAL        r_identity;
                std::string s_identity;
                FloatArray  gcValue;  // for ghost cell only
                Tuple3i     indices; 
                Vector3d    centroidPosition; 
                Vector3d    lowerCorner; 
                Vector3d    upperCorner; 
                REAL        pDirectional[3]; 
                REAL        pFull; 
                REAL        vx[2];
                REAL        vy[2];
                REAL        vz[2];
                REAL        h; // cell size
                REAL        laplacian; 
                std::shared_ptr<GhostCell> ghostCell; 

                inline REAL Divergence() const
                {
                    REAL div = 0.0; 
                    div += (vx[1] - vx[0]) / h; 
                    div += (vy[1] - vy[0]) / h; 
                    div += (vz[1] - vz[0]) / h; 
                    return div; 
                }
        }; 

    private:
        typedef TriangleMesh<REAL>  TriMesh;

        // used for jacobi iteration of the ghost cell values in pressure update
        struct JacobiIterationData 
        {
            int                 cellId; 
            int                 coupleCount = 0; // 0 if fully uncoupled
            std::vector<int>    nnzIndex; 
            std::vector<REAL>   nnzValue; 
            REAL                RHS; 
        };

        struct PML_PressureCell;
        struct PML_VelocityCell; 

        struct PML_PressureCell
        {
            int index; 
            int neighbour_v_left[3];
            int neighbour_v_right[3];
            REAL updateCoefficient[3]; 
            REAL divergenceCoefficient[3]; 
            REAL absorptionCoefficient; 
            Vector3d position;
            bool needRemove; 
            PML_VelocityCell *lNeighbours[3]; 
            PML_VelocityCell *rNeighbours[3]; 
        };

        struct PML_VelocityCell
        {
            int index; 
            int dimension; 
            int neighbourInterior; // 0: no; 1: right is interior; -1: left is interior
            int neighbour_p_left;
            int neighbour_p_right;
            REAL updateCoefficient; 
            REAL gradientCoefficient; 
            REAL velocity; 
            Vector3d position;
            bool needRemove;
        }; 

    private: 
        std::vector<const DistanceField *>   _boundaryFields;
        std::vector<const TriMesh *>         _boundaryMeshes;

        REAL                     _distanceTolerance;
        REAL                     _cellSize; 

        ScalarField              _pressureField;
        ScalarField              _velocityField[ 3 ];

        // isBulkCell and isGhostCell refer to cells in the pressure grid
        BoolArray                _isBulkCell;
        BoolArray                _isGhostCell;
        BoolArray                _isPMLCell;
        BoolArray                _isOneCellCavity;

        // store all the triangle ids in cell for shells (not doing this
        // for rigid bodies)
        std::vector<std::set<TriangleIdentifier, TIComp>> _cellTriangles; 
        // map key: shell object id, map value: grid cell index that collide with this shell object
        std::vector<int> _collidedCellsForShells_v; 

        // isInterfacialCell refers to cells in the three velocity grid
        // interfacial cells are classified as non bulk although its value is
        // valid
        BoolArray                _isVelocityBulkCell[ 3 ];
        BoolArray                _isVelocityInterfacialCell[ 3 ];

        std::vector<PML_PressureCell> _pmlPressureCells; 
        std::vector<PML_VelocityCell> _pmlVelocityCells; 
        std::unordered_map<int, std::shared_ptr<GhostCell_Deprecated> > _ghostCellsCollection; 
        std::unordered_map<GhostCell_Key, GhostCell_UPtr> _ghostCells; 
        std::unordered_map<int, std::set<GhostCell_Key>>  _ghostCellsTree; 
        std::unordered_map<BoundaryGhostCell_Key, GhostCell_UPtr> _boundaryGhostCells; 
        using GhostCellType         = std::unordered_map<GhostCell_Key, GhostCell_UPtr>; 
        using BoundaryGhostCellType = std::unordered_map<BoundaryGhostCell_Key, GhostCell_UPtr>; 
        std::unordered_map<GhostCell_Key, GhostCell_Cache_Ptr> _ghostCellsCached; 
        std::vector<IntArray>    _ghostCellsChildren; 
        BoolArray                _classified; // show if this cell has been classified, used in classifyCellsDynamic_FAST
        BoolArray                _pressureCellHasValidHistory; 
        BoolArray                _velocityCellHasValidHistory[3]; 

        std::map<int,int>        _ghostCellsInverse; 

        IntArray                 _velocityInterfacialCells[ 3 ];
        IntArray                 _interfacialBoundaryIDs[ 3 ];
        FloatArray               _interfacialBoundaryDirections[ 3 ];
        FloatArray               _interfacialBoundaryCoefficients[ 3 ];
        IntArray                 _containingObject;

        // for subdivided ghost cells
        IntArray                 _interfacialGhostCellID[ 3 ]; // map to ghost cell index 
        IntArray                 _ghostCellParents;
        std::vector<Vector3d>    _ghostCellPositions; 
        IntArray                 _ghostCellBoundaryIDs; 

        IntArray                 _ghostCellChildArrayPositions; 
        std::unordered_map<GhostCellId,std::unordered_map<int,TriangleIdentifier>> _ghostCellPreviousTriangles; 

        // Dimensionality of the data we are working with
        int                      _N;

        REAL                     _PML_absorptionWidth;
        REAL                     _PML_absorptionStrength;

        bool                     _cornellBoxBoundaryCondition; // see readme in the solver
        bool                     _useGhostCellBoundary; // see readme in the solver

        // for jacobi iteration on ghost-cell value solve
        std::vector<JacobiIterationData>    _ghostCellCoupledData; 

        // handles all the objects in the scene
        std::shared_ptr<FDTD_Objects>       _objects; 
        PML_WaveSolver_Settings_Ptr         _waveSolverSettings;

        // for finite-volume formulation
        FVMetaData _fvMetaData; 

        // boundary interface 
        std::list<BoundaryInterface_Ptr> _boundaryInterfaces; 

        bool _meshChanged = true; // Used to store if a mesh has changed. If so, do not use cached ghost cell data.
        int _meshChangedCounter = 0; // classifyCells_FAST needs to run twice after new obj is loaded and has intersection with boundary. this keeps track of that. 

    public:
        std::string *grid_id = nullptr; 

        MAC_Grid(){}
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

        // initialize the grid by a boundingbox and cellsize
        MAC_Grid(const BoundingBox &bbox, PML_WaveSolver_Settings_Ptr settings, std::shared_ptr<FDTD_Objects> objects);

        void Reinitialize_MAC_Grid(const BoundingBox &bbox, const REAL &cellSize); 

        // Destructor
        virtual ~MAC_Grid();

        // Initializes field using a rasterized version of the boundary
        void initFieldRasterized( bool useBoundary );

        // Width is the number of cells we wish to absorb in
        void setPMLBoundaryWidth( REAL width, REAL strength );

        void pressureFieldLaplacian(const MATRIX &value, MATRIX &laplacian) const; 
        void pressureFieldLaplacianGhostCell(const MATRIX &value, const FloatArray &ghostCellValue, MATRIX &laplacian) /* const */; 

        void setMeshChanged() {_meshChanged = true; _meshChangedCounter = 1;}

        // Performs a velocity update in the given direction, as detailed
        // by Liu et al. (equation (14))
        void PML_velocityUpdate( const MATRIX &p, const FloatArray &pGC, 
                                 MATRIX &v, int dimension,
                                 REAL t, REAL timeStep, REAL density);

        void PML_velocityUpdateCollocated(const REAL &simulationTime, const MATRIX (&pDirectional)[3], const MATRIX &pFull, MATRIX (&v)[3]); 
        void PML_pressureUpdateCollocated(const REAL &simulationTime, const MATRIX (&v)[3], MATRIX (&_pDirectional)[3], MATRIX &_pLast, MATRIX &_pCurr, MATRIX &_pNext, MATRIX &currLaplacian, MATRIX &pastLaplacian);

        // Performs a pressure update for the given pressure direction,
        // as detailed by Liu et al. (equation (16))
        void PML_pressureUpdate( const MATRIX &v, MATRIX &pDirectional, MATRIX &pFull, int dimension,
                                 REAL timeStep, REAL c, 
                                 const ExternalSourceEvaluator *sourceEvaluator, const REAL simulationTime,
                                 REAL density);

        // Performs a pressure update for the given pressure direction,
        // as detailed by Liu et al. (equation (16))
        void PML_pressureUpdateFull( const MATRIX *v, MATRIX &p,
                                               const REAL &timeStep, const REAL &c, 
                                               const ExternalSourceEvaluator *sourceEvaluator, const REAL &simulationTime,
                                               const REAL &density);

        void UpdatePMLPressure(MATRIX (&pDirectional)[3], MATRIX &pFull); 

        // Performs a pressure update for the ghost cells. 
        void PML_pressureUpdateGhostCells(MATRIX &p, FloatArray &pGC, const REAL &timeStep, const REAL &c, const REAL &simulationTime, const REAL density); 
        void PML_pressureUpdateGhostCells_Coupled(MATRIX &p, FloatArray &pGC, const REAL &timeStep, const REAL &c, const REAL &simulationTime, const REAL density); 
#ifdef USE_FV
        void UpdateGhostCells_FV(MATRIX &p, const REAL &simulationTime); 
#endif

        // Samples data from a z slice of the finite difference grid and
        // puts it in to a matrix
        void sampleZSlice( int slice, const MATRIX &p, MATRIX &sliceData );
        void SampleAxisAlignedSlice(const int &dim, const REAL &offset, MATRIX const (&pDirectional)[3], const MATRIX &pFull, const FloatArray &pGC, const MATRIX (&v)[3], std::vector<Cell> &sampledCells) const; 

        // Smooth field given weights
        void SmoothFieldInplace(MATRIX &p1, MATRIX &p2, MATRIX &p3, REAL w1, REAL w2, REAL w3);

        // Compute the ghost cell index from the flatten index 
        void ComputeGhostCellInverseMap();

        void FreshCellInterpolate(MATRIX &p, const REAL &simulationTime, const REAL &density); 

        // interpolate the fresh cell value
        // deal with the problem for rasterized boundary, when
        // solid->interfacial for velocity cells, or when solid->bulk for
        // pressure cells. In these two situations the cell does not have valid
        // history. 
        // see 2008 Mittal paper section 2.2.3
        void InterpolateFreshPressureCell(MATRIX &p, const REAL &timeStep, const REAL &simulationTime, const REAL &density);
        void InterpolateFreshVelocityCell(MATRIX &v, const int &dimension, const REAL &timeStep, const REAL &simulationTime);

        // set the PML effect region flag
        inline void SetGhostCellBoundary(const bool &isOn) { _useGhostCellBoundary = isOn; }
        inline int N() const { return _N; }
        inline int numPressureCells() const { return _pressureField.numCells(); }
        inline int numVelocityCellsX() const { return _velocityField[ 0 ].numCells(); }
        inline int numVelocityCellsY() const { return _velocityField[ 1 ].numCells(); }
        inline int numVelocityCellsZ() const { return _velocityField[ 2 ].numCells(); }
        inline REAL fieldDiameter() const { return _pressureField.bbox().maxlength(); }
        inline BoundingBox PressureBoundingBox() const {return _pressureField.bbox();}
        inline BoundingBox VelocityBoundingBox(const int &dim) const {return _velocityField[dim].bbox();}
        inline const ScalarField &pressureField() const { return _pressureField; }
        inline ScalarField &pressureField() { return _pressureField; }
        inline const ScalarField &velocityField(const int &ind) const { return _velocityField[ind]; }
        inline ScalarField &velocityField(const int &ind) { return _velocityField[ind]; }
        inline Vector3d pressureFieldPosition(const Tuple3i &index) const { return _pressureField.cellPosition( index ); }
        inline Vector3d pressureFieldPosition(int index) const { return _pressureField.cellPosition( index ); }
        inline Vector3d velocityFieldPosition(const Tuple3i &index, const int &dim) const { return _velocityField[dim].cellPosition( index ); }
        inline Vector3d velocityFieldPosition(int index, const int &dim) const { return _velocityField[dim].cellPosition( index ); }
        inline int pressureFieldVertexIndex(const Tuple3i &index) const { return _pressureField.cellIndex( index ); }
        inline Tuple3i pressureFieldVertexIndex( int index ) const { return _pressureField.cellIndex( index ); }
        inline int velocityFieldVertexIndex(const Tuple3i &index, const int &dim) const { return _velocityField[dim].cellIndex( index ); }
        inline Tuple3i velocityFieldVertexIndex( int index, const int &dim ) const { return _velocityField[dim].cellIndex( index ); }
        inline const Tuple3i &pressureFieldDivisions() const { return _pressureField.cellDivisions(); }
        inline const Tuple3i &velocityFieldDivisions(const int &dim) const { return _velocityField[dim].cellDivisions(); }
        //inline const IntArray &ghostCells() const { return _ghostCells; }
        inline const vector<const TriMesh *> &meshes() const { return _boundaryMeshes; }
        inline const bool IsVelocityCellSolid(const int &cell_idx, const int &dim) { return !_isVelocityInterfacialCell[dim].at(cell_idx) && !_isVelocityBulkCell[dim].at(cell_idx); }
        inline const bool IsPressureCellGhost(const int &cell_idx) const {return _isGhostCell.at(cell_idx);}
        inline const bool IsPressureCellBulk(const int &cell_idx) const {return _isBulkCell.at(cell_idx);}
        inline const bool IsPressureCellSolid(const int &cell_idx) const {return !_isBulkCell.at(cell_idx) && !_isGhostCell.at(cell_idx);}
        inline const FVMetaData &GetFVMetaData(){return _fvMetaData;}
        inline const auto &GetGhostCells(){return _ghostCells;}
        inline const auto &GetBoundaryGhostCells(){return _boundaryGhostCells;}
        inline void ClearGhostCellPreviousTriangles(){_ghostCellPreviousTriangles.clear();}

        void classifyCells_FAST(MATRIX (&pCollocated)[3], const bool &verbose=false); 
        void classifyCellsDynamic(MATRIX &pFull, MATRIX (&p)[3], FloatArray &pGCFull, FloatArray (&pGC)[3], MATRIX (&v)[3], const bool &useBoundary, const bool &verbose=false);
        void classifyCellsDynamic_FAST(MATRIX &pFull, MATRIX (&p)[3], FloatArray &pGCFull, FloatArray (&pGC)[3], MATRIX (&v)[3], const bool &useBoundary, const bool &verbose=false);
        void classifyCellsFV(MATRIX &pFull, MATRIX (&p)[3], FloatArray &pGCFull, FloatArray (&pGC)[3], MATRIX (&v)[3], const bool &useBoundary, const bool &verbose=false);
        void ComputeGhostCellSolveResidual(const FloatArray &p, REAL &minResidual, REAL &maxResidual, int &maxResidualEntry, REAL &maxOffDiagonalEntry); 
        REAL PressureCellType(const int &idx, const BoundingBox *sceneBox=nullptr) const;
        void ResetCellHistory(const bool &valid); 
        void ResetClassified(){std::fill(_classified.begin(), _classified.end(), false);}
        void GetCell(const int &cellIndex, MATRIX const (&pDirectional)[3], const MATRIX &pFull, const FloatArray &pGC, const MATRIX (&v)[3], Cell &cell) const; 
        void SetClassifiedSubset(const ScalarField &field, const int &N, const std::vector<ScalarField::RangeIndices> &indices, const bool &state);
        void CheckClassified(); 
        void Push_Back_GhostCellInfo(const int &gcIndex, const GhostCellInfo &info, FloatArray &pGCFull, FloatArray (&pGC)[3]); 
        std::unique_ptr<GhostCell> MakeGhostCell(const int cell, const int neighbour, const int topology, const int ghostcellType); 
        int InPressureCell(const Vector3d &position); 
        void RemoveOldPML(const BoundingBox &sceneBox); 
        void UpdatePMLAbsorptionCoeffs(const BoundingBox &sceneBox); 
        void UpdatePML(const BoundingBox &sceneBox);
        void FillBoundaryFreshCellGrid(const int &dim, const int &ind, MATRIX &pCurr, const MATRIX &pLast); 
        // if sign>0, grab all +dimension face cell indices and positions
        // if sign<0, grab all -dimension ...
        void GetAllBoundaryCells(const int &dimension, const int &sign, std::vector<int> &indices, std::vector<Vector3d> &positions); 
        // return pressure at cell_b for this solver (solver_b)
        bool BoundaryGhostCellPressure(const std::string &solver_a, const int &cell_a, const int &cell_b, REAL &pressure_b) const; 

        //// debug methods //// 
        void PrintFieldExtremum(const MATRIX &field, const std::string &fieldName); 
        void PrintGhostCellTreeInfo(); 
        void visualizeClassifiedCells(); 
        bool ExamineJacobiMatrix(); 
        REAL EstimateEnergy(const MATRIX &pCurr, const MATRIX &pLast);
        void ToFile_GhostCellCoupledMatrix(const std::string &filename); 
        void ToFile_GhostCellLinearSystem(const char *filename); 

        inline void ClearBoundaryInterface()
        {_boundaryInterfaces.clear();}
        inline void AddBoundaryInterface(BoundaryInterface_Ptr interface)
        {_boundaryInterfaces.push_back(interface);} 

    private:
        // Classifies cells as either a bulk cell, ghost cell, or
        // interfacial cell
        void classifyCells( bool useBoundary );

        bool InsidePML(const Vector3d &x, const REAL &absorptionWidth, int &flag, const BoundingBox *sceneBox=nullptr); 

        // Returns the absorption coefficient along a certain
        // dimension for a point in space.
        inline REAL PML_absorptionCoefficient(const Vector3d &x, REAL absorptionWidth, int dimension, 
                                              const BoundingBox *sceneBox=nullptr) const;

        // Scaling factor for the initial velocity update in each time step
        //
        // We need the absorption coefficient in this dimension
        //
        // Equation (15) in Liu et al. (f_{1,\nu})
        static inline REAL PML_velocityUpdateCoefficient( REAL absorptionCoefficient, REAL timeStep )
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
        static inline REAL PML_pressureGradientCoefficient( REAL absorptionCoefficient, REAL timeStep, REAL cellSize, REAL density )
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
        static inline REAL PML_directionalCoefficient( REAL absorptionCoefficient, REAL timeStep )
        {
            return ( 1.0 / timeStep + absorptionCoefficient / 2.0 );
        }

        // Pressure scaling for pressure update
        //
        // f_{3,\nu} in Liu et al.
        static inline REAL PML_pressureUpdateCoefficient( REAL absorptionCoefficient, REAL timeStep, REAL directionalCoefficient )
        {
            return ( 1.0 / timeStep - absorptionCoefficient / 2.0 ) / directionalCoefficient;
        }

        // Scaling for the velocity divergence term in the pressure update
        //
        // f_{5,\nu} in Liu et al.
        static inline REAL PML_divergenceCoefficient( REAL density, REAL c, REAL directionalCoefficient )
        {
            return ( -1.0 * density * c * c / directionalCoefficient );
        }


        // find image point for the ghost-cell method
        inline void FindImagePoint(const Vector3d &cellPosition, const int &boundaryObjectID, Vector3d &closestPoint, Vector3d &imagePoint, Vector3d &erectedNormal); 

        // fill the Vandermonde matrix 
        inline void FillVandermondeRegularS(const int &row, const Vector3d &cellPosition, Eigen::MatrixXd &V, const Vector3d &origin, const REAL &h);
        inline void FillVandermondeRegular (const Vector3d &cellPosition, Eigen::VectorXd &V);
        inline void FillVandermondeRegular (const int &row, const Vector3d &cellPosition, Eigen::MatrixXd &V);
        inline void FillVandermondeBoundary(const int &row, const Vector3d &boundaryPosition, const Vector3d &boundaryNormal, Eigen::MatrixXd &V);
        inline void FillVandermondeBoundaryS(const int &row, const Vector3d &boundaryPosition, const Vector3d &boundaryNormal, Eigen::MatrixXd &V, const Vector3d &origin, const REAL &h);
        void ClearUnusedCache(); 

    friend std::ostream &operator <<(std::ostream &os, const MAC_Grid &grid); 
    friend std::ostream &operator <<(std::ostream &os, const Cell &cell); 
};

#endif
