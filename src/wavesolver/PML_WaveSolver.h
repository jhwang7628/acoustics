//////////////////////////////////////////////////////////////////////
// PML_WaveSolver.h: Interface for the PML_WaveSolver class
//
//////////////////////////////////////////////////////////////////////

#ifndef PML_WAVE_SOLVER_H
#define PML_WAVE_SOLVER_H

#include <distancefield/distanceField.h>

#include <field/ScalarField.h>

#include <math/MLSModeInterpolator.hpp>
#include <geometry/BoundingBox.h>
#include <geometry/TriangleMesh.hpp>

#include <linearalgebra/MATRIX.h>
#include <linearalgebra/VECTOR.h>
#include <linearalgebra/Vector3.hpp>

#include <utils/Evaluator.h>
#include <utils/timer.hpp>

#include "MAC_Grid.h"
#include "WaveSolver.h"
#include <wavesolver/PML_WaveSolver_Settings.h> 

#include <TYPES.h>

#include <boost/function.hpp>

//////////////////////////////////////////////////////////////////////
// PML_WaveSolver class
//
// Solves the wave equation exterior to an object with given boundary
// conditions.  This version uses a staggered MAC grid discretization
// of the pressure and particle velocity equations.
//
//////////////////////////////////////////////////////////////////////
class PML_WaveSolver : public Solver 
{
    public:
        typedef boost::function<void (const vector<vector<FloatArray> >&w)> WriteCallback;

    private:
        REAL                     _waveSpeed;
        REAL                     _density; 

        // Discretization for the domain
        MAC_Grid                 _grid;

        // Pretty self-explanatory
        const REAL               _cellSize; 
        int                      _subSteps;
        REAL                     _endTime;
        REAL                     _timeStep;
        REAL                     _currentTime;
        int                      _timeIndex;

        // Number of acceleration fields to solve for
        int                      _N;

        // Leapfrog variables for a MAC grid
        // PML requires a separate pressure for each direction
        MATRIX                   _v[ 3 ];
        MATRIX                   _p[ 3 ];
        MATRIX                   _pFull;
        FloatArray               _pGhostCells[ 3 ]; 
        FloatArray               _pGhostCellsFull;

        MATRIX                   _pLastTimestep;      // for restarting
        MATRIX                   _pThisTimestep;      // for restarting
        MATRIX                   _vThisTimestep[ 3 ]; // for restarting

        MATRIX                   _pLaplacian; 
        MATRIX                   _pCollocated[3]; // rotating buffer for pressure: [old, current, new]
        int                      _pCollocatedInd; // point to new

        Timer<false>             _gradientTimer;
        Timer<false>             _divergenceTimer;
        Timer<false>             _stepTimer;
        Timer<false>             _algebraTimer;
        Timer<false>             _memoryTimer;
        Timer<false>             _writeTimer;
        Timer<false>             _ghostCellTimer;
        Timer<false>             _cellClassifyTimer;
        Timer<false>             _freshCellTimer;

        // Optionally write a 2D slice out of the finite difference grid
        int                      _zSlice;
        MATRIX                   _sliceData;
        VECTOR                   _listenerOutput; // not in use
        const Vector3Array      *_listeningPositions; // not in use
        const char              *_outputFile; // not in use
        WriteCallback           *_callback; // not in use
        std::vector<std::vector<FloatArray> > _waveOutput;

        // if on then ghost cell boundary treatment, otherwise rasterized
        bool                     _useGhostCellBoundary; 

        // volumetric pressure point source
        ExternalSourceEvaluator *_sourceEvaluator; 

        // objects in the scene 
        std::shared_ptr<FDTD_Objects> _objects; 
        PML_WaveSolver_Settings_Ptr _waveSolverSettings;

    public: 
        PML_WaveSolver() 
            : _cellSize(0.0), _listeningPositions(nullptr), _outputFile(nullptr)
        {}


        // Provide the size of the domain (bbox), finite difference division
        // size, and a signed distance function for the interior boundary.
        PML_WaveSolver( REAL timeStep,
                        const BoundingBox &bbox, REAL cellSize,
                        const TriMesh &mesh,
                        const DistanceField &distanceField,
                        REAL waveSpeed, 
                        REAL density, 
                        REAL distanceTolerance = 0.0,
                        bool useBoundary = true,
                        const Vector3Array *listeningPositions = NULL,
                        const char *outputFile = NULL,
                        WriteCallback *callback = NULL,
                        int subSteps = 1,
                        int N = 1,
                        REAL endTime = -1.0
                        );

        // Provide the size of the domain (bbox), finite difference division
        // size, and a list of SDFs for the interior boundary
        PML_WaveSolver( REAL timeStep,
                        const BoundingBox &bbox, REAL cellSize,
                        std::vector<const TriMesh *> &meshes,
                        std::vector<const DistanceField *> &boundaryFields,
                        REAL waveSpeed, 
                        REAL density,
                        REAL distanceTolerance = 0.0,
                        bool useBoundary = true,
                        const Vector3Array *listeningPositions = NULL,
                        const char *outputFile = NULL,
                        WriteCallback *callback = NULL,
                        int subSteps = 1,
                        int N = 1,
                        REAL endTime = -1.0
                        );

        // initialize from settings parsed from xml
        PML_WaveSolver(PML_WaveSolver_Settings_Ptr settings, std::shared_ptr<FDTD_Objects> objects);

        // to prevent repeated lines in constructor.
        void Reinitialize_PML_WaveSolver(const bool &useBoundary, const REAL &startTime); 

        // Destructor
        virtual ~PML_WaveSolver(){};

        inline bool GetGhostCellBoundary() { return _useGhostCellBoundary; } 
        inline MAC_Grid &GetGrid() {return _grid;}
        inline void setPMLBoundaryWidth( REAL width, REAL strength ){ _grid.setPMLBoundaryWidth( width, strength ); }
        inline void setZSlice( int slice ) { _zSlice = slice; }
        inline void SetExternalSource(ExternalSourceEvaluator *sourceEvaluator){ _sourceEvaluator = sourceEvaluator; } 

        virtual inline int numCells() const { return _grid.numPressureCells(); }
        virtual inline int N() const { return _N; }
        virtual inline REAL fieldDiameter() const { return _grid.fieldDiameter(); }
        virtual inline REAL currentSimTime() const { return _timeStep * (REAL)_timeIndex; }
        virtual inline Vector3d fieldPosition(const Tuple3i &index) const { return _grid.pressureFieldPosition( index ); }
        virtual inline Vector3d fieldPosition(int index) const { return _grid.pressureFieldPosition( index ); }
        virtual inline Vector3d velocityFieldPosition(const Tuple3i &index, const int &dim) const { return _grid.velocityFieldPosition( index, dim ); }
        virtual inline Vector3d velocityFieldPosition(int index, const int &dim) const { return _grid.velocityFieldPosition( index, dim ); }
        virtual inline Vector3d sceneCenter() const { return _grid.pressureField().bbox().center(); }
        virtual inline const Tuple3i &fieldDivisions() const { return _grid.pressureFieldDivisions(); }
        virtual inline const Tuple3i &velocityFieldDivisions(const int &dim) const { return _grid.velocityFieldDivisions(dim); }
        virtual inline const Vector3Array *listeningPositions() const { return _listeningPositions; }
        virtual inline const vector<const TriMesh *> &meshes() const { return _grid.meshes(); }

        int numVelocityCells(const int &dim) const;
        void SetGhostCellBoundary(const bool &isOn);

        void initSystem( REAL startTime );
        // Initialize the field data using non-zero initial conditions. 
        void initSystemNontrivial( const REAL startTime, const InitialConditionEvaluator * ic_eval ); 

        // fetch the pressure data using interpolation
        void FetchScalarData(const MATRIX &scalar, const ScalarField &field, const Vector3Array &listeningPoints, Eigen::MatrixXd &data); 
        void FetchPressureData(const Vector3Array &listeningPoints, Eigen::MatrixXd &data, const int dim=-1);
        void FetchVelocityData(const Vector3Array &listeningPoints, const int &dimension, Eigen::MatrixXd &data);
        void FetchPressureCellType(const Vector3Array &listeningPoints, Eigen::MatrixXd &data);
        void FetchCell(const int &cellIndex, MAC_Grid::Cell &cell) const; 
        void SampleAxisAlignedSlice(const int &dim, const REAL &offset, std::vector<MAC_Grid::Cell> &sampledCells) const; 
        void GetSolverDomain(Vector3d &minBound, Vector3d &maxBound) const;

        // Takes a single time step
        virtual bool stepSystem(const BoundaryEvaluator &bcEvaluator);
        virtual bool stepSystem();
        virtual bool stepSystemHalf(const int &flag);
        // Takes a single time step with restarting steps controlled by
        // N_restart. internally, smoothing is done using weighted average
        // method described in the paper: 
        // L.F.Shampine, Stability of the leapfrog/midpoint method
        //
        // The idea is if timeindex % N_restart == 0, then average operation is performed according 
        // to equation 8 in the paper
        //
        // p_i <- 1/4 * p_{i-1} + 1/2 * p_i + 1/4 * p_{i+1}
        //
        // TODO this currently produce bad results, maybe need to smooth velocity field as well
        //
        virtual bool stepSystemWithRestart(const int &N_restart); 

        // Get vertex pressure for each field
        virtual void vertexPressure( const Tuple3i &index, VECTOR &pressure ) const;
        virtual void vertexVelocity( const Tuple3i &index, const int &dim, VECTOR &velocity ) const;
        virtual void writeWaveOutput() const;

        //// debugging/testing methods ////
        REAL GetMaxCFL();
        void PrintAllFieldExtremum();

    private:
        void stepLeapfrog();
        void stepCollocated();

    friend std::ostream &operator <<(std::ostream &os, const PML_WaveSolver &solver); 
};

#endif
