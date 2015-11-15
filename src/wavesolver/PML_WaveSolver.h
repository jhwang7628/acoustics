//////////////////////////////////////////////////////////////////////
// PML_WaveSolver.h: Interface for the PML_WaveSolver class
//
//////////////////////////////////////////////////////////////////////

#ifndef PML_WAVE_SOLVER_H
#define PML_WAVE_SOLVER_H

#include <distancefield/distanceField.h>

#include <field/ScalarField.h>

#include <geometry/BoundingBox.h>
#include <geometry/TriangleMesh.hpp>

#include <linearalgebra/MATRIX.h>
#include <linearalgebra/VECTOR.h>
#include <linearalgebra/Vector3.hpp>

#include <utils/Evaluator.h>
#include <utils/timer.hpp>

#include "MAC_Grid.h"
#include "WaveSolver.h"

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
class PML_WaveSolver : public Solver {
    public:
        typedef boost::function<void (const vector<vector<FloatArray> >&w)> WriteCallback;
        typedef boost::function<void (const vector<REAL>&w, const REAL &t, const int &i)> WriteCallbackIndividual;

        // Provide the size of the domain (bbox), finite difference division
        // size, and a signed distance function for the interior boundary.
        // and the write function pointer.
        PML_WaveSolver( REAL timeStep,
                        const BoundingBox &bbox, REAL cellSize,
                        const TriMesh &mesh,
                        const DistanceField &distanceField,
                        REAL distanceTolerance = 0.0,
                        bool useBoundary = true,
                        const Vector3Array *listeningPositions = NULL,
                        const char *outputFile = NULL,
                        WriteCallback *callback = NULL,
                        WriteCallbackIndividual *callbacki = NULL,
                        int subSteps = 1,
                        int N = 1,
                        REAL endTime = -1.0 );

        // Provide the size of the domain (bbox), finite difference division
        // size, and a signed distance function for the interior boundary.
        PML_WaveSolver( REAL timeStep,
                        const BoundingBox &bbox, REAL cellSize,
                        const TriMesh &mesh,
                        const DistanceField &distanceField,
                        WaveSolverPointData * rawData,
                        REAL distanceTolerance = 0.0,
                        bool useBoundary = true,
                        const Vector3Array *listeningPositions = NULL,
                        const char *outputFile = NULL,
                        WriteCallback *callback = NULL,
                        int subSteps = 1,
                        int N = 1,
                        REAL endTime = -1.0 );


        // Provide the size of the domain (bbox), finite difference division
        // size, and a list of SDFs for the interior boundary
        PML_WaveSolver( REAL timeStep,
                        const BoundingBox &bbox, REAL cellSize,
                        std::vector<const TriMesh *> &meshes,
                        std::vector<const DistanceField *> &boundaryFields,
                        REAL distanceTolerance = 0.0,
                        bool useBoundary = true,
                        const Vector3Array *listeningPositions = NULL,
                        const char *outputFile = NULL,
                        WriteCallback *callback = NULL,
                        WriteCallbackIndividual *callbacki = NULL,
                        int subSteps = 1,
                        int N = 1,
                        REAL endTime = -1.0 );


        // Destructor
        virtual ~PML_WaveSolver();

        void setPMLBoundaryWidth( REAL width, REAL strength );
        void setHarmonicSource( const HarmonicSourceEvaluator * hs_eval ); 

        void SetWaveSolverPointDataPosition(); 
        void AddCurrentPressureToWaveSolverPointData(); 


        // Time steps the system in the given interval
        void solveSystem( REAL startTime, REAL endTime,
                const BoundaryEvaluator &bcEvaluator );

        void initSystem( REAL startTime );

        void initSystemNontrivial( const REAL startTime, const InitialConditionEvaluator * ic_eval ); 

        // Takes a single time step
        virtual bool stepSystem( const BoundaryEvaluator &bcEvaluator );
        // Takes a single time step with source function
        virtual bool stepSystem( const BoundaryEvaluator &bcEvaluator, const HarmonicSourceEvaluator *hsEval );

        // Get vertex pressure for each field
        virtual void vertexPressure( const Tuple3i &index, VECTOR &pressure );

        virtual void writeWaveOutput() const;

        virtual  const Tuple3i  &fieldDivisions() const
        {
            return _grid.pressureFieldDivisions();
        }

        // Test the max CFL condition based on particle velocity in the grid. 
        virtual REAL GetMaxCFL();
        


        virtual Vector3d         fieldPosition( const Tuple3i &index ) const
        {
            return _grid.pressureFieldPosition( index );
        }

        virtual Vector3d         fieldPosition( int index ) const
        {
            return _grid.pressureFieldPosition( index );
        }

#if 0
        inline const TriMesh    &mesh() const
        {
            return _grid.mesh();
        }
#endif
        virtual const vector<const TriMesh *> &meshes() const
        {
            return _grid.meshes();
        }

        virtual int              numCells() const
        {
            return _grid.numPressureCells();
        }

        virtual const Vector3Array *listeningPositions() const
        {
            return _listeningPositions;
        }

        virtual REAL             fieldDiameter() const
        {
            return _grid.fieldDiameter();
        }

        virtual  int             N() const
        {
            return _N;
        }

        virtual REAL             currentSimTime() const
        {
            return _timeStep * (REAL)_timeIndex;
        }

        virtual Vector3d         sceneCenter() const
        {
            return _grid.pressureField().bbox().center();
        }

        void                     setZSlice( int slice )
        {
            _zSlice = slice;
        }

    protected:

    private:
        void                     stepLeapfrog( const BoundaryEvaluator &bcEvaluator );
        void                     stepLeapfrog( const BoundaryEvaluator &bcEvaluator, const HarmonicSourceEvaluator *hsEval );

    private:
        static constexpr REAL        WAVE_SPEED = 343.0;

        // Discretization for the domain
        MAC_Grid                 _grid;

        // Pretty self-explanatory
        REAL                     _timeStep;
        REAL                     _currentTime;
        int                      _timeIndex;
        int                      _subSteps;

        // Leapfrog variables for a MAC grid
        //
        // PML requires a separate pressure for each direction
        MATRIX                   _v[ 3 ];
        MATRIX                   _p[ 3 ];
        MATRIX                   _pFull;

        // Number of fields to solve for
        int                      _N;

        REAL                     _endTime;

        const REAL               _cellSize; 
        const Vector3Array      *_listeningPositions;
        const char              *_outputFile;

        std::vector<std::vector<FloatArray> >
            _waveOutput;

        // cells indexed by hsIndex is source cell that emits harmonic waves 
        std::vector<int> hsIndex; 

        VECTOR                   _listenerOutput;

        WriteCallback           *_callback;
        WriteCallbackIndividual *_callbackInd;
        WaveSolverPointData     *_rawData;

        Timer<false>             _gradientTimer;
        Timer<false>             _divergenceTimer;
        Timer<false>             _stepTimer;
        Timer<false>             _algebraTimer;
        Timer<false>             _memoryTimer;
        Timer<false>             _writeTimer;

        // Optionally write a 2D slice out of the finite difference grid
        int                      _zSlice;
        MATRIX                   _sliceData;

};

#endif
