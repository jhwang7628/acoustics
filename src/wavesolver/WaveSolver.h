//////////////////////////////////////////////////////////////////////
// WaveSolver.h: Interface for the WaveSolver class
//
//////////////////////////////////////////////////////////////////////

#ifndef WAVE_SOLVER_H
#define WAVE_SOLVER_H

#include <distancefield/distanceField.h>

#include <field/ScalarField.h>

#include <linearalgebra/MATRIX.h>
#include <linearalgebra/VECTOR.h>

#include <geometry/BoundingBox.h>
#include <geometry/TriangleMesh.hpp>

#include <utils/Evaluator.h>
#include <utils/timer.hpp>

#include "Laplacian.h"

#include <TYPES.h>

#include "WaveSolverPointData.h"

#include <boost/function.hpp>

//////////////////////////////////////////////////////////////////////
// Solver interface
//
//////////////////////////////////////////////////////////////////////
class Solver {
    public:
        typedef TriangleMesh<REAL>  TriMesh;

    public:
        Solver()
        {
        }

        virtual ~Solver()
        {
        }

        virtual const Tuple3i &fieldDivisions() const = 0;
        virtual int N() const = 0;
        virtual const vector<const TriMesh *> &meshes() const = 0;
        virtual bool stepSystem( const BoundaryEvaluator &bcEvaluator ) = 0;
        virtual void writeWaveOutput() const = 0;
        virtual Vector3d fieldPosition( const Tuple3i &index ) const = 0;
        virtual Vector3d fieldPosition( int index ) const = 0;
        virtual const Vector3Array *listeningPositions() const = 0;
        virtual REAL fieldDiameter() const = 0;
        virtual void vertexPressure( const Tuple3i &index, VECTOR &pressure ) = 0;
        virtual int numCells() const = 0;
        virtual REAL currentSimTime() const = 0;
        virtual Vector3d sceneCenter() const = 0;

};

//////////////////////////////////////////////////////////////////////
// WaveSolver class
//
// Solves the wave equation exterior to an object with given boundary
// conditions.
//
// TODO: Eventually should implement perfectly matched layers, or
//       some form of damping.
//////////////////////////////////////////////////////////////////////
class WaveSolver : public Solver {
    public:
        typedef boost::function<void (const vector<vector<FloatArray> >&w)>
            WriteCallback;
        typedef boost::function<void (const vector<REAL>&w, const REAL &t)> WriteCallbackIndividual;

        // Provide the size of the domain (bbox), finite difference division
        // size, and a signed distance function for the interior boundary.
        WaveSolver( REAL timeStep,
                    const BoundingBox &bbox, REAL cellSize,
                    const TriMesh &mesh,
                    const DistanceField &distanceField,
                    bool useLeapfrog,
                    REAL distanceTolerance = 0.0,
                    bool useBoundary = true,
                    bool rasterize = false,
                    const Vector3Array *listeningPositions = NULL,
                    const char *outputFile = NULL,
                    WriteCallback *callback = NULL,
                    int subSteps = 1,
                    int N = 1 );

        // Provide the size of the domain (bbox), finite difference division
        // size, and a list of SDFs for the interior boundary
        WaveSolver( REAL timeStep,
                const BoundingBox &bbox, REAL cellSize,
                std::vector<const TriMesh *> &meshes,
                std::vector<const DistanceField *> &boundaryFields,
                bool useLeapfrog,
                REAL distanceTolerance = 0.0,
                bool useBoundary = true,
                bool rasterize = false,
                const Vector3Array *listeningPositions = NULL,
                const char *outputFile = NULL,
                WriteCallback *callback = NULL,
                int subSteps = 1,
                int N = 1 );

        // Destructor
        virtual ~WaveSolver();

        // Time steps the system in the given interval
        void solveSystem( REAL startTime, REAL endTime,
                const BoundaryEvaluator &bcEvaluator );

        void initSystem( REAL startTime,
                VECTOR *initialPressure = NULL,
                VECTOR *initialVelocity = NULL );

        // Takes a single time step
        virtual bool stepSystem( const BoundaryEvaluator &bcEvaluator );

        virtual void vertexPressure( const Tuple3i &index, VECTOR &pressure );

        virtual void writeWaveOutput() const;

        virtual const Tuple3i   &fieldDivisions() const
        {
            return _laplacian.fieldDivisions();
        }

        virtual Vector3d         fieldPosition( const Tuple3i &index ) const
        {
            return _laplacian.fieldPosition( index );
        }

        virtual  Vector3d        fieldPosition( int index ) const
        {
            return _laplacian.fieldPosition( index );
        }

#if 0
        inline const TriMesh    &mesh() const
        {
            return _laplacian.mesh();
        }
#endif
        virtual const vector<const TriMesh *> &meshes() const
        {
            return _laplacian.meshes();
        }

        virtual int              numCells() const
        {
            return _laplacian.numCells();
        }

        inline const IntArray   &ghostCells() const
        {
            return _laplacian.ghostCells();
        }

        inline const IntArray   &interfacialCells() const
        {
            return _laplacian.interfacialCells();
        }

        virtual const Vector3Array *listeningPositions() const
        {
            return _listeningPositions;
        }

        virtual REAL             fieldDiameter() const
        {
            return _laplacian.fieldDiameter();
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
            return _laplacian.field().bbox().center();
        }

    protected:

    private:
        void                     stepLeapfrog(
                const BoundaryEvaluator &bcEvaluator );

        // Time stepping based on standard centered difference of 2nd
        // time derivative
        void                     stepCenteredDifference(
                const BoundaryEvaluator &bcEvaluator );

    private:
        static const REAL        WAVE_SPEED = 343.0;

        // Laplacian discretization for this domain
        Laplacian                _laplacian;

        // Pretty self-explanatory
        REAL                     _timeStep;
        REAL                     _currentTime;
        int                      _timeIndex;
        int                      _subSteps;

        bool                     _useLeapfrog;

        bool                     _useMACGrid;

        // Store the last two time steps for a second order time differencing
        // scheme.  We form p2 (the current time step) using p0 and p1
        MATRIX                   _p0;
        MATRIX                   _p1;
        MATRIX                   _p2;

        MATRIX                   _v;
        MATRIX                   _p;
        MATRIX                   _a;

        // Number of fields to solve for
        int                      _N;

        // Workspace for Laplacian application
        MATRIX                   _workspace;

        const Vector3Array      *_listeningPositions;
        const char              *_outputFile;

        std::vector<std::vector<FloatArray> >
            _waveOutput;

        VECTOR                   _listenerOutput;

        WriteCallback           *_callback;

        Timer<false>             _laplacianTimer;
        Timer<false>             _boundaryTimer;
        Timer<false>             _stepTimer;
        Timer<false>             _algebraTimer;
        Timer<false>             _memoryTimer;
        Timer<false>             _writeTimer;

};

#endif
