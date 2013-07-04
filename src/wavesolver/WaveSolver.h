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

#include <mesh/TriMesh.h>

#include <util/BoundingBox.h>
#include <util/Evaluator.h>
#include <util/timer.h>

#include "Laplacian.h"

#include <SETTINGS.h>
#include <TYPES.h>

#include <boost/function.hpp>

//////////////////////////////////////////////////////////////////////
// Solver interface
//
//////////////////////////////////////////////////////////////////////
class Solver {
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
    virtual VEC3F fieldPosition( const Tuple3i &index ) const = 0;
    virtual VEC3F fieldPosition( int index ) const = 0;
    virtual const Vector3Array *listeningPositions() const = 0;
    virtual Real fieldDiameter() const = 0;
    virtual void vertexPressure( const Tuple3i &index, VECTOR &pressure ) = 0;
    virtual int numCells() const = 0;
    virtual Real currentSimTime() const = 0;
    virtual VEC3F sceneCenter() const = 0;

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

        // Provide the size of the domain (bbox), finite difference division
        // size, and a signed distance function for the interior boundary.
		WaveSolver( Real timeStep,
                    const BoundingBox &bbox, Real cellSize,
                    const TriMesh &mesh,
                    const DistanceField &distanceField,
                    bool useLeapfrog,
                    Real distanceTolerance = 0.0,
                    bool useBoundary = true,
                    bool rasterize = false,
                    const Vector3Array *listeningPositions = NULL,
                    const char *outputFile = NULL,
                    WriteCallback *callback = NULL,
                    int subSteps = 1,
                    int N = 1 );

        // Provide the size of the domain (bbox), finite difference division
        // size, and a list of SDFs for the interior boundary
        WaveSolver( Real timeStep,
                    const BoundingBox &bbox, Real cellSize,
                    std::vector<const TriMesh *> &meshes,
                    std::vector<const DistanceField *> &boundaryFields,
                    bool useLeapfrog,
                    Real distanceTolerance = 0.0,
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
        void solveSystem( Real startTime, Real endTime,
                          const BoundaryEvaluator &bcEvaluator );

        void initSystem( Real startTime,
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

        virtual VEC3F            fieldPosition( const Tuple3i &index ) const
                                 {
                                   return _laplacian.fieldPosition( index );
                                 }

        virtual  VEC3F           fieldPosition( int index ) const
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
        
        virtual Real             fieldDiameter() const
                                 {
                                   return _laplacian.fieldDiameter();
                                 }

        virtual  int             N() const
                                 {
                                   return _N;
                                 }

        virtual Real             currentSimTime() const
                                 {
                                   return _timeStep * (Real)_timeIndex;
                                 }
        
        virtual VEC3F            sceneCenter() const
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
        static const Real        WAVE_SPEED = 343.0;

        // Laplacian discretization for this domain
        Laplacian                _laplacian;

        // Pretty self-explanatory
        Real                     _timeStep;
        Real                     _currentTime;
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

        Timer                    _laplacianTimer;
        Timer                    _boundaryTimer;
        Timer                    _stepTimer;
        Timer                    _algebraTimer;
        Timer                    _memoryTimer;
        Timer                    _writeTimer;

};

#endif
