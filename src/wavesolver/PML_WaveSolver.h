//////////////////////////////////////////////////////////////////////
// PML_WaveSolver.h: Interface for the PML_WaveSolver class
//
//////////////////////////////////////////////////////////////////////

#ifndef PML_WAVE_SOLVER_H
#define PML_WAVE_SOLVER_H

#include <distancefield/distanceField.h>

#include <field/ScalarField.h>

#include <linearalgebra/MATRIX.h>
#include <linearalgebra/VECTOR.h>

#include <mesh/TriMesh.h>

#include <util/BoundingBox.h>
#include <util/Evaluator.h>
#include <util/timer.h>

#include "MAC_Grid.h"
#include "WaveSolver.h"

#include <SETTINGS.h>
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
    typedef boost::function<void (const vector<vector<FloatArray> >&w)>
                                                                WriteCallback;

    // Provide the size of the domain (bbox), finite difference division
    // size, and a signed distance function for the interior boundary.
		PML_WaveSolver( Real timeStep,
                    const BoundingBox &bbox, Real cellSize,
                    const TriMesh &mesh,
                    const DistanceField &distanceField,
                    Real distanceTolerance = 0.0,
                    bool useBoundary = true,
                    const Vector3Array *listeningPositions = NULL,
                    const char *outputFile = NULL,
                    WriteCallback *callback = NULL,
                    int subSteps = 1,
                    int N = 1,
                    Real endTime = -1.0 );

    // Provide the size of the domain (bbox), finite difference division
    // size, and a list of SDFs for the interior boundary
		PML_WaveSolver( Real timeStep,
                    const BoundingBox &bbox, Real cellSize,
                    std::vector<const TriMesh *> &meshes,
                    std::vector<const DistanceField *> &boundaryFields,
                    Real distanceTolerance = 0.0,
                    bool useBoundary = true,
                    const Vector3Array *listeningPositions = NULL,
                    const char *outputFile = NULL,
                    WriteCallback *callback = NULL,
                    int subSteps = 1,
                    int N = 1,
                    Real endTime = -1.0 );

		// Destructor
		virtual ~PML_WaveSolver();

    void setPMLBoundaryWidth( Real width, Real strength );

    // Time steps the system in the given interval
    void solveSystem( Real startTime, Real endTime,
                      const BoundaryEvaluator &bcEvaluator );

    void initSystem( Real startTime );

    // Takes a single time step
    virtual bool stepSystem( const BoundaryEvaluator &bcEvaluator );

    virtual void vertexPressure( const Tuple3i &index, VECTOR &pressure );

    virtual void writeWaveOutput() const;

    virtual  const Tuple3i  &fieldDivisions() const
                             {
                               return _grid.pressureFieldDivisions();
                             }

    virtual VEC3F            fieldPosition( const Tuple3i &index ) const
                             {
                               return _grid.pressureFieldPosition( index );
                             }

    virtual VEC3F            fieldPosition( int index ) const
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
    
    virtual Real             fieldDiameter() const
                             {
                               return _grid.fieldDiameter();
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
                               return _grid.pressureField().bbox().center();
                             }

    void                     setZSlice( int slice )
                             {
                               _zSlice = slice;
                             }

	protected:

  private:
    void                     stepLeapfrog(
                                        const BoundaryEvaluator &bcEvaluator );

	private:
    static const Real        WAVE_SPEED = 343.0;

    // Discretization for the domain
    MAC_Grid                 _grid;

    // Pretty self-explanatory
    Real                     _timeStep;
    Real                     _currentTime;
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

    Real                     _endTime;

    const Vector3Array      *_listeningPositions;
    const char              *_outputFile;

    std::vector<std::vector<FloatArray> >
                             _waveOutput;
  
    VECTOR                   _listenerOutput;

    WriteCallback           *_callback;

    Timer                    _gradientTimer;
    Timer                    _divergenceTimer;
    Timer                    _stepTimer;
    Timer                    _algebraTimer;
    Timer                    _memoryTimer;
    Timer                    _writeTimer;

    // Optionally write a 2D slice out of the finite difference grid
    int                      _zSlice;
    MATRIX                   _sliceData;

};

#endif
