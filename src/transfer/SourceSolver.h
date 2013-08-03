//////////////////////////////////////////////////////////////////////
// SourceSolver.h: Interface
//
//////////////////////////////////////////////////////////////////////

#ifndef SOURCE_SOLVER_H
#define SOURCE_SOLVER_H

#include <distancefield/distanceField.h>

#include <linearalgebra/MATRIX.h>
#include <linearalgebra/MATRIX3.h>
#include <linearalgebra/VECTOR.h>

#include <SETTINGS.h>
#include <TYPES.h>

#include <util/IO.h>
#include <util/Evaluator.h>
#include <util/MERSENNETWISTER.h>

#include <mesh/Triangle.h>

#include <math/Function.h>

class TriMesh;

//////////////////////////////////////////////////////////////////////
// SourceSolver class
//
// Used to solve rigid acceleration radiation problems for
// a triangle mesh.
//////////////////////////////////////////////////////////////////////
class SourceSolver {
	public:
		SourceSolver( const TriMesh &mesh );

		// Destructor
		virtual ~SourceSolver();

    // Principle acceleration directions
    enum Direction {
      TRANS_X = 0,
      TRANS_Y,
      TRANS_Z,
      ROT_X,
      ROT_Y,
      ROT_Z
    };

    // Assembles a system to solve for source coefficients.
    // Surface points with area weights are used.
    //
    // nTerms: Number of spatial source terms (total number
    //         of terms depends on the associated time scales)
    // timeScale: The interpolation width used for Bspline
    //            pulses emitted from the source
    // timeStep: Time integration step
    // T: Pulse is supported in interval [-T,T]
    // c: Speed of sound
    // acceleration: Time-dependent acceleration function
    // pulseInterval: Can optionally choose to send out
    //                pulses at an interval different from
    //                timeScale (the time scale associated
    //                with the pulses themselves remains the
    //                same)
    void buildCollocationSystem( Direction direction,
                                 int nTerms,
                                 Real timeScale,
                                 Real timeStep,
                                 Real T, Real c, Real density,
                                 RealFunction &acceleration,
                                 MATRIX &system,
                                 VECTOR &rhs,
                                 Real pulseInterval );

    // Builds a set of sources to approximate a rigid acceleration
    // boundary condition to within the given tolerance.
    void trainSourceList( Direction direction,
                          Real timeScale,
                          Real timeStep,
                          Real T, Real c, Real density,
                          RealFunction &acceleration,
                          // Distance field parameters
                          int resolution, const char *filePrefix,
                          int pointsPerVoxel, bool randomScatter,
                          // Training parameters
                          int sourcesPerIteration,
                          int candidatesPerIteration,
                          int maxSourcesPerBatch,
                          Real tolerance );

    // Used to compare sources
    struct SourceRanking {
      SourceRanking()
        : _index( -1 ),
          _quality( 0.0 )
      {
      }

      SourceRanking( int index, Real quality )
        : _index( index ),
          _quality( quality )
      {
      }

      bool operator<( const SourceRanking &s )
      {
        return _quality > s._quality;
      }

      int                    _index;
      Real                   _quality;
    };

	protected:

  private:
    void clear();

    // Stores the spatial point index, multipole term index, and
    // time emission index.
    struct SourceIndex {
      SourceIndex()
        : _sourcePointIdx( -1 ),
          _multipoleIdx( -1 ),
          _timeIdx( -1 )
      {
      }

      SourceIndex( int sourcePointIdx, int multipoleIdx, int timeIdx )
        : _sourcePointIdx( sourcePointIdx ),
          _multipoleIdx( multipoleIdx ),
          _timeIdx( timeIdx )
      {
      }

      void print() const
      {
        cout << "SourceIndex:" << endl;
        cout << SDUMP( _sourcePointIdx ) << endl;
        cout << SDUMP( _multipoleIdx ) << endl;
        cout << SDUMP( _timeIdx ) << endl << endl;
      }

      int                    _sourcePointIdx;
      int                    _multipoleIdx;
      int                    _timeIdx;
    };

    typedef std::vector<SourceIndex>   SourceList;

    // Stores a list of active (unused) pulses for a given source
    class SourcePulseList {
      public:
        SourcePulseList( const VEC3F &sourcePosition,
                         const Vector3Array &meshPoints,
                         // Needed to figure out pulse range
                         Real timeScale, int S, Real T, Real c )
          : _pos( sourcePosition ),
            _pulseActive( 0 )
        {
          SourceSolver::pulseIndexRange( _pos, meshPoints, timeScale, S, T, c,
                                         _pulseMin, _pulseMax );

          _activePulses = _pulseMax - _pulseMin + 1;

          TRACE_ASSERT( _activePulses >= 1,
                        "Invalid number of active pulses" );

          _pulseActive.resize( _activePulses, true );
        }

        virtual ~SourcePulseList()
        {
        }

        inline const VEC3F &pos() const
        {
          return _pos;
        }

        inline int numPulses() const
        {
          return _pulseMax - _pulseMin + 1;
        }

        inline int numActivePulses() const
        {
          return _activePulses;
        }

        inline bool pulseActive( int pulse_idx ) const
        {
          return _pulseActive[ pulse_idx - _pulseMin ];
        }

        inline int pulseMin() const
        {
          return _pulseMin;
        }

        inline int pulseMax() const
        {
          return _pulseMax;
        }

        // By choosing a pulse for this source location, we
        // deactivate it
        inline void selectPulse( int pulse_idx )
        {
          if ( _pulseActive[ pulse_idx - _pulseMin ] )
          {
            _activePulses -= 1;
          }

          _pulseActive[ pulse_idx - _pulseMin ] = false;
        }

        inline void deselectPulse( int pulse_idx )
        {
          if ( !_pulseActive[ pulse_idx - _pulseMin ] )
          {
            _activePulses += 1;
          }

          _pulseActive[ pulse_idx - _pulseMin ] = true;
        }

      private:
        VEC3F                _pos;
        BoolArray            _pulseActive;
        int                  _activePulses;

        int                  _pulseMin, _pulseMax;

    };

    // Unit acceleration vector in the appropriate direction
    static VEC3F getAccelerationVector( Direction direction );

    // Gets the directional acceleration term at a point on the mesh
    Real getAccelerationTerm( Direction direction,
                              const VEC3F &accelerationVector,
                              const VEC3F &position,
                              const VEC3F &normal );

    // Find closest and farthest points relative to the mesh's
    // center of mass
    void findPositionRange();

    // Given pulse time scale timeScale, support width S, acceleration pulse
    // interval [-T,T] and speed of sound c, computes the minimum
    // and maximum indices for pulses which can lead to acceleration
    // boundary conditions.
    void pulseIndexRange( Real timeScale, int S, Real T, Real c,
                          int &pulseMin, int &pulseMax );

    // In this version, we can specify a distinct interval between
    // pulses (other then the time scale)
    void pulseIndexRange( Real timeScale, Real pulseInterval,
                          int S, Real T, Real c,
                          int &pulseMin, int &pulseMax );

    // Builds a column in a collocation system
    //
    // timeScale: The interpolation width used for Bspline pulses
    //            emitted froma source
    // timeStep: Time integation step
    // timePoints: Number of points to use in time integration
    // spatialTerm: Index of the multipole to use
    // pulseIndex: Index for the pulse function
    // pulseMin: Minimum pulse function index
    //
    // Optionally we may build a row, rather than a column.  This
    // may be more convenient for matrix multiplication, given
    // row-major ordering.
    void buildCollocationColumn( const SourceIndex &sourceTerm,
                                 Real timeScale,
                                 Real timeStep,
                                 Real tStart,
                                 int timePoints,
                                 Real c,
                                 Function &pulseFunction,
                                 Vector3Array &meshPoints,
                                 Vector3Array &meshNormals,
                                 FloatArray &meshWeights,
                                 MATRIX &system,
                                 int systemColumn,
                                 bool buildRow = false );

    // Builds a column in the collocation system
    //
    // timeScale: The interpolation width used for Bspline
    //            pulses emitted from the source
    // timeStep: Time integration step
    // tStart: Time integration start time
    // timePoints: Number of points to use in time integration
    // spatialTerm: Index of the multipole to use
    // pulseIndex: Index for the pulse function
    // pulseMin: Minimum pulse function index
    void buildCollocationRow( Direction direction,
                              int nTerms,
                              Real timeScale,
                              Real timeStep,
                              Real tStart,
                              int timePoints,
                              int spatialTerm,
                              int pulseIndex, int pulseMin,
                              Real c,
                              Function &pulseFunction,
                              const vector<Triangle *> &triangles,
                              Vector3Array &normals,
                              MATRIX &system,
                              Real pulseInterval );

    // Builds the right hand side of the collocation system
    //
    // timeStep: Time integration step
    // tStart: Time integration start time
    // timePoints: Number of points to use in time integration
    // acceleration: Acceleration function being considered
    void buildCollocationRHS( Direction direction,
                              Real timeStep,
                              Real tStart,
                              int timePoints,
                              RealFunction &acceleration,
                              const vector<Triangle *> &triangles,
                              Vector3Array &normals,
                              const VEC3F &accelerationVector,
                              Real density,
                              VECTOR &rhs );

    // Same as the above, but just takes a list of points instead of
    // triangles
    void buildCollocationRHS( Direction direction,
                              Real timeStep,
                              Real tStart,
                              int timePoints,
                              RealFunction &acceleration,
                              const Vector3Array &meshPoints,
                              const Vector3Array &meshNormals,
                              const FloatArray &meshWeights,
                              const VEC3F &accelerationVector,
                              Real density,
                              VECTOR &rhs );

    //////////////////////////////////////////////////////////////////////
    // Functions for adaptively generating source locations in space
    // and time to try to fit a boundary condition
    //////////////////////////////////////////////////////////////////////

    void initTrainingWorkspace( int candidatesPerIteration,
                                int nRows );

    // Generates a set of collocation points at which we will
    // attempt to enforce a boundary condition
    void generateCollocationPoints( Vector3Array &points,
                                    Vector3Array &normals,
                                    FloatArray &weights );

    // For each posible source position, build a list of source
    // pulses that have meaningful contributions when produced from
    // that source location.
    //
    // resolution, filePrefix, pointsPerVoxel and randomScatter
    // have the same meaning as the function below.
    void initSourcePulseList( int resolution, const char *filePrefix,
                              int pointsPerVoxel, bool randomScatter,
                              const Vector3Array &meshPoints,
                              Real timeScale, int S, Real T, Real c );

    // Builds a set of candidate source positions inside of an object
    // using a signed distance field.
    //
    // resolution: Resolution of the field in voxels
    // filePrefix: SDFs can take a while to build, so it can be useful
    //             to save them to disk.
    // pointsPerVoxel/randomScatter: Either randomly scatter pointsPerVoxel
    //                               points in each interior voxel, or
    //                               position pointsPerVoxel^3 points in
    //                               a regular grid in the voxel.
    void generateSourcePositions( int resolution, const char *filePrefix,
                                  int pointsPerVoxel, bool randomScatter,
                                  Vector3Array &sourcePoints );

    // Generates points inside a voxel of an SDF.  Only inserts points
    // inside the object
    void generateVoxelPoints( int pointsPerVoxel, bool randomScatter,
                              DistanceField *field,
                              int voxel_idx_x, int voxel_idx_y, int voxel_idx_z,
                              Vector3Array &sourcePoints, Real distance );

    // Scatters points randomly inside a voxel
    void generateRandomVoxelPoints(
                              int pointsPerVoxel,
                              DistanceField *field,
                              int voxel_idx_x, int voxel_idx_y, int voxel_idx_z,
                              Vector3Array &sourcePoints, Real distance );

    // Generates points uniformly inside a voxel
    void generateUniformVoxelPoints(
                              int pointsPerVoxel,
                              DistanceField *field,
                              int voxel_idx_x, int voxel_idx_y, int voxel_idx_z,
                              Vector3Array &sourcePoints, Real distance );

    // Randomly generates a set of candidate source points.  We
    // will choose Npos source positions, then Npulse pulse start
    // times for each position.
    void generateCandidatePoints( const Vector3Array &meshPoints,
                                  int numCandidates,
                                  SourceList &sources );

    // Chooses a unique candidate point not yet included in the system
    // and not yet included in the current set of candidate sources
    IntPair generateUniqueCandidate( const Vector3Array &meshPoints,
                                     MERSENNETWISTER &generator ) const;

    // For a set of sources, generate the coefficients for the source
    // system and place them in the workspace
    void buildSourceCoefficients( const SourceList &sources,
                                  Real timeScale,
                                  Real timeStep,
                                  Real tStart,
                                  int timePoints,
                                  Real c,
                                  Function &pulseFunction,
                                  Vector3Array &meshPoints,
                                  Vector3Array &meshNormals,
                                  FloatArray &meshWeights );

    // Ranks sources in order of quality.
    // Here, we defined "best", as the source whose rows are most
    // orthogonal to the current residual.  That is, columns having
    // the smallest dot product with the residual.
    //
    // This assumes that buildSourceCoefficients has been called already
    void rankCandidateSources( const SourceList &sources,
                               std::vector<SourceRanking> &ranking );

    // Given a ranking of candidate source columns, append some number of
    // columns to the adaptive system.
    //
    // Also clears any changes made to _sourcePulses for columns
    // not selected.
    void appendCandidateColumns( const SourceList &sources,
                                 const std::vector<SourceRanking> &ranking,
                                 int sourcesPerIteration,
                                 SourceList &finalSourceList );

    // Updates normal equations for a least squares solve
    void updateNormalEquations( const SourceList &sources,
                                const std::vector<SourceRanking> &ranking,
                                int sourcesPerIteration );

    // Enlarges systems if necessary
    void checkSystemSize( int sourcesPerIteration );

    // Solves the current least squares system and updates the residual.
    void updateResidual();

    // Assuming that pulses are emitted from a source at position
    // sourcePos at interval timeScale, and with support support * timeScale,
    // and given a point on the mesh, figure out the index range of
    // pulses that must be considered in order to consider all contributions
    // to that point.
    void pulseIndexRange( const VEC3F &sourcePoint, const VEC3F &meshPoint,
                          Real timeScale, int S, Real T, Real c,
                          int &pulseMin, int &pulseMax );

    void writeAdaptiveSystem( const char *fileName )
    {
      MATRIX                 adaptiveSystem( _adaptiveSystemRows,
                                             _sourceWorkspace.cols(),
                                             _adaptiveSystem );

      adaptiveSystem.write( fileName );
    }

    // For testing
    friend int main( int argc, char ** argv );

    //friend bool operator<( const SourceRanking &s1, const SourceRanking &s2 );

  private:
    // Similar to the above, but figures out the index range if we want
    // to look at a single source point interacting with *all* points
    // on the mesh.
    static void pulseIndexRange( const VEC3F &sourcePoint,
                                 const Vector3Array &meshPoints,
                                 Real timeScale, int S, Real T, Real c,
                                 int &pulseMin, int &pulseMax );

    // For a given pulse, figure out the range of times over which it
    // interacts with the mesh
    static void pulseTimeRange( const VEC3F &sourcePoint,
                                const SourceIndex &idx,
                                const Vector3Array &meshPoints,
                                Real timeScale, int S, Real c,
                                Real &tMin, Real &tMax );

    // Same as the above, but only for a specific receiver point
    static void pulseTimeRange( const VEC3F &sourcePoint,
                                const SourceIndex &idx,
                                const VEC3F &meshPoint,
                                Real timeScale, int S, Real c,
                                Real &tMin, Real &tMax );

    // Finds the closest and farthest points in the given list relative
    // to the given source location
    static void distanceBounds( const Vector3Array &meshPoints,
                                const VEC3F &sourcePoint,
                                Real &rMin, Real &rMax );

	private:
    const TriMesh               &_mesh;

    // Minimum and maximum distances to the surface of the
    // mesh from its center of mass
    Real                         _rMin;
    Real                         _rMax;

    // Least squares system and right-hand side
    MATRIX                       _system;
    VECTOR                       _rhs;
    VECTOR                       _rhsFull;

    // For adaptive generation of source points, store a list
    // of currently active space/time pairs
    std::vector<SourcePulseList> _sourcePulses;

    // Workspace for generating source coefficients during
    // an adaptive training process
    MATRIX                       _sourceWorkspace;

    // We also want the norms of each row in the workspace
    VECTOR                       _sourceWorkspaceNorms;

    // Workspace for accumulating adaptive system.  We need
    // to keep a copy around since LAPACK's least squares
    // solver overwrites the system.
    //
    // Note: For convenience (due to row-major ordering) we
    //       will store the transpose of the system, so we
    //       should keep this in mind when using the least
    //       squares solver.
    Real                        *_adaptiveSystem;
    Real                        *_adaptiveSystemCopy;

    // We may only use a subset of the system rows
    int                          _adaptiveSystemRows;
    int                          _adaptiveSystemCapacity;

    // Store the normal equations for the system.
    Real                        *_adaptiveNormalEquations;
    Real                        *_adaptiveNormalEquationsCopy;

    // Leading dimension for the normal equations (so that we don't
    // have to reinitialize every time we want to regenerate them)
    int                          _adaptiveNormalLDA;

    // We will also need to keep a copy of the right hand
    // side, and a residual vector
    VECTOR                       _rhsCopy;
    VECTOR                       _residual;

    const static IndexRange      DIPOLE_RANGE;

    MERSENNETWISTER              _generator;

};

static bool operator<( const SourceSolver::SourceRanking &s1,
                       const SourceSolver::SourceRanking &s2 )
{
  return s1._quality > s2._quality;
}

#endif
