//////////////////////////////////////////////////////////////////////
// SourceSolver.cpp: Implementation
//
//////////////////////////////////////////////////////////////////////

#include "SourceSolver.h"

#include <math/BSpline.h>
#include <math/MultipoleFunction.h>

#include <mesh/Triangle.h>
#include <mesh/TriMesh.h>

#include <util/IO.h>

#include <algorithm>

#include <float.h>
#include <limits.h>

#ifdef USE_OMP
#include <omp.h>
#endif

const IndexRange SourceSolver::DIPOLE_RANGE( 1, 3 );
//const IndexRange SourceSolver::DIPOLE_RANGE( 4, 9 );

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
SourceSolver::SourceSolver( const TriMesh &mesh )
  : _mesh( mesh ),
    _adaptiveSystem( NULL ),
    _adaptiveSystemCopy( NULL ),
    _adaptiveNormalEquations( NULL ),
    _adaptiveSystemRows( 0 ),
    _generator( false )
{
  findPositionRange();
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
SourceSolver::~SourceSolver()
{
  clear();
}

//////////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////////
void SourceSolver::buildCollocationSystem( Direction direction,
                                           int nTerms,
                                           Real timeScale,
                                           Real timeStep,
                                           Real T, Real c, Real density,
                                           RealFunction &acceleration,
                                           MATRIX &system,
                                           VECTOR &rhs,
                                           Real pulseInterval )
{
  VEC3F                      accelerationVector;
  Vector3Array               normals;
  const vector<Triangle *>  &triangles = _mesh.triangles();
  int                        pulseMin, pulseMax, numPulses;
  Function                   pulseFunction;
  int                        timePoints;
  Real                       timeLength;
  Real                       tStart;

  accelerationVector = getAccelerationVector( direction );

  _mesh.computeVertexNormals( normals );

  // FIXME: For now, just get a 9th order spline
  pulseFunction = BSpline::getSpline( 9, timeScale );

  if ( pulseInterval > 0.0 )
  {
    pulseIndexRange( timeScale, pulseInterval,
                     pulseFunction.supportSamples(), T, c,
                     pulseMin, pulseMax );
  }
  else
  {
    pulseIndexRange( timeScale, pulseFunction.supportSamples(), T, c,
                     pulseMin, pulseMax );
  }

  numPulses = pulseMax - pulseMin + 1;

  // Figure out how many time integration points to use
#if 0
  timeLength = timeScale
    * (Real)( numPulses + 2 * pulseFunction.supportSamples() );
#endif
  timeLength = 2 * T + 2.0 * timeScale * (Real)pulseFunction.supportSamples();
#if 0
  timeLength = 2.0 * T;
#endif

  cout << SDUMP( numPulses ) << endl;
  cout << SDUMP( pulseMin ) << endl;
  cout << SDUMP( pulseMax ) << endl;

#if 0
  tStart = timeScale
    * (Real)( pulseMin - pulseFunction.supportSamples() );
#endif
  tStart = -1.0 * T - timeScale * (Real)pulseFunction.supportSamples();
#if 0
  tStart = -1.0 * T;
#endif

  cout << SDUMP( tStart ) << endl;

  timePoints = (int)ceil( timeLength / timeStep ) + 1;

  cout << "In SourceSolver::buildCollocationSystem" << endl;
  cout << SDUMP( numPulses ) << endl;
  cout << SDUMP( nTerms * numPulses ) << endl;
  cout << endl;
  cout << SDUMP( timePoints ) << endl;
  cout << SDUMP( triangles.size() ) << endl;
  cout << SDUMP( triangles.size() * timePoints ) << endl;

  // Build the initial system
  system.resizeAndWipe( triangles.size() * timePoints, nTerms * numPulses );
  rhs.resizeAndWipe( triangles.size() * timePoints );

  for ( int pulseIdx = pulseMin; pulseIdx <= pulseMax; pulseIdx++ )
  for ( int multipoleIdx = 0; multipoleIdx < nTerms; multipoleIdx++ )
  {
    buildCollocationRow( direction, nTerms, timeScale, timeStep, tStart,
                         timePoints, multipoleIdx, pulseIdx, pulseMin,
                         c, pulseFunction, triangles, normals, system,
                         pulseInterval );
  }

  printf( "\n" );

  buildCollocationRHS( direction, timeStep, tStart, timePoints,
                       acceleration, triangles, normals, accelerationVector,
                       density, rhs );
}

//////////////////////////////////////////////////////////////////////
// Builds a set of sources to approximate a rigid acceleration
// boundary condition to within the given tolerance.
//////////////////////////////////////////////////////////////////////
void SourceSolver::trainSourceList( Direction direction,
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
                                    Real tolerance )
{
  VEC3F                      accelerationVector;
  Vector3Array               meshPoints;
  Vector3Array               meshNormals;
  FloatArray                 meshWeights;
  Function                   pulseFunction;
  int                        timePoints;
  Real                       timeLength;
  Real                       tStart;
  Real                       residualError;
  int                        S;
  SourceList                 candidateSources;
  SourceList                 finalSources;
  vector<SourceRanking>      sourceRanking;
  int                        iteration = 0;

  accelerationVector = getAccelerationVector( direction );

  generateCollocationPoints( meshPoints, meshNormals, meshWeights );

  // FIXME: For now, just get a 9th order spline
  pulseFunction = BSpline::getSpline( 5, timeScale );

  S = pulseFunction.supportSamples();

#if 0
  timeLength = 2 * T + 2.0 * timeScale * (Real)pulseFunction.supportSamples();
  tStart = -1.0 * T - timeScale * (Real)pulseFunction.supportSamples();
#endif
  timeLength = 2.0 * T + 8.0 * timeScale * (Real)pulseFunction.supportSamples();
  tStart = -1.0 * T - 4.0 * timeScale * (Real)pulseFunction.supportSamples();
#if 0
  timeLength = 4.0 * timeScale * (Real)pulseFunction.supportSamples() / c;
  timeLength += 4.0 * T;
  tStart = -2.0 * timeScale * (Real)pulseFunction.supportSamples() / c;
  tStart -= 2.0 * T;
#endif

  timePoints = (int)ceil( timeLength / timeStep ) + 1;

  // Initialize workspaces
  initTrainingWorkspace( candidatesPerIteration,
                         timePoints * meshPoints.size() );

  initSourcePulseList( resolution, filePrefix, pointsPerVoxel, randomScatter,
                       meshPoints, timeScale, S, T, c );

  candidateSources.reserve( candidatesPerIteration );
  sourceRanking.reserve( candidatesPerIteration );

#if 0
  cout << "In SourceSolver::trainSourceList" << endl;
  cout << SDUMP( tStart ) << endl;
  cout << SDUMP( timePoints ) << endl;
  cout << SDUMP( meshPoints.size() ) << endl;
  cout << SDUMP( meshPoints.size() * timePoints ) << endl;
#endif

  buildCollocationRHS( direction, timeStep, tStart, timePoints,
                       acceleration, meshPoints, meshNormals, meshWeights,
                       accelerationVector, density, _rhs );

  _rhsFull.copyInplace( _rhs );

  // Set up the initial residual
  _residual.copyInplace( _rhs );

  residualError = _residual.norm2() / _rhs.norm2();

  _rhs.write( "trainingRHS.vector" );

  while ( residualError > tolerance )
  {
    // Reset the system if it has become too large
    if ( _adaptiveSystemRows + sourcesPerIteration > maxSourcesPerBatch )
    {
      printf( "\n\n*** Resetting system with error %f ***\n\n", residualError );

      // The residual is the new RHS
      _rhs.copyInplace( _residual );

      _adaptiveSystemRows = 0;
    }

    iteration += 1;
    
    // Select a set of candidates
    generateCandidatePoints( meshPoints, candidatesPerIteration,
                             candidateSources );

    cout << "Building source coefficients" << endl;
    // Build the matrix columns for those candidates
    buildSourceCoefficients( candidateSources,
                             timeScale, timeStep, tStart, timePoints, c,
                             pulseFunction,
                             meshPoints, meshNormals, meshWeights );

#if 0
    char buf[ 1024 ];
    sprintf( buf, "candidates_%02d.matrix", iteration );
    _sourceWorkspace.write( buf );
    break;
#endif

    cout << "Ranking" << endl;
    // Rank the sources from best to worst
    rankCandidateSources( candidateSources, sourceRanking );

    checkSystemSize( sourcesPerIteration );

    updateNormalEquations( candidateSources, sourceRanking,
                           sourcesPerIteration );

    cout << "Appending sources" << endl;
    // Choose the best ones and append them to the system
    appendCandidateColumns( candidateSources, sourceRanking,
                            sourcesPerIteration, finalSources );

    // Do a least squares solve with the new system and update
    // the residual
    cout << "Least squares and residual update" << endl;
    updateResidual();

    residualError = _residual.norm2() / _rhsFull.norm2();

    printf( "Training iteration %d:\n", iteration );
    printf( "\t\tnumSources: %d\n", (int)finalSources.size() );
    printf( "\t\terror: %f\n\n", residualError );

    if ( iteration % 10 == 0 )
    {
      char buf[ 1024 ];

      sprintf( buf, "residual_iteration_%04d.matrix", iteration );

      _residual.write( buf );
    }

#if 0
    if ( iteration % 10 == 0 )
    {
      char buf[ 1024 ];
      
      sprintf( buf, "trainingSystem_iteration_%03d.matrix", iteration );

      MATRIX tmp( _adaptiveSystemRows, _sourceWorkspace.cols(),
                  _adaptiveSystem );

      tmp.write( buf );

      sprintf( buf, "normalSystem_iteration_%03d.matrix", iteration );

      MATRIX normalEqn( _adaptiveSystemRows, _adaptiveSystemRows,
                        _adaptiveNormalEquations, _adaptiveNormalLDA );

      normalEqn.write( buf );
    }
#endif
  }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void SourceSolver::clear()
{
  delete[] _adaptiveSystem;
  delete[] _adaptiveSystemCopy;
  delete[] _adaptiveNormalEquations;

  _adaptiveSystem = NULL;
  _adaptiveSystemCopy = NULL;
  _adaptiveNormalEquations = NULL;
}

//////////////////////////////////////////////////////////////////////
// Builds a column in a collocation system
//
// timeScale: The interpolation width used for Bspline pulses
//            emitted froma source
// timeStep: Time integation step
// timePoints: Number of points to use in time integration
// spatialTerm: Index of the multipole to use
// pulseIndex: Index for the pulse function
// pulseMin: Minimum pulse function index
//////////////////////////////////////////////////////////////////////
void SourceSolver::buildCollocationColumn( const SourceIndex &sourceTerm,
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
                                           bool buildRow )
{
  Real                       r;
  VEC3F                      gradient;
  Real                       emissionTime;
  Real                       tMin, tMax;
  Real                       weight;
  Real                       time;
  int                        t_idx_start, t_idx_end;
  int                        S = pulseFunction.supportSamples();
  int                        spatialTerm = sourceTerm._multipoleIdx;
  int                        numPoints = meshPoints.size();

  const VEC3F &sourcePoint
                    = _sourcePulses[ sourceTerm._sourcePointIdx ].pos();

  emissionTime = timeScale * (Real)sourceTerm._timeIdx;
  //emissionTime = timeScale * (Real)sourceTerm._timeIdx / c;

  for ( int point_idx = 0; point_idx < meshPoints.size(); point_idx++ )
  {
    const VEC3F             &meshPoint = meshPoints[ point_idx ];
    const VEC3F             &normal = meshNormals[ point_idx ];

    weight = meshWeights[ point_idx ];

    // Figure out the time range we need to worry about
    pulseTimeRange( sourcePoint, sourceTerm, meshPoint,
                    timeScale, S, c, tMin, tMax );

#if 0
    cout << "Pulse time range: " << endl;
    cout << SDUMP( tMin ) << endl;
    cout << SDUMP( tMax ) << endl << endl;
#endif

    t_idx_start = (int)( ( tMin - tStart ) / timeStep );
    t_idx_end = (int)ceil( ( tMax - tStart ) / timeStep );

    TRACE_ASSERT( t_idx_start >= 0 && t_idx_end <= timePoints,
                  "Pulse time out of range" );

    for ( int t_idx = t_idx_start; t_idx < t_idx_end; t_idx++ )
    {
      time = tStart + timeStep * (Real)t_idx;

      gradient = MultipoleFunction::evaluateMultipoleGradient( 
                                spatialTerm, meshPoint, sourcePoint,
                                pulseFunction.function(),
                                pulseFunction.derivatives(),
                                c,
                                time - emissionTime );

#if 0
      if ( norm2( gradient ) > 0.0 )
      {
        cout << "Found a non-zero gradient" << endl;
        cout << SDUMP( gradient ) << endl;
        cout << SDUMP( weight ) << endl;
        cout << SDUMP( normal ) << endl;
        cout << SDUMP( weight * ( gradient * normal ) ) << endl;
      }
#endif

      if ( buildRow )
      {
        system( systemColumn, t_idx * (int)numPoints + point_idx )
          = weight * ( gradient * normal );
      }
      else
      {
        system( t_idx * (int)numPoints + point_idx, systemColumn )
          = weight * ( gradient * normal );
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Builds a column in the collocation system
// timeScale: The interpolation width used for Bspline
//            pulses emitted from the source
// timeStep: Time integration step
// tStart: Time integration start time
// timePoints: Number of points to use in time integration
// spatialTerm: Index of the multipole to use
// pulseIndex: Index for the pulse function
// pulseMin: Minimum pulse function index
//////////////////////////////////////////////////////////////////////
void SourceSolver::buildCollocationRow( Direction direction,
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
                                        Real pulseInterval )
{
  Real                       pulseStart, pulseEnd, pulseCenter;
  int                        tStartIdx, tEndIdx;
  int                        matrixColumn;
  const VEC3F               &x_cm = _mesh.getCenterOfMass();
  VEC3F                      triangleCentroid;
  Real                       triangleArea;
  VEC3F                      gradient;

  matrixColumn = ( pulseIndex - pulseMin ) * nTerms + spatialTerm;

  // Determine starting and ending times for the actual pulse
  // we are considering.  This depends on whether or not the pulse
  // interval is defined by the time scale, or some other
  // quantity
  if ( pulseInterval > 0.0 )
  {
    pulseStart = pulseInterval * (Real)pulseIndex
        - timeScale * (Real)pulseFunction.supportSamples();

    pulseEnd = pulseStart
        + 2.0 * timeScale * (Real)pulseFunction.supportSamples();

    pulseCenter = pulseInterval * (Real)pulseIndex;
  }
  else
  {
    pulseStart = timeScale
      * (Real)( pulseIndex - pulseFunction.supportSamples() );

    pulseEnd = pulseStart
        + 2.0 * timeScale * (Real)pulseFunction.supportSamples();

    pulseCenter = timeScale * (Real)pulseIndex;
  }

  printf( "Column: M %02d, pIdx %04d, "
          "pMin %04d, "
          " tStart = [ %f ], "
          " pStart = [ %f ], "
          " pEnd = [ %f ]\r",
          //" pCenter = [ %f ]\r",
          spatialTerm, pulseIndex, pulseMin, tStart,
          pulseStart, pulseEnd );

#if 0
  cout << SDUMP( tStart ) << endl;
  cout << SDUMP( pulseStart ) << endl;
  cout << SDUMP( pulseEnd ) << endl;
  cout << SDUMP( pulseCenter ) << endl;
#endif

#if 0
  TRACE_ASSERT( pulseStart >= tStart,
                "Pulse starts outside of integration range" );
  TRACE_ASSERT( pulseEnd <= tStart + timeStep * (Real)timePoints,
                "Pulse ends outside of integration range" );
#endif

  // Determine which time samples we need to worry about
  tStartIdx = (int)( ( pulseStart - tStart ) / timeStep );
  tEndIdx = (int)ceil( ( pulseEnd - tStart ) / timeStep );

  //printf( "Considering time indices [%d, %d]\n", tStartIdx, tEndIdx );

#if 0
  TRACE_ASSERT( tStartIdx >= 0 && tStartIdx < timePoints,
                "Start time index out of range" );
  TRACE_ASSERT( tEndIdx >= 0 && tEndIdx < timePoints,
                "End time index out of range" );

  printf( "Start integration time: %f\n", tStart + timeStep * (Real)tStartIdx );
  printf( "End integration time: %f\n", tStart + timeStep * (Real)tEndIdx );
#endif

  //for ( int tIdx = tStartIdx; tIdx <= tEndIdx; tIdx++ )
  for ( int tIdx = 0; tIdx < timePoints; tIdx++ )
  for ( int spaceIdx = 0; spaceIdx < triangles.size(); spaceIdx++ )
  {
    triangleCentroid = triangles[ spaceIdx ]->getCentroid();
    triangleArea = triangles[ spaceIdx ]->getArea();

    const VEC3F &triangleNormal = triangles[ spaceIdx ]->getNormal();

    gradient = MultipoleFunction::evaluateMultipoleGradient( 
                              spatialTerm, triangleCentroid, x_cm,
                              pulseFunction.function(),
                              pulseFunction.derivatives(),
                              c,
                              tStart + timeStep * (Real)tIdx - pulseCenter );

#if 0
    if ( norm( gradient ) > 0.0 )
    {
      cout << "Non-zero norm!" << endl;
    }
#endif

    system( tIdx * (int)triangles.size() + spaceIdx, matrixColumn )
      = triangleArea * ( gradient * triangleNormal );
  }
}

//////////////////////////////////////////////////////////////////////
// Builds the right hand side of the collocation system
//
// timeStep: Time integration step
// tStart: Time integration start time
// timePoints: Number of points to use in time integration
// acceleration: Acceleration function being considered
//////////////////////////////////////////////////////////////////////
void SourceSolver::buildCollocationRHS( Direction direction,
                                        Real timeStep,
                                        Real tStart,
                                        int timePoints,
                                        RealFunction &acceleration,
                                        const vector<Triangle *> &triangles,
                                        Vector3Array &normals,
                                        const VEC3F &accelerationVector,
                                        Real density,
                                        VECTOR &rhs )
{
  VEC3F                      triangleCentroid;
  Real                       triangleArea;
  Real                       accelerationValue;

  for ( int tIdx = 0; tIdx < timePoints; tIdx++ )
  for ( int spaceIdx = 0; spaceIdx < triangles.size(); spaceIdx++ )
  {
    triangleCentroid = triangles[ spaceIdx ]->getCentroid();
    triangleArea = triangles[ spaceIdx ]->getArea();

    const VEC3F &triangleNormal = triangles[ spaceIdx ]->getNormal();

    accelerationValue = acceleration( tStart + timeStep * (Real)tIdx );

    rhs( tIdx * (int)triangles.size() + spaceIdx )
      = accelerationValue
      * getAccelerationTerm( direction, accelerationVector,
                             triangleCentroid, triangleNormal );
  }
}

//////////////////////////////////////////////////////////////////////
// Same as the above, but just takes a list of points instead of
// triangles
//////////////////////////////////////////////////////////////////////
void SourceSolver::buildCollocationRHS( Direction direction,
                                        Real timeStep,
                                        Real tStart,
                                        int timePoints,
                                        RealFunction &acceleration,
                                        const Vector3Array &meshPoints,
                                        const Vector3Array &meshNormals,
                                        const FloatArray &meshWeights,
                                        const VEC3F &accelerationVector,
                                        Real density,
                                        VECTOR &rhs )
{
  Real                       weight;
  Real                       accelerationValue;

  for ( int tIdx = 0; tIdx < timePoints; tIdx++ )
  for ( int spaceIdx = 0; spaceIdx < meshPoints.size(); spaceIdx++ )
  {
    const VEC3F             &pos = meshPoints[ spaceIdx ];

    weight = meshWeights[ spaceIdx ];

    const VEC3F &normal = meshNormals[ spaceIdx ];

    accelerationValue = acceleration( tStart + timeStep * (Real)tIdx );

    rhs( tIdx * (int)meshPoints.size() + spaceIdx )
      = accelerationValue * weight
      * getAccelerationTerm( direction, accelerationVector, pos, normal );
  }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void SourceSolver::initTrainingWorkspace( int candidatesPerIteration,
                                          int nRows )
{
  int                        work_size;
  int                        normalSize;

  work_size = candidatesPerIteration;
  work_size *= ( DIPOLE_RANGE.second - DIPOLE_RANGE.first  + 1 );

  normalSize = candidatesPerIteration;
  normalSize *= normalSize;

  // Initialize variables
  clear();

  _rhs.resizeAndWipe( nRows );
  _rhsCopy.resizeAndWipe( nRows );
  _rhsFull.resizeAndWipe( nRows );
  _residual.resizeAndWipe( nRows );

  _sourceWorkspace.resizeAndWipe( candidatesPerIteration, nRows );

  _adaptiveSystem = new Real[ candidatesPerIteration * nRows ];
  _adaptiveSystemCopy = new Real[ candidatesPerIteration * nRows ];
  _adaptiveNormalEquations = new Real[ normalSize ];
  _adaptiveNormalEquationsCopy = new Real[ normalSize ];

  _adaptiveSystemRows = 0;
  _adaptiveSystemCapacity = candidatesPerIteration;

  _adaptiveNormalLDA = _adaptiveSystemCapacity;
}

//////////////////////////////////////////////////////////////////////
// Generates a set of collocation points at which we will
// attempt to enforce a boundary condition
//////////////////////////////////////////////////////////////////////
void SourceSolver::generateCollocationPoints( Vector3Array &points,
                                              Vector3Array &normals,
                                              FloatArray &weights )
{
  VEC3F                      triangleCentroid;
  Real                       triangleArea;

  for ( int point_idx = 0; point_idx < _mesh.triangles().size(); point_idx++ )
  {
    triangleCentroid = _mesh.triangles()[ point_idx ]->getCentroid();
    triangleArea = _mesh.triangles()[ point_idx ]->getArea();

    const VEC3F &triangleNormal = _mesh.triangles()[ point_idx ]->getNormal();

    points.push_back( triangleCentroid );
    normals.push_back( triangleNormal );
    weights.push_back( triangleArea );
  }
}

//////////////////////////////////////////////////////////////////////
// For each posible source position, build a list of source
// pulses that have meaningful contributions when produced from
// that source location.
//
// resolution, filePrefix, pointsPerVoxel and randomScatter
// have the same meaning as the function below.
//////////////////////////////////////////////////////////////////////
void SourceSolver::initSourcePulseList( int resolution, const char *filePrefix,
                                        int pointsPerVoxel, bool randomScatter,
                                        const Vector3Array &meshPoints,
                                        Real timeScale, int S, Real T, Real c )
{
  Vector3Array               sourcePoints;
  int                        totalPulses = 0;
  int                        maxPulses = 0;

  generateSourcePositions( resolution, filePrefix, pointsPerVoxel,
                           randomScatter, sourcePoints );

  _sourcePulses.clear();

  for ( int source_idx = 0; source_idx < sourcePoints.size(); source_idx++ )
  {
    _sourcePulses.push_back( SourcePulseList( sourcePoints[ source_idx ],
                                              meshPoints,
                                              timeScale, S, T, c ) );

    totalPulses += _sourcePulses.back().numPulses();
    maxPulses = max( _sourcePulses.back().numPulses(), maxPulses );
  }

  printf( "Candidate pool has %d pulses\n", totalPulses );
  printf( "Maximum pulse length in pool is %d\n", maxPulses );
}

//////////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////////
void SourceSolver::generateSourcePositions( int resolution,
                                            const char *filePrefix,
                                            int pointsPerVoxel,
                                            bool randomScatter,
                                            Vector3Array &sourcePoints )
{
  DistanceField             *field;
  int                        voxelsFound = 0;
  //Real                       distance = -0.0125;
  Real                       distance = -0.00125;

  field = _mesh.buildDistanceField( resolution, filePrefix );

  TRACE_ASSERT( field, "Unable to build distance field" );

  // Check each voxel in the grid to see if part of it lies inside the
  // object
  for ( int voxel_idx_x = 0; voxel_idx_x < resolution - 1; voxel_idx_x++ )
  for ( int voxel_idx_y = 0; voxel_idx_y < resolution - 1; voxel_idx_y++ )
  for ( int voxel_idx_z = 0; voxel_idx_z < resolution - 1; voxel_idx_z++ )
  {
    if ( field->voxelInside( voxel_idx_x, voxel_idx_y, voxel_idx_z,
                             false /* Only partial inclusion needed */,
                             // FIXME: Try only points sufficiently inside
                             //        the object
                             distance ) )
    {
      generateVoxelPoints( pointsPerVoxel, randomScatter, field,
                           voxel_idx_x, voxel_idx_y, voxel_idx_z,
                           sourcePoints, distance );

      voxelsFound++;
    }
  }

  printf( "SourceSolver::generateSourcePositions: Found %d voxels inside\n",
          voxelsFound );
  printf( "SourceSolver::generateSourcePositions: Generated %d sources\n",
          (int)sourcePoints.size() );

  delete field;
}

//////////////////////////////////////////////////////////////////////
// Generates points inside a voxel of an SDF.  Only inserts points
// inside the object
//////////////////////////////////////////////////////////////////////
void SourceSolver::generateVoxelPoints(
                              int pointsPerVoxel, bool randomScatter,
                              DistanceField *field,
                              int voxel_idx_x, int voxel_idx_y, int voxel_idx_z,
                              Vector3Array &sourcePoints, Real distance )
{
  if ( randomScatter )
  {
    generateRandomVoxelPoints( pointsPerVoxel, field,
                               voxel_idx_x, voxel_idx_y, voxel_idx_z,
                               sourcePoints, distance );
  }
  else
  {
    generateUniformVoxelPoints( pointsPerVoxel, field,
                                voxel_idx_x, voxel_idx_y, voxel_idx_z,
                                sourcePoints, distance );
  }
}

//////////////////////////////////////////////////////////////////////
// Scatters points randomly inside a voxel
//////////////////////////////////////////////////////////////////////
void SourceSolver::generateRandomVoxelPoints(
                          int pointsPerVoxel,
                          DistanceField *field,
                          int voxel_idx_x, int voxel_idx_y, int voxel_idx_z,
                          Vector3Array &sourcePoints, Real distance )
{
  MERSENNETWISTER            generator;
  VEC3F                      lowerBound, diagonal;
  VEC3F                      pos;

  lowerBound = field->gridPosition( voxel_idx_x,
                                    voxel_idx_y,
                                    voxel_idx_z );

  diagonal = field->gridPosition( voxel_idx_x + 1,
                                  voxel_idx_y + 1,
                                  voxel_idx_z + 1 );
  diagonal -= lowerBound;

  for ( int point_idx = 0; point_idx < pointsPerVoxel; point_idx++ )
  {
    // Build a random point position
    for ( int i = 0; i < 3; i++ )
    {
      pos[ i ] = lowerBound[ i ] + generator.rand( diagonal[ i ] );
    }

    if ( field->distance( pos ) < distance )
    {
      sourcePoints.push_back( pos );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Generates points uniformly inside a voxel
//////////////////////////////////////////////////////////////////////
void SourceSolver::generateUniformVoxelPoints(
                          int pointsPerVoxel,
                          DistanceField *field,
                          int voxel_idx_x, int voxel_idx_y, int voxel_idx_z,
                          Vector3Array &sourcePoints, Real distance )
{
  VEC3F                      lowerBound, pointDivision;
  VEC3F                      pos;

  lowerBound = field->gridPosition( voxel_idx_x,
                                    voxel_idx_y,
                                    voxel_idx_z );

  pointDivision = field->gridPosition( voxel_idx_x + 1,
                                       voxel_idx_y + 1,
                                       voxel_idx_z + 1 );

  pointDivision -= lowerBound;

  pointDivision /= (Real)( pointsPerVoxel + 1 );

  for ( int point_idx_x = 0; point_idx_x < pointsPerVoxel; point_idx_x++ )
  for ( int point_idx_y = 0; point_idx_y < pointsPerVoxel; point_idx_y++ )
  for ( int point_idx_z = 0; point_idx_z < pointsPerVoxel; point_idx_z++ )
  {
    pos = lowerBound;

    for ( int i = 0; i < 3; i++ )
    {
      pos[ i ] += pointDivision[ i ] * (Real)( i + 1 );
    }

    if ( field->distance( pos ) < distance )
    {
      sourcePoints.push_back( pos );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Randomly generates a set of candidate source points.  We
// will choose Npos source positions, then Npulse pulse start
// times for each position.  For each position we will generate
// three dipole sources, for a total of 3 * Npos * Npulse sources.
//////////////////////////////////////////////////////////////////////
void SourceSolver::generateCandidatePoints( const Vector3Array &meshPoints,
                                            int numCandidates,
                                            SourceList &sources )
{
  // Use this set to guarantee that no duplicates are chosen
  IntPair                    newSource;

  sources.clear();

  for ( int source_idx = 0; source_idx < numCandidates; source_idx++ )
  {
    newSource = generateUniqueCandidate( meshPoints, _generator );

    // Temporarily make this source inactive so that we don't
    // pick it twice
    _sourcePulses[ newSource.first ].selectPulse( newSource.second );

    // Add all dipole sources for this space/time pair
    for ( int dipole_idx = DIPOLE_RANGE.first;
          dipole_idx <= DIPOLE_RANGE.second; dipole_idx++ )
    {
      sources.push_back( SourceIndex( newSource.first,
                                      dipole_idx,
                                      newSource.second ) );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Chooses a unique candidate point not yet included in the system
// and not yet included in the current set of candidate sources
//////////////////////////////////////////////////////////////////////
IntPair SourceSolver::generateUniqueCandidate(
                                    const Vector3Array &meshPoints,
                                    MERSENNETWISTER &generator ) const
{
  bool                       pointFound;
  int                        point_idx, time_idx;
  int                        positionRange = (int)_sourcePulses.size() - 1;

  pointFound = false;

  // Start by finding a source location in space with available pulses
  while ( !pointFound )
  {
    point_idx = generator.randInt( positionRange );

    if ( _sourcePulses[ point_idx ].numActivePulses() > 0 )
    {
      pointFound = true;
    }
  }

  const SourcePulseList     &sourcePoint = _sourcePulses[ point_idx ];

  pointFound = false;

  // Find a suitable time index
  while ( !pointFound )
  {
    time_idx = generator.randInt( sourcePoint.numPulses() - 1 );
    time_idx += sourcePoint.pulseMin();

    if ( sourcePoint.pulseActive( time_idx ) )
    {
      pointFound = true;
    }
  }

  return IntPair( point_idx, time_idx );
}

//////////////////////////////////////////////////////////////////////
// For a set of sources, generate the coefficients for the source
// system and place them in the workspace
//////////////////////////////////////////////////////////////////////
void SourceSolver::buildSourceCoefficients( const SourceList &sources,
                                            Real timeScale,
                                            Real timeStep,
                                            Real tStart,
                                            int timePoints,
                                            Real c,
                                            Function &pulseFunction,
                                            Vector3Array &meshPoints,
                                            Vector3Array &meshNormals,
                                            FloatArray &meshWeights )
{
  if ( _sourceWorkspace.rows() < sources.size()
    || _sourceWorkspace.cols() != timePoints * meshPoints.size() )
  {
    _sourceWorkspace.resizeAndWipe( sources.size(),
                                    timePoints * meshPoints.size() );
  }
  else
  {
    _sourceWorkspace.clear();
  }

  if ( _sourceWorkspaceNorms.size() < sources.size() )
  {
    _sourceWorkspaceNorms.resizeAndWipe( sources.size() );
  }

#pragma omp parallel for schedule(static) default(shared)
  for ( int source_idx = 0; source_idx < sources.size(); source_idx++ )
  {
#if 0
    cout << "Evaluating column " << source_idx << endl;
#endif

    buildCollocationColumn( sources[ source_idx ],
                            timeScale, timeStep, tStart, timePoints, c,
                            pulseFunction,
                            meshPoints, meshNormals, meshWeights,
                            _sourceWorkspace, source_idx,
                            true /* Fill in a row, rather than a column */ );

    _sourceWorkspaceNorms[ source_idx ] = VECTOR::norm2(
                                            _sourceWorkspace.row( source_idx ),
                                            _sourceWorkspace.cols() );
  }
}

//////////////////////////////////////////////////////////////////////
// Ranks sources in order of quality.
// Here, we defined "best", as the source whose rows are most
// orthogonal to the current residual.
//
// This assumes that buildSourceCoefficients has been called already
//////////////////////////////////////////////////////////////////////
void SourceSolver::rankCandidateSources( const SourceList &sources,
                                         vector<SourceRanking> &ranking )
{
  if ( _rhsCopy.size() != _rhs.size() )
  {
    _rhsCopy.resizeAndWipe( _rhs.size() );
  }

  ranking.clear();
  ranking.resize( sources.size() );

  // Multiply the transpose of the candidate source system (stored in
  // _sourceWorkspace) by the current residual.  This will give us
  // the dot product of each source column with the residual.
  MATRIX::gemv( _sourceWorkspace.data(), _residual.data(), _rhsCopy.data(),
                sources.size(), _sourceWorkspace.cols(),
                false /* _sourceWorkspace already stores the transpose */ );

  // Normalize each entry by the magnitude of the corresponding
  // column, and store the results in the ranking list
  for ( int source_idx = 0; source_idx < sources.size(); source_idx++ )
  {
    ranking[ source_idx ]._index = source_idx;

    ranking[ source_idx ]._quality = abs( _rhsCopy( source_idx ) );

#if 0
    printf( "Adding ranking for index %d with quality %f, norm %f.\n",
            ranking[ source_idx ]._index, ranking[ source_idx ]._quality,
            _sourceWorkspaceNorms[ source_idx ] );
#endif

    ranking[ source_idx ]._quality /= _sourceWorkspaceNorms[ source_idx ];
  }

  std::sort( ranking.begin(), ranking.end() );
}

//////////////////////////////////////////////////////////////////////
// Given a ranking of candidate source columns, append some number of
// columns to the adaptive system
//////////////////////////////////////////////////////////////////////
void SourceSolver::appendCandidateColumns( const SourceList &sources,
                                           const vector<SourceRanking> &ranking,
                                           int sourcesPerIteration,
                                           SourceList &finalSourceList )
{
  int                        source_idx;

  for ( int rank_idx = 0; rank_idx < sourcesPerIteration; rank_idx++ )
  {
    source_idx = ranking[ rank_idx ]._index;

    finalSourceList.push_back( sources[ source_idx ] );

    MATRIX::copyRow( _sourceWorkspace.data(), _adaptiveSystem,
                     source_idx, _adaptiveSystemRows,
                     _sourceWorkspace.cols(), _sourceWorkspace.cols() );

#if 0
    cout << "Choosing source " << endl;
    sources[ source_idx ].print();
#endif

    _adaptiveSystemRows += 1;
  }

  for ( int rank_idx = sourcesPerIteration;
        rank_idx < ranking.size(); rank_idx++ )
  {
    source_idx = ranking[ rank_idx ]._index;

    // Make this pulse available again for future iterations
    _sourcePulses[ sources[ source_idx ]._sourcePointIdx ].deselectPulse(
                                              sources[ source_idx ]._timeIdx );
  }
}

//////////////////////////////////////////////////////////////////////
// Updates normal equations for a least squares solve
//////////////////////////////////////////////////////////////////////
void SourceSolver::updateNormalEquations( const SourceList &sources,
                                          const vector<SourceRanking> &ranking,
                                          int sourcesPerIteration )
{
  int                        source_idx;
  Real                      *baseOutputData;
  Real                      *baseInputData;

  // Build the matrix of new rows to add
  for ( int rank_idx = 0; rank_idx < sourcesPerIteration; rank_idx++ )
  {
    source_idx = ranking[ rank_idx ]._index;

    MATRIX::copyRow( _sourceWorkspace.data(), _adaptiveSystemCopy,
                     source_idx, rank_idx,
                     _sourceWorkspace.cols(), _sourceWorkspace.cols() );
  }

  MATRIX tmp1( sourcesPerIteration, _sourceWorkspace.cols(),
               _adaptiveSystemCopy );

#if 0
  tmp1.write( "test1.matrix" );
#endif

  baseOutputData = _adaptiveNormalEquations;
  baseOutputData += _adaptiveSystemRows * _adaptiveNormalLDA;
  baseOutputData += _adaptiveSystemRows;

  // Form the lower triangular part (this should be all we need
  // for a symmetric solver
  MATRIX::syrk( _adaptiveSystemCopy, baseOutputData,
                sourcesPerIteration, _sourceWorkspace.cols(),
                _sourceWorkspace.cols() /* Leading dimension of new data */,
                _adaptiveNormalLDA /* Leading dimension of nml equations */ );

  MATRIX tmp2( sourcesPerIteration, sourcesPerIteration, baseOutputData,
               _adaptiveNormalLDA );

#if 0
  tmp2.write( "testNormal.matrix" );
  abort();
#endif

  if ( _adaptiveSystemRows > 0 )
  {
    baseOutputData = _adaptiveNormalEquations;
    baseOutputData += _adaptiveSystemRows * _adaptiveNormalLDA;

    MATRIX::gemm( _adaptiveSystemCopy, _adaptiveSystem, baseOutputData,
                  sourcesPerIteration, _sourceWorkspace.cols(),
                  _adaptiveSystemRows, _sourceWorkspace.cols(),
                  false /* Do not transpose */, true /* Transpose */,
                  _sourceWorkspace.cols(), _sourceWorkspace.cols(), /* LDA */
                  _adaptiveNormalLDA );
  }
}

//////////////////////////////////////////////////////////////////////
// Enlarges systems if necessary
//////////////////////////////////////////////////////////////////////
void SourceSolver::checkSystemSize( int sourcesPerIteration )
{
  int                        oldSize;

  // If the system is too big, resize it
  if ( _adaptiveSystemRows + sourcesPerIteration > _adaptiveSystemCapacity )
  {
    oldSize = _adaptiveSystemCapacity;

    _adaptiveSystemCapacity *= 2;

    Real                    *newSystem;

    newSystem = new Real[ _adaptiveSystemCapacity * _sourceWorkspace.cols() ];

    // Copy old contents to new system
    MATRIX::copy( newSystem, _adaptiveSystem,
                  _adaptiveSystemRows, _sourceWorkspace.cols() );

    delete[] _adaptiveSystem;

    _adaptiveSystem = newSystem;

    newSystem = new Real[ _adaptiveSystemCapacity * _sourceWorkspace.cols() ];

    delete[] _adaptiveSystemCopy;

    _adaptiveSystemCopy = newSystem;

    newSystem = new Real[ _adaptiveSystemCapacity * _adaptiveSystemCapacity ];

    // Copy all entries from the normal equations.  Need to be
    // careful here about the leading dimension
    for ( int row_idx = 0; row_idx < _adaptiveSystemRows; row_idx++ )
    {
      MATRIX::copyRow( _adaptiveNormalEquations, newSystem,
                       row_idx, row_idx,
                       oldSize, _adaptiveSystemCapacity,
                       _adaptiveSystemRows );
    }

    delete[] _adaptiveNormalEquations;

    _adaptiveNormalEquations = newSystem;

    newSystem = new Real[ _adaptiveSystemCapacity * _adaptiveSystemCapacity ];

    delete[] _adaptiveNormalEquationsCopy;

    _adaptiveNormalEquationsCopy = newSystem;

    _adaptiveNormalLDA = _adaptiveSystemCapacity;
  }
}

//////////////////////////////////////////////////////////////////////
// Solves the current least squares system and updates the residual.
//////////////////////////////////////////////////////////////////////
void SourceSolver::updateResidual()
{
  int                        status;

  if ( _rhsCopy.size() != _rhs.size() )
  {
    _rhsCopy.resizeAndWipe( _rhs.size() );
  }

  if ( _residual.size() != _rhs.size() )
  {
    _residual.resizeAndWipe( _rhs.size() );
  }

  _rhsCopy.copyInplace( _rhs );

#if 0
  // Copy the system so that we don't overwrite its contents
  MATRIX::copy( _adaptiveSystemCopy, _adaptiveSystem,
                _adaptiveSystemRows, _sourceWorkspace.cols() );

  // Solve the least squares problem
  status = MATRIX::leastSquaresFullRank( _adaptiveSystemCopy, _rhsCopy.data(),
                                         _adaptiveSystemRows,
                                         _sourceWorkspace.cols(),
                                         1 /* One right hand side */,
                                         true /* Use transpose */ );
#endif
#if 0
  MATRIX::transposeBLAS( _adaptiveSystemCopy, _adaptiveSystem,
                         _adaptiveSystemRows, _sourceWorkspace.cols() );

  // Solve the least squares problem
  Real *singularValues = new Real[ _sourceWorkspace.cols() ];
  int rank;
  MATRIX::leastSquaresSVD( _adaptiveSystemCopy, _rhsCopy.data(), singularValues,
                           _sourceWorkspace.cols(), _adaptiveSystemRows,
                           1e-12, rank );

  delete[] singularValues;
#endif
  // Try normal equations
  MATRIX::copy( _adaptiveNormalEquationsCopy, _adaptiveNormalEquations,
                _adaptiveNormalLDA, _adaptiveNormalLDA );

  MATRIX::gemv( _adaptiveSystem, _rhs.data(), _rhsCopy.data(),
                _adaptiveSystemRows, _sourceWorkspace.cols(), false );

  MATRIX::solveSPD( _adaptiveNormalEquationsCopy, _rhsCopy.data(),
                    _adaptiveSystemRows, _adaptiveNormalLDA );

  // Form the residual
  _residual.copyInplace( _rhs );

  MATRIX::gemv( _adaptiveSystem, _rhsCopy.data(), _residual.data(),
                _adaptiveSystemRows, _sourceWorkspace.cols(),
                true /* Use transpose */,
                -1.0, 1.0 /* Subtract from residual */ );
}

//////////////////////////////////////////////////////////////////////
// Assuming that pulses are emitted from a source at position
// sourcePos at interval timeScale, and with support support * timeScale,
// and given a point on the mesh, figure out the index range of
// pulses that must be considered in order to consider all contributions
// to that point.
//////////////////////////////////////////////////////////////////////
void SourceSolver::pulseIndexRange( const VEC3F &sourcePoint,
                                    const VEC3F &meshPoint,
                                    Real timeScale, int S, Real T, Real c,
                                    int &pulseMin, int &pulseMax )
{
  Real                       r = norm( meshPoint - sourcePoint );
  Real                       supp = timeScale * (Real)S;

#if 0
  // FIXME: Problem here?

  // Smallest index whose influence has not yet passed the
  // receiver point
  pulseMin = (int)floor( ( -T - ( r + supp ) / c ) / timeScale );

  // Largest index whose influence has reached the receiver point
  pulseMax = (int)ceil( ( T - ( r - supp ) / c ) / timeScale );
#endif
#if 0
  pulseMin = (int)floor( ( -supp - T - r / c ) / timeScale );

  pulseMax = (int)ceil( ( T + supp - r / c ) / timeScale );
#endif
  //pulseMin = (int)floor( ( -1.0 * T - 2.0 * supp - r / c ) / timeScale );
  pulseMin = (int)floor( ( -1.0 * T - r / c ) / timeScale );
  //pulseMax = (int)floor( ( T + 2.0 * supp - r / c ) / timeScale );
  pulseMax = (int)floor( ( T - r / c ) / timeScale );
}

//////////////////////////////////////////////////////////////////////
// Similar to the above, but figures out the index range if we want
// to look at a single source point interacting with *all* points
// on the mesh.
//////////////////////////////////////////////////////////////////////
void SourceSolver::pulseIndexRange( const VEC3F &sourcePoint,
                                    const Vector3Array &meshPoints,
                                    Real timeScale, int S, Real T, Real c,
                                    int &pulseMin, int &pulseMax )
{
  Real                       rMin, rMax;
  Real                       supp = timeScale * (Real)S;

  distanceBounds( meshPoints, sourcePoint, rMin, rMax );

  TRACE_ASSERT( rMin < rMax, "Distance bounds are invalid" );

#if 0
  // FIXME: Problem here?

  // Smallest index whose influence has not yet moved outside
  // of the mesh (ie. has not exceeded rMax)
  //pulseMin = (int)floor( ( -T - ( rMax + supp ) / c ) / timeScale );
  pulseMin = (int)floor( ( -T - ( rMax + supp ) / c ) / timeScale );

  // Largest index whose influence has reached the closest point
  // on the mesh (ie. has reached rMin)
  //pulseMax = (int)ceil( ( T - ( rMin - supp ) / c ) / timeScale );
  pulseMax = (int)ceil( ( T - ( rMin - supp ) / c ) / timeScale );
#endif
#if 0
  pulseMin = (int)floor( ( -supp - T - rMax / c ) / timeScale );

  pulseMax = (int)ceil( ( T + supp - rMin / c ) / timeScale );
#endif
  //pulseMin = (int)floor( ( -1.0 * T - 2.0 * supp - rMin / c ) / timeScale );
  pulseMin = (int)floor( ( -1.0 * T - rMin / c ) / timeScale );
  //pulseMax = (int)floor( ( T + 2.0 * supp - rMax / c ) / timeScale );
  pulseMax = (int)floor( ( T  - rMax / c ) / timeScale );

  // FIXME
  Real emissionTime = timeScale * (Real)pulseMin;
  Real tMin = emissionTime - supp + rMin / c;
  Real tMax = emissionTime + supp + rMax / c;

  TRACE_ASSERT( tMin >= -1.0 * ( T + 2.0 * supp ), "Min index out of bounds" );
  TRACE_ASSERT( tMax <= ( T + 2.0 * supp ), "Min index out of bounds" );

  emissionTime = timeScale * (Real)pulseMax;
  tMin = emissionTime - supp + rMin / c;
  tMax = emissionTime + supp + rMax / c;

  TRACE_ASSERT( tMin >= -1.0 * ( T + 2.0 * supp ), "Max index out of bounds" );
  TRACE_ASSERT( tMax <= ( T + 2.0 * supp ), "Max index out of bounds" );
}

//////////////////////////////////////////////////////////////////////
// Unit acceleration vector in the appropriate direction
//////////////////////////////////////////////////////////////////////
VEC3F SourceSolver::getAccelerationVector( Direction direction )
{
  switch ( direction )
  {
    case TRANS_X:
    default:
    {
      return VEC3F( 1.0, 0.0, 0.0 );
    }
    case TRANS_Y:
    {
      return VEC3F( 0.0, 1.0, 0.0 );
    }
    case TRANS_Z:
    {
      return VEC3F( 0.0, 0.0, 1.0 );
    }
    case ROT_X:
    {
      return VEC3F( 1.0, 0.0, 0.0 );
    }
    case ROT_Y:
    {
      return VEC3F( 0.0, 1.0, 0.0 );
    }
    case ROT_Z:
    {
      return VEC3F( 0.0, 0.0, 1.0 );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Gets the directional acceleration term at a point on the mesh
//////////////////////////////////////////////////////////////////////
Real SourceSolver::getAccelerationTerm( Direction direction,
                                        const VEC3F &accelerationVector,
                                        const VEC3F &position,
                                        const VEC3F &normal )
{
  if ( direction <= TRANS_Z )
  {
    return accelerationVector * normal;
  }
  else
  {
#if 0
    VEC3F crossTerm = cross( position - _mesh.getCenterOfMass(), normal );
    cout << "Got cross term " << crossTerm << endl;
#endif
    return accelerationVector * cross( position - _mesh.getCenterOfMass(),
                                       normal );
  }
}

//////////////////////////////////////////////////////////////////////
// Find closest and farthest points relative to the mesh's
// center of mass
//////////////////////////////////////////////////////////////////////
void SourceSolver::findPositionRange()
{
  Real             r;
  const VEC3F     &x_cm = _mesh.getCenterOfMass();

  _rMin = FLT_MAX;
  _rMax = FLT_MIN;

  for ( int i = 0; i < _mesh.vertices().size(); i++ )
  {
    r = norm( _mesh.vertices()[ i ] - x_cm );

    _rMin = min( r, _rMin );
    _rMax = max( r, _rMax );
  }
}

//////////////////////////////////////////////////////////////////////
// Given pulse time scale timeScale, support width S, acceleration pulse
// interval [-T,T] and speed of sound c, computes the minimum
// and maximum indices for pulses which can lead to acceleration
// boundary conditions.
//////////////////////////////////////////////////////////////////////
void SourceSolver::pulseIndexRange( Real timeScale, int S, Real T, Real c,
                                    int &pulseMin, int &pulseMax )
{
  pulseMin = (int)( -1.0 * ( _rMax / c + T ) / timeScale ) - S;
  pulseMax = (int)( ( T - _rMin / c ) / timeScale ) + S;

  TRACE_ASSERT( pulseMin <= pulseMax, "Invalid pulse range" );
}

//////////////////////////////////////////////////////////////////////
// In this version, we can specify a distinct interval between
// pulses (other then the time scale)
//////////////////////////////////////////////////////////////////////
void SourceSolver::pulseIndexRange( Real timeScale, Real pulseInterval,
                                    int S, Real T, Real c,
                                    int &pulseMin, int &pulseMax )
{
  int boundingPulses;

  // Figure out how many pulses intervals are needed to accomodate
  // the width of the pulse
  boundingPulses = (int)ceil( (Real)S * timeScale / pulseInterval );

  pulseMin = (int)( -1.0 * ( _rMax / c + T ) / pulseInterval ) - boundingPulses;
  pulseMax = (int)( ( T - _rMin / c ) / pulseInterval ) + boundingPulses;

  TRACE_ASSERT( pulseMin <= pulseMax, "Invalid pulse range" );
}

//////////////////////////////////////////////////////////////////////
// For a given pulse, figure out the range of times over which it
// interacts with the mesh
//////////////////////////////////////////////////////////////////////
void SourceSolver::pulseTimeRange( const VEC3F &sourcePoint,
                                   const SourceIndex &idx,
                                   const Vector3Array &meshPoints,
                                   Real timeScale, int S, Real c,
                                   Real &tMin, Real &tMax )
{
  Real                       rMin, rMax;
  Real                       emissionTime;
  Real                       supp = timeScale * (Real)S;

  emissionTime = timeScale * (Real)idx._timeIdx;
  //emissionTime = timeScale * (Real)idx._timeIdx / c;

  distanceBounds( meshPoints, sourcePoint, rMin, rMax );

#if 0
  cout << "In pulseTimeRange" << endl;
  cout << SDUMP( emissionTime ) << endl;
  cout << SDUMP( supp ) << endl;
  cout << SDUMP( rMin ) << "   " << SDUMP( rMax ) << endl;
#endif

#if 0
  // Pulse first reaches rMin:
  tMin = emissionTime + ( rMin - supp ) / c;

  // Pulse last interacts with rMax:
  tMax = emissionTime + ( rMax + supp ) / c;
#endif
  // Pulse first reaches rMin
  tMin = emissionTime - supp + rMin / c;

  // Pulse last interacts with rMax
  tMax = emissionTime + supp + rMax / c;
}

//////////////////////////////////////////////////////////////////////
// Same as the above, but only for a specific receiver point
//////////////////////////////////////////////////////////////////////
void SourceSolver::pulseTimeRange( const VEC3F &sourcePoint,
                                   const SourceIndex &idx,
                                   const VEC3F &meshPoint,
                                   Real timeScale, int S, Real c,
                                   Real &tMin, Real &tMax )
{
  Real                       r;
  Real                       emissionTime;
  Real                       supp = timeScale * (Real)S;

  emissionTime = timeScale * (Real)idx._timeIdx;
  //emissionTime = timeScale * (Real)idx._timeIdx / c;

  r = norm( meshPoint - sourcePoint );

#if 0
  // Pulse first reaches rMin:
  tMin = emissionTime + ( r - supp ) / c;

  // Pulse last interacts with rMax:
  tMax = emissionTime + ( r + supp ) / c;
#endif
  // Pulse first reaches rMin
  tMin = emissionTime - supp + r / c;

  // Pulse last interacts with rMax
  tMax = emissionTime + supp + r / c;
}

//////////////////////////////////////////////////////////////////////
// Finds the closest and farthest points in the given list relative
// to the given source location
//////////////////////////////////////////////////////////////////////
void SourceSolver::distanceBounds( const Vector3Array &meshPoints,
                                   const VEC3F &sourcePoint,
                                   Real &rMin, Real &rMax )
{
  Real                       r;

  rMin = FLT_MAX;
  rMax = 0.0;

  for ( int point_idx = 0; point_idx < meshPoints.size(); point_idx++ )
  {
    r = norm( meshPoints[ point_idx ] - sourcePoint );

    rMin = min( r, rMin );
    rMax = max( r, rMax );
  }
}
