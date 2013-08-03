//////////////////////////////////////////////////////////////////////
// apply_impulses.cpp: Using a set of scene objects and a file with
//                     impulse information, generates an acceleration
//                     noise signal
//
//////////////////////////////////////////////////////////////////////

#include <TYPES.h>

#include <datastructure/ThreadSpecificData.h>

#include <math/SampledFunction.h>

#include <mesh/ClosestPointMesh.h>
#include <mesh/GTS_TriMesh.h>
#include <mesh/RigidMesh.h>
#include <mesh/TriMesh.h>

#include <parser/Parser.h>

#include <transfer/AccelerationNoiseModel.h>
#include <transfer/ProxyManager.h>

#include <util/IO.h>
#include <util/MERSENNETWISTER.h>
#include <util/timer.h>
#include <util/trace.h>

#include <fstream>

using namespace std;

//////////////////////////////////////////////////////////////////////
// For impulse randomization
//////////////////////////////////////////////////////////////////////
bool                         RANDOMIZE_IMPULSES = false;
Real                         SIM_TIM_STEP = 0.001;
MERSENNETWISTER              GENERATOR;

//////////////////////////////////////////////////////////////////////
// For representing an impulse between two objects
//////////////////////////////////////////////////////////////////////
struct ObjectImpulseTriplet {
  int                                      _objA_id;
  int                                      _objB_id;

  AccelerationNoiseModel::TwoObjectImpact  _impactData;
};

//////////////////////////////////////////////////////////////////////
// For representing an impulse between an object and the ground
// plane
//////////////////////////////////////////////////////////////////////
struct PlaneImpulseTriplet {
  int                                      _objA_id;

  AccelerationNoiseModel::PlaneImpact      _impactData;
};

//////////////////////////////////////////////////////////////////////
// Attempts to read impulses
//////////////////////////////////////////////////////////////////////
bool ReadImpulses( const char *filename,
                   vector<ObjectImpulseTriplet> &objectImpulses,
                   vector<PlaneImpulseTriplet> &planeImpulses,
                   Real &maxObjectImpulse, Real &maxPlaneImpulse );
                   
//////////////////////////////////////////////////////////////////////
// Writes impulses to binary files
//////////////////////////////////////////////////////////////////////
void WriteImpulses( const char *filename,
                    vector<ObjectImpulseTriplet> &objectImpulses,
                    vector<PlaneImpulseTriplet> &planeImpulses );

//////////////////////////////////////////////////////////////////////
// Samples impulses from a file stream
//////////////////////////////////////////////////////////////////////
void GatherImpulses( ifstream &fin,
                     vector<ObjectImpulseTriplet> &objectImpulses,
                     vector<PlaneImpulseTriplet> &planeImpulses,
                     Real &maxObjectImpulse,
                     Real &maxPlaneImpulse,
                     // Default plane material parameters are the
                     // material parameters of concrete
                     Real planeYoungsModulus = 3.0e10,
                     Real planePoissonRatio = 0.2 );

//////////////////////////////////////////////////////////////////////
// Samples an inter-object impulse from a file stream
//////////////////////////////////////////////////////////////////////
void ReadImpulse( ifstream &fin, ObjectImpulseTriplet &objectImpulse );

//////////////////////////////////////////////////////////////////////
// Samples a plane-object impulse from a file stream
//////////////////////////////////////////////////////////////////////
void ReadImpulse( ifstream &fin, PlaneImpulseTriplet &planeImpulse );

//////////////////////////////////////////////////////////////////////
// Adds together the signals for a bunch of sampled functions
//////////////////////////////////////////////////////////////////////
void AddThreadSignals( ThreadSpecificData<SampledFunction> &pressures,
                       FloatArray &output );

//////////////////////////////////////////////////////////////////////
// Prunes impulses below a certain threshold
//
// Returns a list of active impulse indices, and fills in the
// pruned spots in the impulse length and scales array
//////////////////////////////////////////////////////////////////////
void PruneImpulses( const vector<PlaneImpulseTriplet> &planeImpulses,
                    const vector<ObjectImpulseTriplet> &objectImpulses,
                    bool truncate, Real planeThreshold, Real objectThreshold,
                    FloatArray &planeImpulseLengths,
                    FloatArray &planeImpulseScales,
                    FloatArray &objectImpulseLengths,
                    FloatArray &objectImpulseScales,
                    IntArray &activePlaneImpulses,
                    IntArray &activeObjectImpulses );

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
int main( int argc, char **argv )
{
  string                         configFileName( "default.xml" );
  string                         impulseFileName;
  Parser                        *parser = NULL;
  Parser::SceneObjectParms       parms;

  vector<ObjectImpulseTriplet>   objectImpulses;
  vector<PlaneImpulseTriplet>    planeImpulses;

  Timer                          synthTimer( "synthesis" );

#if 0
  SampledFunction                pressure( 192000 );
  SampledFunction                interObjectPressure( 192000 );
#endif

  ThreadSpecificData<SampledFunction> threadPressures(
                                            SampledFunction( 192000 ) );
  ThreadSpecificData<SampledFunction> threadInterObjectPressures(
                                            SampledFunction( 192000 ) );

  VEC3F                          listeningPosition( 1.0, 1.8, 1.0 );
  VEC3F                          imageListeningPosition( 1.0, -1.8, 1.0 );

  Real                           groundYoungsModulus = 9.0e9;
  Real                           groundPoissonRatio = 0.2;

  ofstream                      *planeTimeScaleFile = NULL;
  ofstream                      *objectTimeScaleFile = NULL;

  ProxyManager                  *proxyManager = NULL;

  bool                           useImage = false;
  bool                           noPlane = false;
  bool                           truncate = true;

  bool                           useCompressedFields = false;

  if ( argc < 2 )
  {
    cerr << "Usage: " << argv[ 0 ] << " <impulse file> "
         << "[<impulse output prefix>]" << endl;

    return 1;
  }

  impulseFileName = argv[ 1 ];

  if ( argc >= 3 )
  {
    useCompressedFields = ( atoi( argv[ 2 ] ) != 0 );

    if ( useCompressedFields ) {
      cout << "Using compressed fields" << endl << endl;
    }
  }

  if ( argc >= 4 )
  {
    char                     buf[ 1024 ];

    sprintf( buf, "%s_plane.txt", argv[ 3 ] );
    planeTimeScaleFile = new ofstream( buf );

    sprintf( buf, "%s_object.txt", argv[ 3 ] );
    objectTimeScaleFile = new ofstream( buf );
  }

  if ( argc >= 5 )
  {
    configFileName = argv[ 4 ];
  }

  parser = Parser::buildParser( configFileName );

  if ( !parser )
  {
    cerr << "ERROR: Could not build parser from " << configFileName << endl;
    return 1;
  }

  parms = parser->getSceneObjectParms( useCompressedFields );

  if ( parms._randomizeImpulseTimes ) {
    RANDOMIZE_IMPULSES = true;
    SIM_TIM_STEP = parms._simTimeStep;
  }

  ifstream                       impulseStream( impulseFileName.c_str() );

  if ( parms._useProxies ) {
    proxyManager = new ProxyManager( parms._rigidMeshes,
                                     parms._proxyMinScale,
                                     parms._proxyIncrements,
                                     parms._proxyPathPrefix,
                                     parms._useSymmetricProxies ?
                                        parms._symmetricProxyPrefix
                                      : parms._proxyFilePrefix,
                                     parms._useSymmetricProxies,
                                     parms._useSymmetricProxies,
                                     parms._matchWaveletSamplingRates,
                                     parms._proxyMeasure,
                                     parms._useDirectionSetProxies,
                                     parms._proxyDirectionSetFile,
                                     parms._proxyDirectionSetSuffix );
  }

  if ( impulseStream.fail() )
  {
    cerr << "Failed to open " << impulseFileName << endl;
    return 1;
  }

  cout << SDUMP( parms._collisionTimeScale ) << endl;

  Real                           maxObjectImpulse, maxPlaneImpulse;

  printf( "Reading impulses\n" );
  if ( !ReadImpulses( impulseFileName.c_str(),
                      objectImpulses, planeImpulses,
                      maxObjectImpulse, maxPlaneImpulse ) ) {
    GatherImpulses( impulseStream, objectImpulses, planeImpulses,
                    maxObjectImpulse, maxPlaneImpulse );

    WriteImpulses( impulseFileName.c_str(), objectImpulses, planeImpulses );
  }
                  //groundYoungsModulus, groundPoissonRatio );
  printf( "Done\n" );

  impulseStream.close();

  cout << SDUMP( objectImpulses.size() ) << endl;
  cout << SDUMP( planeImpulses.size() ) << endl;

  cout << SDUMP( maxObjectImpulse ) << endl;
  cout << SDUMP( maxPlaneImpulse ) << endl;

  Real maxImpulse = max( maxPlaneImpulse, maxObjectImpulse );

  FloatArray                 planeImpulseLengths( planeImpulses.size() );
  FloatArray                 planeImpulseScales( planeImpulses.size() );

  int numProcessed = 0;

  synthTimer.tick();
#pragma omp parallel for schedule(static) default(shared)
  for ( int i = 0; i < planeImpulses.size(); i++ )
  {
#if 0
    if ( i % 16 == 0 )
    {
      printf( "Impulse %08d of %08d\r", i, (int)planeImpulses.size() );
      fflush( stdout );
    }
#endif
#if 0
#pragma omp critical
    {
      numProcessed += 1;

      if ( numProcessed % 64 == 0 )
      {
        printf( "Impulse %08d of %08d\r", numProcessed,
                (int)planeImpulses.size() );
        fflush( stdout );
      }
    }
#endif
#if 0
#pragma omp parallel
  while ( numProcessed < planeImpulses.size() )
  {
    int                              i;
#pragma omp critical (getNextPlaneIndex)
    {
      i = numProcessed;
      numProcessed++;
    }
#endif

    PlaneImpulseTriplet             &hit = planeImpulses[ i ];
    bool                             addedSound;

    TRACE_ASSERT( hit._impactData._impulseMagnitude >= 0.0 );

    // Truncate impulses with a relative speed of less then 1cm/s
    if ( truncate && abs( hit._impactData._relativeSpeed ) < 3e-2 )
    //if ( truncate && abs( hit._impactData._relativeSpeed ) < 0.2 * 3e-2 )
    //if ( truncate && abs( hit._impactData._relativeSpeed ) < 3e-3 )
    //if ( hit._impactData._impulseMagnitude < 0.00005 * maxImpulse )
    //if ( hit._impactData._impulseMagnitude < 0.01 * maxPlaneImpulse )
    {
      planeImpulseLengths[ i ] = -1.0;
      planeImpulseScales[ i ] = -1.0;

      continue;
    }

    SampledFunction                 &pressure = threadPressures.get();

    int                              objectID = hit._objA_id;
    int                              meshID = parms._objectMeshIDs[ objectID ];

    AccelerationNoiseModel::MeshSet  objectA( *parms._meshes[ meshID ],
                                              *parms._curvatureMeshes[ meshID ],
                                              *parms._rigidMeshes[ meshID ],
                                              *parms._distanceMeshes[ meshID ],
                                              meshID, parms._useProxies );

    addedSound = AccelerationNoiseModel::AddImpactNoise(
                                            objectA, hit._impactData,
                                            listeningPosition,
                                            /*
                                            *parms._pulseModels[ meshID ],
                                            */
                                            parms._compactPulseModels[ meshID ],
                                            pressure,
                                            //planeTimeScaleFile );
                                            NULL,
                                            &planeImpulseLengths[ i ],
                                            &planeImpulseScales[ i ],
                                            proxyManager,
                                            parms._collisionTimeScale );

#if 0
    if ( !addedSound )
    {
      cout << "Truncating force scale" << endl;
      ( *planeTimeScaleFile ) << "-1.0 -1.0" << endl;
    }
#endif

    if ( useImage )
    {
      AccelerationNoiseModel::AddImpactNoise(
                                        objectA, hit._impactData,
                                        imageListeningPosition,
                                        /*
                                        *parms._pulseModels[ meshID ],
                                        */
                                        parms._compactPulseModels[ meshID ],
                                        pressure );
    }
  }
  printf( "\n" );

  FloatArray                 objectImpulseLengths( objectImpulses.size() );
  FloatArray                 objectImpulseScales( objectImpulses.size() );

  numProcessed = 0;

#pragma omp parallel for schedule(static) default(shared)
  for ( int i = 0; i < objectImpulses.size(); i++ )
  {
#if 0
    if ( i % 16 == 0 )
    {
      printf( "Impulse %08d of %08d\r", i, (int)objectImpulses.size() );
      fflush( stdout );
    }
#endif

#if 0
#pragma omp critical
    {
      numProcessed += 1;

      if ( numProcessed % 64 == 0 )
      {
        printf( "Impulse %08d of %08d\r", numProcessed,
                (int)objectImpulses.size() );
        fflush( stdout );
      }
    }
#endif
#if 0
#pragma omp parallel
  while ( numProcessed < objectImpulses.size() )
  {
    int                              i;
#pragma omp critical (getNextObjIndex)
    {
      i = numProcessed;
      numProcessed++;
    }
#endif

    ObjectImpulseTriplet            &hit = objectImpulses[ i ];
    int                              objectID_A = hit._objA_id;
    int                              objectID_B = hit._objB_id;
    int                              meshID_A;
    int                              meshID_B;
    bool                             addedSound;

    // Truncate impulses with a relative speed of less then 1cm/s
    if ( truncate && abs( hit._impactData._relativeSpeed ) < 5e-2 )
    //if ( truncate && abs( hit._impactData._relativeSpeed ) < 0.2 * 5e-2 )
    //if ( truncate && abs( hit._impactData._relativeSpeed ) < 5e-3 )
    //if ( hit._impactData._impulseMagnitude < 0.00005 * maxImpulse )
    //if ( hit._impactData._impulseMagnitude < 0.01 * maxObjectImpulse )
    {
      objectImpulseLengths[ i ] = -1.0;
      objectImpulseScales[ i ] = -1.0;

#if 0
      if ( objectTimeScaleFile )
      {
        ( *objectTimeScaleFile ) << "-1.0 -1.0" << endl;
      }
#endif

      continue;
    }

    SampledFunction                 &interObjectPressure
                                            = threadInterObjectPressures.get();
    
    meshID_A = parms._objectMeshIDs[ objectID_A ];
    meshID_B = parms._objectMeshIDs[ objectID_B ];

    AccelerationNoiseModel::MeshSet  objectA(
                                        *parms._meshes[ meshID_A ],
                                        *parms._curvatureMeshes[ meshID_A ],
                                        *parms._rigidMeshes[ meshID_A ],
                                        *parms._distanceMeshes[ meshID_A ],
                                        meshID_A, parms._useProxies );

    AccelerationNoiseModel::MeshSet  objectB(
                                        *parms._meshes[ meshID_B ],
                                        *parms._curvatureMeshes[ meshID_B ],
                                        *parms._rigidMeshes[ meshID_B ],
                                        *parms._distanceMeshes[ meshID_B ],
                                        meshID_B, parms._useProxies );

    addedSound = AccelerationNoiseModel::AddImpactNoise(
                                        objectA, objectB,
                                        hit._impactData,
                                        listeningPosition,
                                        /*
                                        *parms._pulseModels[ meshID_A ],
                                        *parms._pulseModels[ meshID_B ],
                                        */
                                        parms._compactPulseModels[ meshID_A ],
                                        parms._compactPulseModels[ meshID_B ],
                                        interObjectPressure,
                                        //objectTimeScaleFile );
                                        NULL,
                                        &objectImpulseLengths[ i ],
                                        &objectImpulseScales[ i ],
                                        proxyManager,
                                        parms._collisionTimeScale );

#if 0
    if ( !addedSound )
    {
      cout << "Truncating force scale" << endl;
      ( *objectTimeScaleFile ) << "-1.0 -1.0" << endl;
    }
#endif

    if ( useImage )
    {
      AccelerationNoiseModel::AddImpactNoise(
                                        objectA, objectB,
                                        hit._impactData,
                                        imageListeningPosition,
                                        /*
                                        *parms._pulseModels[ meshID_A ],
                                        *parms._pulseModels[ meshID_B ],
                                        */
                                        parms._compactPulseModels[ meshID_A ],
                                        parms._compactPulseModels[ meshID_B ],
                                        interObjectPressure );
    }
  }
  synthTimer.tock();

  // FIXME
  {
    printf( "Synthesis took %f s\n", synthTimer.getTotalSecs() );
    //abort();
  }

  int truncatedImpulses = 0;

  if ( planeTimeScaleFile )
  {
    for ( int i = 0; i < planeImpulseLengths.size(); i++ )
    {
      ( *planeTimeScaleFile ) << planeImpulseLengths[ i ] << ' '
                              << planeImpulseScales[ i ] << endl;
    }
  }

  if ( objectTimeScaleFile )
  {
    for ( int i = 0; i < objectImpulseLengths.size(); i++ )
    {
      if ( objectImpulseLengths[ i ] > 0.0 ) {
        truncatedImpulses += 1;
      }

      ( *objectTimeScaleFile ) << objectImpulseLengths[ i ] << ' '
                               << objectImpulseScales[ i ] << endl;
    }
  }

  for ( int i = 0; i < planeImpulseLengths.size(); i++ ) {
    if ( planeImpulseLengths[ i ] > 0.0 ) {
      truncatedImpulses += 1;
    }
  }

  for ( int i = 0; i < objectImpulseLengths.size(); i++ ) {
    if ( objectImpulseLengths[ i ] > 0.0 ) {
      truncatedImpulses += 1;
    }
  }

  FloatArray pressureOutput;
  FloatArray interObjectPressureOutput;

  printf( "Found %d truncated impulses\n", truncatedImpulses );

#if 0
  writeVector( "test.vector", pressure.data() );
  writeVector( "objectTest.vector", interObjectPressure.data() );
#endif

  cout << "Adding thread signals" << endl;

  AddThreadSignals( threadPressures, pressureOutput );
  AddThreadSignals( threadInterObjectPressures, interObjectPressureOutput );

  cout << "Writing" << endl;

  if ( useCompressedFields ) {
    writeVector( "test_compressed.vector", pressureOutput );
    writeVector( "planeTest_compressed.vector", pressureOutput );
    writeVector( "objectTest_compressed.vector", interObjectPressureOutput );
  } else if ( parms._useProxies ) {
    char                     fileName[ 1024 ];
    char                     planeFileName[ 1024 ];
    char                     objectFileName[ 1024 ];
    const char              *measureSuffix;
    const char              *compressionSuffix;
    const char              *directionSetSuffix;

    if ( parms._proxyMeasure == VOLUME ) {
      measureSuffix = "_V";
    } else if ( parms._proxyMeasure == MATCH_RIGID ) {
      measureSuffix = "_M";
    } else {
      measureSuffix = "_A";
    }

    if ( parms._useSymmetricProxies ) {
      compressionSuffix = "_comp";
    } else {
      compressionSuffix = "";
    }

    if ( parms._useDirectionSetProxies ) {
      directionSetSuffix = parms._proxyDirectionSetSuffix.c_str();
    } else {
      directionSetSuffix = "";
    }

    sprintf( fileName, "proxy%s%s%s_test.vector",
             measureSuffix, compressionSuffix, directionSetSuffix );
    sprintf( planeFileName, "proxy%s%s%s_planeTest.vector",
             measureSuffix, compressionSuffix, directionSetSuffix );
    sprintf( objectFileName, "proxy%s%s%s_objectTest.vector",
             measureSuffix, compressionSuffix, directionSetSuffix );

    writeVector( fileName, pressureOutput );
    writeVector( planeFileName, pressureOutput );
    writeVector( objectFileName, interObjectPressureOutput );

#if 0
    if ( parms._useSymmetricProxies ) {
      if ( parms._proxyMeasure == VOLUME ) {
        writeVector( "proxy_V_comp_test.vector", pressureOutput );
        writeVector( "proxy_V_comp_planeTest.vector", pressureOutput );
        writeVector( "proxy_V_comp_objectTest.vector",
                     interObjectPressureOutput );
      } else if ( parms._proxyMeasure == MATCH_RIGID ) {
        writeVector( "proxy_M_comp_test.vector", pressureOutput );
        writeVector( "proxy_M_comp_planeTest.vector", pressureOutput );
        writeVector( "proxy_M_comp_objectTest.vector",
                     interObjectPressureOutput );
      } else {
        writeVector( "proxy_A_comp_test.vector", pressureOutput );
        writeVector( "proxy_A_comp_planeTest.vector", pressureOutput );
        writeVector( "proxy_A_comp_objectTest.vector",
                     interObjectPressureOutput );
      }
    } else {
      if ( parms._proxyMeasure == VOLUME ) {
        writeVector( "proxy_V_test.vector", pressureOutput );
        writeVector( "proxy_V_planeTest.vector", pressureOutput );
        writeVector( "proxy_V_objectTest.vector", interObjectPressureOutput );
      } else if ( parms._proxyMeasure == MATCH_RIGID ) {
        writeVector( "proxy_M_test.vector", pressureOutput );
        writeVector( "proxy_M_planeTest.vector", pressureOutput );
        writeVector( "proxy_M_objectTest.vector", interObjectPressureOutput );
      } else {
        writeVector( "proxy_A_test.vector", pressureOutput );
        writeVector( "proxy_A_planeTest.vector", pressureOutput );
        writeVector( "proxy_A_objectTest.vector", interObjectPressureOutput );
      }
    }
#endif
  } else if ( RANDOMIZE_IMPULSES ) {
    writeVector( "randomized_test.vector", pressureOutput );
    writeVector( "randomized_planeTest.vector", pressureOutput );
    writeVector( "randomized_objectTest.vector", interObjectPressureOutput );
  } else {
    writeVector( "new_test.vector", pressureOutput );
    writeVector( "new_planeTest.vector", pressureOutput );
    writeVector( "new_objectTest.vector", interObjectPressureOutput );
  }

  parms.clear();

  delete planeTimeScaleFile;
  delete objectTimeScaleFile;

  return 0;
}

//////////////////////////////////////////////////////////////////////
// Attempts to read impulses
//////////////////////////////////////////////////////////////////////
bool ReadImpulses( const char *filename,
                   vector<ObjectImpulseTriplet> &objectImpulses,
                   vector<PlaneImpulseTriplet> &planeImpulses,
                   Real &maxObjectImpulse, Real &maxPlaneImpulse )
{
  char                       binaryFilename[ 1024 ];
  FILE                      *f;

  int                        objectSize;
  int                        planeSize;

  size_t                     bytesRead;

  maxObjectImpulse = 0.0;
  maxPlaneImpulse = 0.0;

  sprintf( binaryFilename, "%s.dat", filename );

  f = fopen( binaryFilename, "rb" );
  if ( !f ) {
    cerr << "Could not load " << binaryFilename << endl;
    return false;
  }

  bytesRead = fread( &objectSize, sizeof( int ), 1, f );
  objectImpulses.clear();
  objectImpulses.resize( objectSize );
  bytesRead = fread( objectImpulses.data(), sizeof( ObjectImpulseTriplet ),
                     objectSize, f );

  bytesRead = fread( &planeSize, sizeof( int ), 1, f );
  planeImpulses.clear();
  planeImpulses.resize( planeSize );
  bytesRead = fread( planeImpulses.data(), sizeof( PlaneImpulseTriplet ),
                     planeSize, f );

  fclose( f );

  for ( int i = 0; i < objectImpulses.size(); i++ ) {
    maxObjectImpulse = max( maxObjectImpulse,
                            objectImpulses[ i ]._impactData._impulseMagnitude );
  }

  for ( int i = 0; i < planeImpulses.size(); i++ ) {
    maxPlaneImpulse = max( maxPlaneImpulse,
                           planeImpulses[ i ]._impactData._impulseMagnitude );
  }

  return true;
}
                   
//////////////////////////////////////////////////////////////////////
// Writes impulses to binary files
//////////////////////////////////////////////////////////////////////
void WriteImpulses( const char *filename,
                    vector<ObjectImpulseTriplet> &objectImpulses,
                    vector<PlaneImpulseTriplet> &planeImpulses )
{
  char                       binaryFilename[ 1024 ];
  FILE                      *f;

  int                        objectSize = objectImpulses.size();
  int                        planeSize = planeImpulses.size();

  size_t                     bytesWritten;

  sprintf( binaryFilename, "%s.dat", filename );

  f = fopen( binaryFilename, "wb" );
  if ( !f ) {
    cerr << "Could not open " << binaryFilename << endl;
    abort();
  }

  bytesWritten = fwrite( &objectSize, sizeof( int ), 1, f );
  bytesWritten = fwrite( objectImpulses.data(), sizeof( ObjectImpulseTriplet ),
                         objectSize, f );

  bytesWritten = fwrite( &planeSize, sizeof( int ), 1, f );
  bytesWritten = fwrite( planeImpulses.data(), sizeof( PlaneImpulseTriplet ),
                         planeSize, f );

  fclose( f );
}

//////////////////////////////////////////////////////////////////////
// Samples inter-object impulses from a file stream
//////////////////////////////////////////////////////////////////////
void GatherImpulses( ifstream &fin,
                     vector<ObjectImpulseTriplet> &objectImpulses,
                     vector<PlaneImpulseTriplet> &planeImpulses,
                     Real &maxObjectImpulse,
                     Real &maxPlaneImpulse,
                     Real planeYoungsModulus,
                     Real planePoissonRatio )
{
  char                       impulseType;
  char                       lineBuffer[ 4096 ];

  ObjectImpulseTriplet       objectImpulse;
  PlaneImpulseTriplet        planeImpulse;

  maxObjectImpulse = 0.0;
  maxPlaneImpulse = 0.0;

  planeImpulse._impactData._planeYoungsModulus = planeYoungsModulus;
  planeImpulse._impactData._planePoissonRatio = planePoissonRatio;

  int impulseNum = 0;

  while ( !fin.eof() )
  {
    // Process the next impulse

    // Get the impulse type
    fin >> impulseType;

    if ( impulseType == 'P' )
    {
      ReadImpulse( fin, objectImpulse );

      maxObjectImpulse = max( maxObjectImpulse,
                              objectImpulse._impactData._impulseMagnitude );

      objectImpulses.push_back( objectImpulse );
    }
    else if ( impulseType == 'C' )
    {
      ReadImpulse( fin, planeImpulse );

      maxPlaneImpulse = max( maxPlaneImpulse,
                             planeImpulse._impactData._impulseMagnitude );

      planeImpulses.push_back( planeImpulse );
    }
    else
    {
      cout << SDUMP( impulseType ) << endl;
      TRACE_ASSERT( NULL, "Should never get here" );
    }

    impulseNum++;
#if 0
    if ( impulseType != 'P' )
    {
      TRACE_ASSERT( impulseType == 'C' );

      // Must have a plane impulse, so skip this line
      fin.getline( lineBuffer, 4096 );
    }
#endif
  }
}

//////////////////////////////////////////////////////////////////////
// Samples an inter-object impulse from a file stream
//////////////////////////////////////////////////////////////////////
void ReadImpulse( ifstream &fin, ObjectImpulseTriplet &impulse )
{
  fin >> impulse._impactData._impulseTime;

  if ( RANDOMIZE_IMPULSES ) {
    impulse._impactData._impulseTime
      += GENERATOR.rand( SIM_TIM_STEP ) - SIM_TIM_STEP / 2.0;
  }
  
  // Object A information
  fin >> impulse._objA_id;
  fin >> impulse._impactData._posA[ 0 ];
  fin >> impulse._impactData._posA[ 1 ];
  fin >> impulse._impactData._posA[ 2 ];
  fin >> impulse._impactData._centerOfMassA[ 0 ];
  fin >> impulse._impactData._centerOfMassA[ 1 ];
  fin >> impulse._impactData._centerOfMassA[ 2 ];
  fin >> impulse._impactData._inverseRotationA.w();
  fin >> impulse._impactData._inverseRotationA.x();
  fin >> impulse._impactData._inverseRotationA.y();
  fin >> impulse._impactData._inverseRotationA.z();

  // Object B information
  fin >> impulse._objB_id;
  fin >> impulse._impactData._posB[ 0 ];
  fin >> impulse._impactData._posB[ 1 ];
  fin >> impulse._impactData._posB[ 2 ];
  fin >> impulse._impactData._centerOfMassB[ 0 ];
  fin >> impulse._impactData._centerOfMassB[ 1 ];
  fin >> impulse._impactData._centerOfMassB[ 2 ];
  fin >> impulse._impactData._inverseRotationB.w();
  fin >> impulse._impactData._inverseRotationB.x();
  fin >> impulse._impactData._inverseRotationB.y();
  fin >> impulse._impactData._inverseRotationB.z();

  // Other impulse information
  fin >> impulse._impactData._relativeSpeed;
  fin >> impulse._impactData._impulseMagnitude;
  fin >> impulse._impactData._impulseDirection[ 0 ];
  fin >> impulse._impactData._impulseDirection[ 1 ];
  fin >> impulse._impactData._impulseDirection[ 2 ];
}

//////////////////////////////////////////////////////////////////////
// Samples a plane-object impulse from a file stream
//////////////////////////////////////////////////////////////////////
void ReadImpulse( ifstream &fin, PlaneImpulseTriplet &impulse )
{
  fin >> impulse._impactData._impulseTime;

  if ( RANDOMIZE_IMPULSES ) {
    impulse._impactData._impulseTime
      += GENERATOR.rand( SIM_TIM_STEP ) - SIM_TIM_STEP / 2.0;
  }
  
  // Object A information
  fin >> impulse._objA_id;
  fin >> impulse._impactData._posA[ 0 ];
  fin >> impulse._impactData._posA[ 1 ];
  fin >> impulse._impactData._posA[ 2 ];
  fin >> impulse._impactData._centerOfMassA[ 0 ];
  fin >> impulse._impactData._centerOfMassA[ 1 ];
  fin >> impulse._impactData._centerOfMassA[ 2 ];
  fin >> impulse._impactData._inverseRotationA.w();
  fin >> impulse._impactData._inverseRotationA.x();
  fin >> impulse._impactData._inverseRotationA.y();
  fin >> impulse._impactData._inverseRotationA.z();

  // Other impulse information
  fin >> impulse._impactData._relativeSpeed;
  fin >> impulse._impactData._impulseMagnitude;
  fin >> impulse._impactData._impulseDirection[ 0 ];
  fin >> impulse._impactData._impulseDirection[ 1 ];
  fin >> impulse._impactData._impulseDirection[ 2 ];
}

//////////////////////////////////////////////////////////////////////
// Adds together the signals for a bunch of sampled functions
//////////////////////////////////////////////////////////////////////
void AddThreadSignals( ThreadSpecificData<SampledFunction> &pressures,
                       FloatArray &output )
{
  const map<int, SampledFunction> &allPressures = pressures.getAll();

  for ( map<int, SampledFunction>::const_iterator itr = allPressures.begin();
        itr != allPressures.end(); itr++ )
  {
    if ( output.size() < itr->second.data().size() )
    {
      output.resize( itr->second.data().size(), 0.0 );
    }
    
    for ( int i = 0; i < itr->second.data().size(); i++ )
    {
      output[ i ] += itr->second.data()[ i ];
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Prunes impulses below a certain threshold
//
// Returns a list of active impulse indices, and fills in the
// pruned spots in the impulse length and scales array
//////////////////////////////////////////////////////////////////////
void PruneImpulses( const vector<PlaneImpulseTriplet> &planeImpulses,
                    const vector<ObjectImpulseTriplet> &objectImpulses,
                    bool truncate, Real planeThreshold, Real objectThreshold,
                    FloatArray &planeImpulseLengths,
                    FloatArray &planeImpulseScales,
                    FloatArray &objectImpulseLengths,
                    FloatArray &objectImpulseScales,
                    IntArray &activePlaneImpulses,
                    IntArray &activeObjectImpulses )
{
  activePlaneImpulses.clear();
  for ( int i = 0; i < planeImpulses.size(); i++ ) {
    const PlaneImpulseTriplet     &hit = planeImpulses[ i ];

    // Truncate impulses with a relative speed of less then 1cm/s
    if ( truncate && abs( hit._impactData._relativeSpeed ) < planeThreshold ) {
      planeImpulseLengths[ i ] = -1.0;
      planeImpulseScales[ i ] = -1.0;
    } else {
      activePlaneImpulses.push_back( i );
    }
  }

  activeObjectImpulses.clear();
  for ( int i = 0; i < objectImpulses.size(); i++ ) {
    const ObjectImpulseTriplet    &hit = objectImpulses[ i ];

    if ( truncate && abs( hit._impactData._relativeSpeed ) < objectThreshold ) {
      objectImpulseLengths[ i ] = -1.0;
      objectImpulseScales[ i ] = -1.0;
    } else {
      activeObjectImpulses.push_back( i );
    }
  }
}
