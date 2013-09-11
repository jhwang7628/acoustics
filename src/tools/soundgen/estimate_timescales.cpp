//////////////////////////////////////////////////////////////////////
// estimate_timescales.cpp: Using a set of scene objects and a file
//                          with impulse information, dumps out a file
//                          containing Hertz time scale estimates as
//                          well as appropriate scalings to apply to
//                          impulse forces, assuming forces are
//                          integrated over half sine pulses with the
//                          given time scales
//
//////////////////////////////////////////////////////////////////////

#include <TYPES.h>

#include <geometry/ClosestPointMesh.h>
#include <geometry/GTS_TriMesh.h>
#include <geometry/RigidMesh.h>
#include <geometry/TriangleMesh.hpp>

#include <parser/Parser.h>

#include <transfer/AccelerationNoiseModel.h>
#include <transfer/ProxyManager.h>

#include <utils/IO.h>
#include <utils/MERSENNETWISTER.h>
#include <utils/timer.hpp>
#include <utils/trace.h>

#include <fstream>

using namespace std;

//////////////////////////////////////////////////////////////////////
// For impulse randomization
//////////////////////////////////////////////////////////////////////
bool                         RANDOMIZE_IMPULSES = false;
REAL                         SIM_TIME_STEP = 0.001;
MERSENNETWISTER              GENERATOR;

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
int main( int argc, char **argv )
{
    string                         configFileName( "default.xml" );
    string                         impulseFileName;
    Parser                        *parser = NULL;
    Parser::SceneObjectParms       parms;

    REAL                           groundYoungsModulus = 9.0e9;
    REAL                           groundPoissonRatio = 0.2;

    ofstream                      *planeTimeScaleFile = NULL;
    ofstream                      *objectTimeScaleFile = NULL;

    bool                           truncate = true;

    if ( argc < 3 )
    {
        cerr << "Usage: " << argv[ 0 ] << " <impulse file> "
            << "[<impulse output prefix>]" << endl;

        return 1;
    }

    impulseFileName = argv[ 1 ];

    char                     buf[ 1024 ];

    sprintf( buf, "%s_plane.txt", argv[ 2 ] );
    planeTimeScaleFile = new ofstream( buf );

    sprintf( buf, "%s_object.txt", argv[ 2 ] );
    objectTimeScaleFile = new ofstream( buf );

    if ( argc >= 4 )
    {
        configFileName = argv[ 3 ];
    }

    parser = Parser::buildParser( configFileName );

    if ( !parser )
    {
        cerr << "ERROR: Could not build parser from " << configFileName << endl;
        return 1;
    }

    parms = parser->getSceneObjectParms( false /* don't load PAN data */,
                                         false, false /* Other PAN booleans (not used) */);

    if ( parms._randomizeImpulseTimes ) {
        RANDOMIZE_IMPULSES = true;
        SIM_TIME_STEP = parms._simTimeStep;
    }

    cout << SDUMP( parms._collisionTimeScale ) << endl;

    printf( "Reading impulses\n" );

    // Read impulse data from disk
    ImpulseIO::ImpulseSet          impulses = ImpulseIO::ReadImpulses( impulseFileName.c_str(),
                                                                       parms._randomizeImpulseTimes,
                                                                       parms._simTimeStep,
                                                                       groundYoungsModulus,
                                                                       groundPoissonRatio,
                                                                       true /* Try to read
                                                                               from binary */ );

    printf( "Done\n" );

    cout << SDUMP( impulses._objectImpulses.size() ) << endl;
    cout << SDUMP( impulses._planeImpulses.size() ) << endl;

    cout << SDUMP( impulses._maxObjectImpulse ) << endl;
    cout << SDUMP( impulses._maxPlaneImpulse ) << endl;

    REAL maxImpulse = max( impulses._maxPlaneImpulse, impulses._maxObjectImpulse );

    FloatArray                 planeImpulseLengths( impulses._planeImpulses.size() );
    FloatArray                 planeImpulseScales( impulses._planeImpulses.size() );

    for ( int i = 0; i < impulses._planeImpulses.size(); i++ )
    {

        ImpulseIO::PlaneImpulse         &hit = impulses._planeImpulses[ i ];

        TRACE_ASSERT( hit._impactData._impulseMagnitude >= 0.0 );

        // Truncate impulses with a relative speed of less then 1cm/s
        //
        // FIXME: This is awful
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

        int                              objectID = hit._objA_id;
        int                              meshID = parms._objectMeshIDs[ objectID ];

        AccelerationNoiseModel::MeshSet  objectA( *parms._meshes[ meshID ],
                                                  *parms._curvatureMeshes[ meshID ],
                                                  *parms._rigidMeshes[ meshID ],
                                                  *parms._distanceMeshes[ meshID ],
                                                  meshID, parms._useProxies );

        AccelerationNoiseModel::HertzImpactData impactData
                        = AccelerationNoiseModel::ImpactTimeScale( objectA, hit._impactData,
                                                                   parms._collisionTimeScale );

        planeImpulseLengths[ i ] = impactData.first;
        planeImpulseScales[ i ] = impactData.second;

        printf( "Impulse time scale = %e s, force scale = %e\n",
                impactData.first, impactData.second );
    }

    FloatArray                 objectImpulseLengths( impulses._objectImpulses.size() );
    FloatArray                 objectImpulseScales( impulses._objectImpulses.size() );

    for ( int i = 0; i < impulses._objectImpulses.size(); i++ )
    {
        ImpulseIO::ObjectImpulse        &hit = impulses._objectImpulses[ i ];
        int                              objectID_A = hit._objA_id;
        int                              objectID_B = hit._objB_id;
        int                              meshID_A;
        int                              meshID_B;

        // Truncate impulses with a relative speed of less then 1cm/s
        //
        // FIXME: this is awful
        if ( truncate && abs( hit._impactData._relativeSpeed ) < 5e-2 )
            //if ( truncate && abs( hit._impactData._relativeSpeed ) < 0.2 * 5e-2 )
            //if ( truncate && abs( hit._impactData._relativeSpeed ) < 5e-3 )
            //if ( hit._impactData._impulseMagnitude < 0.00005 * maxImpulse )
            //if ( hit._impactData._impulseMagnitude < 0.01 * maxObjectImpulse )
        {
            objectImpulseLengths[ i ] = -1.0;
            objectImpulseScales[ i ] = -1.0;

            continue;
        }

        meshID_A = parms._objectMeshIDs[ objectID_A ];
        meshID_B = parms._objectMeshIDs[ objectID_B ];

        AccelerationNoiseModel::MeshSet  objectA( *parms._meshes[ meshID_A ],
                                                  *parms._curvatureMeshes[ meshID_A ],
                                                  *parms._rigidMeshes[ meshID_A ],
                                                  *parms._distanceMeshes[ meshID_A ],
                                                  meshID_A, parms._useProxies );

        AccelerationNoiseModel::MeshSet  objectB( *parms._meshes[ meshID_B ],
                                                  *parms._curvatureMeshes[ meshID_B ],
                                                  *parms._rigidMeshes[ meshID_B ],
                                                  *parms._distanceMeshes[ meshID_B ],
                                                  meshID_B, parms._useProxies );

    }

    int truncatedImpulses = 0;

    // Write plane impulse time scale data
    if ( planeTimeScaleFile )
    {
        for ( int i = 0; i < planeImpulseLengths.size(); i++ )
        {
            ( *planeTimeScaleFile ) << planeImpulseLengths[ i ] << ' '
                << planeImpulseScales[ i ] << endl;
        }
    }

    // Write object-object impulse time scale data
    if ( objectTimeScaleFile )
    {
        for ( int i = 0; i < objectImpulseLengths.size(); i++ )
        {
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

    printf( "Found %d truncated impulses\n", truncatedImpulses );

    parms.clear();

    delete planeTimeScaleFile;
    delete objectTimeScaleFile;

    return 0;
}
