//////////////////////////////////////////////////////////////////////
// ProxyManager.cpp: Implementation
//
//////////////////////////////////////////////////////////////////////

#include "ProxyManager.h"

#include <linearalgebra/Matrix3.hpp>

#include <transfer/CompressedMultiTermApproximation.h>
#include <transfer/MultiTermApproximation.h>

#include <utils/MathUtil.h>
#include <utils/trace.h>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
ProxyManager::ProxyManager( vector<RigidMesh *> &rigidMeshes,
                            REAL proxyMinScale, int proxyIncrements,
                            string proxyPathPrefix, string proxyFilePrefix,
                            bool useWavelets, bool useSymmetricProxies,
                            bool matchSamplingRate,
                            ObjectMeasure measure,
                            bool useDirectionSet,
                            const string &directionSetFile,
                            const string &directionSetSuffix )
    : _rigidMeshes( rigidMeshes )
{
    // Build a proxy ellipsoid for each mesh
    for ( int mesh_idx = 0; mesh_idx < _rigidMeshes.size(); mesh_idx++ ) {
        _rigidMeshes[ mesh_idx ]->fitEllipsoid( measure );
    }

    fitEllipsoids( proxyMinScale, proxyIncrements,
                   proxyPathPrefix, proxyFilePrefix,
                   useWavelets, useSymmetricProxies,
                   matchSamplingRate,
                   measure,
                   useDirectionSet,
                   directionSetFile,
                   directionSetSuffix );
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
ProxyManager::~ProxyManager()
{
    for ( int i = 0; i < _proxyModels.size(); i++ )
        for ( int j = 0; j < _proxyModels[ i ].size(); j++ ) {
            _proxyModels[ i ][ j ]->clear();

            delete _compactProxyModels[ i ][ j ];
        }
}

//////////////////////////////////////////////////////////////////////
// For a given object ID and listening positiong, return it's proxy
// and scaling factor, and also rotate the listening position and
// translational/rotational accelerations in to the proxy space
//////////////////////////////////////////////////////////////////////
void ProxyManager::getProxyData( int objectID,
                                 CompactRadialApproximation * &proxy,
                                 REAL &scale,
                                 Point3d &listeningPosition,
                                 Vector3d &translationalAcceleration,
                                 Vector3d &rotationalAcceleration )
{
    const Matrix3d            &axes = _rigidMeshes[ objectID ]->inertiaAxes();
    IntPair                    proxyID = _proxyIDs[ objectID ];

#if 0
    listeningPosition
        = listeningPosition - _rigidMeshes[ objectID ]->centerOfMass();
#endif

    listeningPosition = axes.transpose() * listeningPosition;

    // Swap coordinates, since the z-axis is the major axis, and the
    // x-axis is the minor axis
    MathUtil::swap( listeningPosition[ 0 ], listeningPosition[ 2 ] );

#if 0
    listeningPosition += _rigidMeshes[ objectID ]->centerOfMass();
#endif

    scale = _rigidMeshes[ objectID ]->ellipseScale();

    translationalAcceleration = axes * translationalAcceleration;
    MathUtil::swap( translationalAcceleration[ 0 ],
            translationalAcceleration[ 2 ] );

    rotationalAcceleration = axes * rotationalAcceleration;
    MathUtil::swap( rotationalAcceleration[ 0 ],
            rotationalAcceleration[ 2 ] );

    proxy = _compactProxyModels[ proxyID.first ][ proxyID.second ];

    TRACE_ASSERT( proxy != NULL );
}

//////////////////////////////////////////////////////////////////////
// Clips each rigid mesh to one of our proxy ellipsoids
//////////////////////////////////////////////////////////////////////
void ProxyManager::fitEllipsoids( REAL proxyMinScale, int proxyIncrements,
                                  string proxyPathPrefix,
                                  string proxyFilePrefix,
                                  bool useWavelets, bool useSymmetricProxies,
                                  bool matchSamplingRate,
                                  ObjectMeasure measure,
                                  bool useDirectionSet,
                                  const string &directionSetFile,
                                  const string &directionSetSuffix )
{
    REAL                       scaleDivision;
    REAL                       xScale;
    REAL                       yScale;

    int                        xIncrement;
    int                        yIncrement;

    initProxyStorage( proxyIncrements );

    scaleDivision = ( 1.0 - proxyMinScale ) / (REAL)( proxyIncrements );

    for ( int mesh_idx = 0; mesh_idx < _rigidMeshes.size(); mesh_idx++ ) {
        xScale = _rigidMeshes[ mesh_idx ]->ellipseScaleX();
        yScale = _rigidMeshes[ mesh_idx ]->ellipseScaleY();

        xIncrement = FitProxyIncrement( xScale, scaleDivision );
        yIncrement = FitProxyIncrement( yScale, scaleDivision );

        xIncrement = min( proxyIncrements, xIncrement );
        yIncrement = min( proxyIncrements, yIncrement );

        xIncrement = max( 0, xIncrement );
        yIncrement = max( 0, yIncrement );

        // FIXME: debugging
        xIncrement = 0;
        yIncrement = 0;

        TRACE_ASSERT( xIncrement >= yIncrement );

        // Figure out new scales for the ellipsoid
        xScale = 1.0 - scaleDivision * (REAL)xIncrement;
        yScale = 1.0 - scaleDivision * (REAL)yIncrement;

        _rigidMeshes[ mesh_idx ]->setEllipsoid( xScale, yScale, measure );

        printf( "Fitting mesh %d to ellipse (%d, %d)\n",
                mesh_idx, xIncrement, yIncrement );
        printf( "Ellipse scale is %e\n", _rigidMeshes[ mesh_idx ]->ellipseScale() );

        loadProxy( xIncrement, yIncrement, proxyPathPrefix, proxyFilePrefix,
                useWavelets, useSymmetricProxies, matchSamplingRate,
                useDirectionSet, directionSetFile, directionSetSuffix );

        _proxyIDs.push_back( IntPair( xIncrement, yIncrement ) );
    }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void ProxyManager::initProxyStorage( int proxyIncrements )
{
    _proxyModels.resize( proxyIncrements + 1 );
    _compactProxyModels.resize( proxyIncrements + 1 );

    for ( int i = 0; i < proxyIncrements + 1; i++ ) {
        _proxyModels[ i ].resize( i + 1 );
        _compactProxyModels[ i ].resize( i + 1 );

        for ( int j = 0; j < _proxyModels[ i ].size(); j++ ) {
            _proxyModels[ i ][ j ] = new RadialApproximation::AccelerationSet();
            _compactProxyModels[ i ][ j ] = NULL;
        }
    }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void ProxyManager::loadProxy( int xIncrement, int yIncrement,
                              std::string proxyPathPrefix,
                              std::string proxyFilePrefix,
                              bool useWavelets, bool useSymmetricProxies,
                              bool matchSamplingRate,
                              bool useDirectionSet,
                              const std::string &directionSetFile,
                              const std::string &directionSetSuffix )
{
    char                       buf[ 1024 ];

    TRACE_ASSERT( xIncrement < _proxyModels.size()
            && xIncrement < _compactProxyModels.size() );
    TRACE_ASSERT( yIncrement < _proxyModels[ xIncrement ].size()
            && yIncrement < _compactProxyModels[ xIncrement ].size() );

    if ( _compactProxyModels[ xIncrement ][ yIncrement ] != NULL ) {
        return;
    }

    // Get the proxy file path
    sprintf( buf, "%s_%03d_%03d/%s", proxyPathPrefix.c_str(),
            xIncrement, yIncrement, proxyFilePrefix.c_str() );

    // FIXME: get wavelet stuff in here
    if ( useWavelets ) {
        CompressedMultiTermApproximation::ReadAccelerationSet(
                buf, *_proxyModels[ xIncrement ][ yIncrement ],
                useSymmetricProxies, matchSamplingRate,
                useDirectionSet ? &directionSetSuffix : NULL );
    } else {
        MultiTermApproximation::ReadAccelerationSet(
                buf, *_proxyModels[ xIncrement ][ yIncrement ] );
    }

    // Get the compact approximation
    _compactProxyModels[ xIncrement ][ yIncrement ]
        = _proxyModels[ xIncrement ][ yIncrement ]->_allFields[ 0 ]
        ->buildCompactApproximation(
                *_proxyModels[ xIncrement ][ yIncrement ],
                useDirectionSet ? &directionSetFile : NULL );
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
int ProxyManager::FitProxyIncrement( REAL scale, REAL scaleDivision )
{
    REAL                       lowerBound;
    REAL                       upperBound;
    REAL                       realIncrement;

    scale = 1.0 - scale;

    if ( scale >= 1.0 ) {
        return (int)( 1.0 / scaleDivision );
    } else if ( scale <= 0 ) {
        return 0;
    }

    realIncrement = scale / scaleDivision;
    lowerBound = floor( realIncrement );
    upperBound = ceil( realIncrement );

    if ( abs( lowerBound - realIncrement ) < abs( upperBound - realIncrement ) ) {
        return (int)lowerBound;
    } else {
        return (int)upperBound;
    }
}
