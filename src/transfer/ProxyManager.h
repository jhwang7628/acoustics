//////////////////////////////////////////////////////////////////////
// ProxyManager.h: Interface
//
//////////////////////////////////////////////////////////////////////

#ifndef PROXY_MANAGER_H
#define PROXY_MANAGER_H

#include <TYPES.h>

#include <geometry/RigidMesh.h>

#include <transfer/RadialApproximation.h>

#include <vector>

//////////////////////////////////////////////////////////////////////
// ProxyManager class
//
// Comments
//////////////////////////////////////////////////////////////////////
class ProxyManager {
    public:
        ProxyManager( std::vector<RigidMesh *> &rigidMeshes,
                      REAL proxyMinScale, int proxyIncrements,
                      std::string proxyPathPrefix,
                      std::string proxyFilePrefix,
                      bool useWavelets, bool useSymmetricProxies,
                      bool matchSamplingRate,
                      ObjectMeasure measure,
                      bool useDirectionSet,
                      const std::string &directionSetFile,
                      const std::string &directionSetSuffix );

        // Destructor
        virtual ~ProxyManager();

        // For a given object ID and listening positiong, return it's proxy
        // and scaling factor, and also rotate the listening position and
        // translational/rotational accelerations in to the proxy space
        void getProxyData( int objectID, CompactRadialApproximation * &proxy,
                           REAL &scale,
                           Point3d &listeningPosition,
                           Vector3d &translationalAcceleration,
                           Vector3d &rotationalAcceleration );

    protected:

    private:
        // Clips each rigid mesh to one of our proxy ellipsoids
        void fitEllipsoids( REAL proxyMinScale, int proxyIncrements,
                            std::string proxyPathPrefix,
                            std::string proxyFilePrefix,
                            bool useWavelets, bool useSymmetricProxies,
                            bool matchSamplingRate,
                            ObjectMeasure measure,
                            bool useDirectionSet,
                            const std::string &directionSetFile,
                            const std::string &directionSetSuffix );

        void initProxyStorage( int proxyIncrements );

        void loadProxy( int xIncrement, int yIncrement,
                        std::string proxyPathPrefix,
                        std::string proxyFilePrefix,
                        bool useWavelets, bool useSymmetricProxies,
                        bool matchSamplingRate,
                        bool useDirectionSet,
                        const std::string &directionSetFile,
                        const std::string &directionSetSuffix );

        static int FitProxyIncrement( REAL scale, REAL scaleDivision );

    private:
        // The rigid meshes to build proxies for
        std::vector<RigidMesh *>  &_rigidMeshes;

        std::vector<IntPair>       _proxyIDs;

        // All proxy approximations
        std::vector< std::vector<RadialApproximation::AccelerationSet *> >
            _proxyModels;

        std::vector< std::vector<CompactRadialApproximation *> >
            _compactProxyModels;

};

#endif
