#ifndef IMPULSE_SERIES_READER_H
#define IMPULSE_SERIES_READER_H
#include <string>
#include <vector>
#include <memory>
#include <config.h>
#include <linearalgebra/Vector3.hpp>
#include <geometry/TriangleMesh.hpp>
#include <wavesolver/FDTD_Objects.h>
#include <modal_model/ImpulseSeriesObject.h>

//##############################################################################
// IO Class that reads the text impulse data.
//
// Note:
//  Read the impulses created by rigidsim tool in the repo.
//  The dumped file that this class deals with usually has the hard-coded name
//  "modalImpulses.txt" from the simulator.
//
// Note:
//  If you wish to read the file "impulses.txt", might want to refer to
//  ImpulseIO class. Its a mess.
//##############################################################################
class ImpulseSeriesReader
{
    public:
        typedef std::shared_ptr<ImpulseSeriesObject> ImpulseSeriesObjectPtr;

    private:
        std::string _impulseFile;
        std::string _rigidsimConfigFile;

        void LoadRigidsimConfig(ImpulseSeriesObjectPtr object);

    public:
        ImpulseSeriesReader(const std::string &impulseFile, const std::string &rigidsimConfigFile)
            : _impulseFile(impulseFile), _rigidsimConfigFile(rigidsimConfigFile)
        {}

        void LoadImpulses(const int &loadObjectID, ImpulseSeriesObjectPtr object);
        void LoadImpulses(const int &loadObjectID, ImpulseSeriesObjectPtr object, std::shared_ptr<FDTD_Objects> objects);
};

#endif
