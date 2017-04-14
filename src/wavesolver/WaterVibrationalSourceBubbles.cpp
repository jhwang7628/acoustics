#include <wavesolver/WaterVibrationalSourceBubbles.h>
#include <wavesolver/FDTD_RigidSoundObject.h>
#include <wavesolver/Wavesolver_ConstantsAndTypes.h>
#include <sndgen/WavReader.hpp>
#include <utils/STL_Wrapper.h>
#include <fstream>

#include "bubbles/FileInput.hpp"

//##############################################################################
//##############################################################################
WaterVibrationalSourceBubbles::
WaterVibrationalSourceBubbles(RigidObjectPtr owner, const std::string &dataDir)
    : VibrationalSource(owner), _surfaceMesh(owner->GetMeshPtr())
{
    Initialize(dataDir);
}

//##############################################################################
//##############################################################################
REAL WaterVibrationalSourceBubbles::
Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time, const int &hintTriangle)
{
    if (normal.dotProduct(_wantedNormal)/normal.length() > _validAngleThreshold)
    {
        // TODO: return acceleration at this point
    }
    else
        return 0.0;
}

//##############################################################################
//##############################################################################
REAL WaterVibrationalSourceBubbles::
Evaluate(const int &vertexID, const Vector3d &vertexNormal, const REAL &time)
{
    throw std::runtime_error("**ERROR** Cannot sample water vibrational source using vertexID");
    return 0.0;
}

//##############################################################################
//##############################################################################
REAL WaterVibrationalSourceBubbles::
EvaluateVelocity(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    throw std::runtime_error("**ERROR** not implemented");
}

//##############################################################################
//##############################################################################
REAL WaterVibrationalSourceBubbles::
EvaluateDisplacement(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    throw std::runtime_error("**ERROR** not implemented");
}

//##############################################################################
//##############################################################################
void WaterVibrationalSourceBubbles::
Initialize(const std::string &dataDir)
{
    std::cout << "Initialize WaterVibrationalSourceBubbles from directory: " << dataDir << std::endl;
    _fileInfo = parseFileNames(dataDir);

    // DEBUGGING
    for (auto iter = _fileInfo.begin(); iter != _fileInfo.end(); ++iter)
    {
        std::cout << iter->first << ":\n"
                  << "   " << iter->second.meshFile << "\n"
                  << "   " << iter->second.datFile << "\n"
                  << "   " << iter->second.freqFile << std::endl;
    }

    exit(1);

    _curTime = -1;
    _t1 = _t2 = -1;
}

//##############################################################################
//##############################################################################
void WaterVibrationalSourceBubbles::
step(REAL time)
{
    if (time >= _t2)
    {
        // Advance
        _t1 = _t2;
        _m1 = _m2;
        _v1 = _v2;
        _kd1 = _kd2;
        _b1 = _b2;

        // New t2
        auto iter = _fileInfo.upper_bound(time);
        if (iter == _fileInfo.end())
        {
            // Past the last solution data time
            _t2 = -1;
        }
        else
        {
            _t2 = iter->first;
            _m2.loadGmsh(iter->second.meshFile);
            _b2 = parseFreqFile(iter->second.freqFile);
            _v2 = loadSurfaceDatFile(_b2,
                                     iter->second.datFile,
                                     _m2);
        }
    }

    // TODO: update all oscillators

    _curTime = time;
}

